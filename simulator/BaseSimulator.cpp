#include "BaseSimulator.h"

#include "./SimulatorParams.h"
#include "./Helper.h"
#include "../file_format/yarnRepr.h"
#include "./threading/threading.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "../easy_profiler_stub.h"

#include <functional>
#include <thread>

namespace simulator {
  BaseSimulator::BaseSimulator(file_format::YarnRepr _yarns,
    SimulatorParams _params):
    thread_pool(std::thread::hardware_concurrency()),
    yarns(_yarns),
    params(_params),
    m(_yarns.yarns[0].points.rows()),
    Q(flatten(_yarns.yarns[0].points)),
    dQ(Eigen::MatrixXd::Zero(3ll * m, 1)),
    F(Eigen::MatrixXd::Zero(3ll * m, 1)),
    constraints(m, &thread_pool)
  {
    log() << "Initializing Simulator" << std::endl;
    log() << "> Found " << m << " control points" << std::endl;

    log() << "Calculating Segment Length" << std::endl;
    segmentLength.resize(m - 1ll);
    for (int i = 0; i < m - 1; i++) {
      segmentLength[i] = params.cInit * (pointAt(Q, i) - pointAt(Q, i + 1)).norm();
    }

    int N = m - 3;
    double totalLength = 0;
    catmullRomLength.resize(N);
    for (int i = 0; i < N; i++) {
      int index = i * 3;
      DECLARE_POINTS2(p, Q, index);

      totalLength += catmullRomLength[i] = params.cInit * integrate<double>([&](double s)->double {
        DECLARE_BASIS_D2(bD, s);
        return POINT_FROM_BASIS(p, bD).norm();
        }, 0, 1);
    }
    log() << "> Total Length: " << totalLength << std::endl;


    log() << "Initializing AABB tree" << std::endl;
    collisionTree = aabb::Tree(3, 0.05, Q.rows() - 4, true);
    std::vector<double> lowerBound;
    std::vector<double> upperBound;
    for (int i = 0; i < m - 3; i++) {
      catmullRomBoundingBox(Q, i, lowerBound, upperBound, yarns.yarns[0].radius);
      collisionTree.insertParticle((unsigned)i, lowerBound, upperBound);
    }
  }

  // Initialize identity mass matrix by default
  void BaseSimulator::constructMassMatrix() {
    M = Eigen::SparseMatrix<double>(3ll * m, 3ll * m);
    M.setIdentity();
    invM = Eigen::SparseMatrix<double>(3ll * m, 3ll * m);
    invM.setIdentity();
  }

  void BaseSimulator::initialize() {
    log() << "Constructing Mass Matrix and Inverse" << std::endl;
    this->constructMassMatrix();
    log() << "Setting-up constraints" << std::endl;
    this->setUpConstraints();
  }

  void BaseSimulator::setPosition(const file_format::YarnRepr& yarn) {
    const Eigen::MatrixXd& v = yarn.yarns[0].points;
    assert(v.cols() == 3 && 3 * v.rows() == dQ.rows());
    Q = flatten(v);
  }

  void BaseSimulator::setVelocity(const file_format::YarnRepr& yarn) {
    const Eigen::MatrixXd& v = yarn.yarns[0].points;
    assert(v.cols() == 3 && 3 * v.rows() == dQ.rows());
    dQ = flatten(v);
  }

  const file_format::YarnRepr& BaseSimulator::getYarns() {
    yarns.yarns[0].points = inflate(Q);
    return yarns;
  }

  file_format::YarnRepr BaseSimulator::getVelocityYarns() {
    file_format::YarnRepr yarn = yarns.createAlike();
    yarn.yarns[0].points = inflate(dQ);
    return yarn;
  }

  ///////////////
  // Stepping
#define WRITE_MATRIX(Q) writeMatrix(#Q"-" + std::to_string(numStep) + ".csv", Q);

  void BaseSimulator::step(const StateGetter& cancelled) {
    EASY_FUNCTION();

    for (int i = 0; i < params.steps; i++) {
      IFDEBUG log() << "Step " << i << std::endl;

      Eigen::MatrixXd originalQ = Q;

      if (cancelled()) break;
      this->stepImpl(cancelled);

      if (cancelled()) break;
      this->fastProjection(cancelled);

      dQ = (Q - originalQ) / params.h;

      if (cancelled()) break;
      this->updateCollisionTree(cancelled);

      if (cancelled()) break;
      this->postStep(cancelled);

      yarns.yarns[0].points = inflate(Q);
    }
  }

  void BaseSimulator::updateCollisionTree(const StateGetter& cancelled) {
    EASY_FUNCTION();

    IFDEBUG log() << "Updating collision tree" << std::endl;
    // Update AABB tree
    std::vector<double> lowerBound;
    std::vector<double> upperBound;
    for (int i = 0; i < m - 3; i++) {
      catmullRomBoundingBox(Q, i, lowerBound, upperBound, yarns.yarns[0].radius);
      collisionTree.updateParticle((unsigned)i, lowerBound, upperBound);
    }
  }

  void BaseSimulator::fastProjection(const StateGetter& cancelled) {
    EASY_FUNCTION();

    int nIter = 0;
    Eigen::MatrixXd constraint;
    Eigen::MatrixXd& Qj = Q;
    double& h = params.h;
    double cValue;
    while ((cValue = maxCoeff(constraint = constraints.calculate(Qj))) > params.fastProjErrorCutoff
      && nIter < params.fastProjMaxIter && !cancelled()) {
      IFDEBUG log() << "- iter: " << nIter << ", constraint: " << cValue << std::endl;

      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

      Eigen::SparseMatrix<double> dConstraint = constraints.getJacobian(Qj);

      EASY_BLOCK("Solving Constraint");

      // Make sure that there's at least one solution
      Eigen::SparseMatrix<double> leftHandSide = dConstraint * invM * dConstraint.transpose();
      for (int i = 0; i < leftHandSide.rows(); i++) {
        leftHandSide.coeffRef(i, i) += 1e-9;
      }

      solver.compute(leftHandSide);
      if (solver.info() != Eigen::Success) {
        log() << "--- solve failed (1)"
          << " iter " << nIter << " STOPPING" << std::endl;
        break;
      }

      Eigen::MatrixXd lambda = solver.solve(constraint);
      if (solver.info() != Eigen::Success) {
        log() << "--- solve failed (2)"
          << " iter " << nIter << " STOPPING" << std::endl;
        break;
      }
      EASY_END_BLOCK;

      Eigen::MatrixXd dQj = invM * dConstraint.transpose() * lambda;
      Qj -= dQj;

      nIter++;
    }
  }

  ///////////////
  // Contact Force

  static inline Eigen::Vector3d contactForce(Eigen::Vector3d direction,
      double distance2, double threshold2) {
    return (-threshold2 / distance2 / distance2 + 1 / threshold2) * direction;
  }

  void BaseSimulator::contactForceBetweenSegments
      (int thread_id,
      std::vector<Eigen::MatrixXd> *forces,
      int ii, int jj) {
    Eigen::MatrixXd &F = (*forces)[thread_id];
    int iIndex = ii * 3;
    int jIndex = jj * 3;

    DECLARE_POINTS2(pi, Q, iIndex);
    DECLARE_POINTS2(pj, Q, jIndex);
    DECLARE_POINTS2(piD, dQ, iIndex);
    DECLARE_POINTS2(pjD, dQ, jIndex);

    double coeffE = params.kContact * catmullRomLength[ii] * catmullRomLength[jj];
    double coeffD = catmullRomLength[ii] * catmullRomLength[jj];
    double kDt = params.kDt;
    double kDn = params.kDn;
    //dataLock.unlock();

    double r = yarns.yarns[0].radius;
    double thresh2 = 4.0 * r * r;

    double step = 1.0 / params.contactForceSamples;
    double halfStep = step / 2;

    Eigen::MatrixXd gradEi = Eigen::MatrixXd::Zero(12, 1);
    Eigen::MatrixXd gradEj = Eigen::MatrixXd::Zero(12, 1);

    Eigen::MatrixXd gradDi = Eigen::MatrixXd::Zero(12, 1);
    Eigen::MatrixXd gradDj = Eigen::MatrixXd::Zero(12, 1);

    for (int i = 0; i < params.contactForceSamples; i++) {
      double si = i * step + halfStep;
      DECLARE_BASIS2(bi, si);
      Eigen::Vector3d Pi = POINT_FROM_BASIS(pi, bi);
      Eigen::Vector3d PiD = POINT_FROM_BASIS(piD, bi);

      for (int j = 0; j < params.contactForceSamples; j++) {
        double sj = j * step + halfStep;
        DECLARE_BASIS2(bj, sj);

        Eigen::Vector3d Pj = POINT_FROM_BASIS(pj, bj);
        Eigen::Vector3d diff = Pj - Pi;

        double norm2 = diff.squaredNorm();
        if (norm2 >= thresh2) continue;

        Eigen::Vector3d PjD = POINT_FROM_BASIS(pjD, bj);
        Eigen::Vector3d diffD = PjD - PiD;

        double tmp = -2 * (kDt - kDn) / sqrt(norm2);

        for (int kk = 0; kk < 4; kk++) {
          // Contact Energy
          gradEi.block<3, 1>(3ll * kk, 0) += contactForce(-2 * bi[kk] * diff, norm2, thresh2);
          gradEj.block<3, 1>(3ll * kk, 0) += contactForce(2 * bj[kk] * diff, norm2, thresh2);

          // Contact Damping
          gradDi.block<3, 1>(3ll * kk, 0) += kDt * (-2 * bi[kk] * diffD) + tmp * (-bi[kk] * diff);
          gradDj.block<3, 1>(3ll * kk, 0) += kDt * (2 * bi[kk] * diffD) + tmp * (bi[kk] * diff);
        }
      }
    }
    F.block<12, 1>(iIndex, 0) -= coeffE * gradEi * step * step;
    F.block<12, 1>(iIndex, 0) -= coeffD * gradDi * step * step;

    F.block<12, 1>(jIndex, 0) -= coeffE * gradEj * step * step;
    F.block<12, 1>(jIndex, 0) -= coeffD * gradDj * step * step;
  }
  
  void BaseSimulator::applyContactForce(const StateGetter& cancelled) {
    EASY_FUNCTION();

    std::vector<Eigen::MatrixXd> forces;
    for (int i = 0; i < thread_pool.size(); i++) {
      forces.push_back(std::move(Eigen::MatrixXd(Q.rows(), Q.cols())));
      forces[i].setZero();
    }

    threading::submitProducerAndWait(thread_pool,
      [this, &forces](int, ctpl::thread_pool *thread_pool){
        int N = m - 3;
        for (int i = 0; i < N; i++) {
          std::vector<unsigned int> intersections = collisionTree.query(i);
          for (int j : intersections) {
            if (j > i + 1) {
              using namespace std::placeholders;
              auto task = std::bind(&BaseSimulator::contactForceBetweenSegments,
                                    this, _1,
                                    &forces,
                                    i, j);
              thread_pool->push(task);
            }
          }
        }
      });

    for (auto force : forces) {
      F += force;
    }
  }

  ///////////////
  // Constraints

  void BaseSimulator::addSegmentLengthConstraint(int i) {
    double length = segmentLength[i];

    Constraints::Func f = [=](const Eigen::MatrixXd& q)->double {
      Eigen::Vector3d p0 = pointAt(q, i);
      Eigen::Vector3d p1 = pointAt(q, i + 1);
      double current = (p1 - p0).norm();
      return current / length - 1;
    };

    Constraints::JacobianFunc fD = [=](const Eigen::MatrixXd& q, Constraints::Referrer ref) {
      Eigen::Vector3d p0 = pointAt(q, i);
      Eigen::Vector3d p1 = pointAt(q, i + 1);
      Eigen::Vector3d diff = p1 - p0;
      double norm = diff.norm();
      for (int ax = 0; ax < 3; ax++) {
        ref(i, ax) = -diff(ax) / length / norm;
        ref(i + 1, ax) = diff(ax) / length / norm;
      }
    };
    constraints.addConstraint(f, fD);
  }

  void BaseSimulator::addCatmullRomLengthConstraint(int i) {
    int index = i * 3;
    double length = catmullRomLength[i];

    Constraints::Func f = [=](const Eigen::MatrixXd& q)->double {
      DECLARE_POINTS2(p, q, index);
      double currentLength = integrate<double>([&](double s)->double {
        DECLARE_BASIS_D2(bD, s);
        return POINT_FROM_BASIS(p, bD).norm();
        }, 0, 1);
      return 1 - currentLength / length;
    };

    using Vec12 = Eigen::Matrix<double, 12, 1>;

    Constraints::JacobianFunc fD = [=](const Eigen::MatrixXd& q, const Constraints::Referrer& ref) {
      DECLARE_POINTS2(p, q, index);
      Vec12 res = integrate<Vec12>([&](double s)->Vec12 {
        Vec12 ans;
        DECLARE_BASIS_D2(bD, s);
        Eigen::Vector3d P = POINT_FROM_BASIS(p, bD);
        double norm = P.norm();

        for (int kk = 0; kk < 4; kk++) {
          ans.block<3, 1>(kk * 3ll, 0) = (bD[kk] / norm) * P;
        }
        return ans;
        }, 0, 1);

      res *= -1.0 / length;

      for (int ii = 0; ii < 12; ii++) {
        ref(i, ii) += res(ii);
      }
    };
    constraints.addConstraint(f, fD);
  }

  void BaseSimulator::addPinConstraint(int i, Eigen::Vector3d point) {
    int index = i * 3;

    for (int ax = 0; ax < 3; ax++) {
      Constraints::Func f = [=](const Eigen::MatrixXd& q)->double {
        return coordAt(q, i, ax) - point(ax);
      };

      Constraints::JacobianFunc fD = [=](const Eigen::MatrixXd& q, const Constraints::Referrer& ref) {
        ref(i, ax) += 1;
      };

      constraints.addConstraint(f, fD);
    }
  }

} // namespace simulator

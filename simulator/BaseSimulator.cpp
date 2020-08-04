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

    initializeContactForceMetaData();
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

      // Solve
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

  void BaseSimulator::initializeContactForceMetaData() {
    catmullRomCoefficient.resize(params.contactForceSamples, 4);

    const double halfStep = 0.5 / params.contactForceSamples;

    Eigen::VectorXd one = Eigen::VectorXd::Ones(params.contactForceSamples);
    // Curve parameter
    Eigen::VectorXd s =
      Eigen::VectorXd::LinSpaced(params.contactForceSamples, 0 + halfStep, 1 - halfStep);
    // s^2
    Eigen::VectorXd s2 = s.cwiseProduct(s);
    // s^3
    Eigen::VectorXd s3 = s2.cwiseProduct(s);

    // Fill in `catmullRomCoefficient`
    catmullRomCoefficient.col(0) = -0.5 * s + s2 - 0.5 * s3;
    catmullRomCoefficient.col(1) = one - 2.5 * s2 + 1.5 * s3;
    catmullRomCoefficient.col(2) = s / 2 + 2 * s2 - 1.5 * s3;
    catmullRomCoefficient.col(3) = -0.5 * s2 + 0.5 * s3;
  }

  void BaseSimulator::contactForceBetweenSegments
      (int thread_id,
      std::vector<Eigen::MatrixXd> *forces,
      int ii, int jj) {
    // EASY_FUNCTION();

    Eigen::MatrixXd &F = (*forces)[thread_id];

    const double step = 1.0 / params.contactForceSamples;

    const double kDt = params.kDt;
    const double kDn = params.kDn;

    // Contact energy gradient
    Eigen::MatrixXd gradEi = Eigen::MatrixXd::Zero(12, 1);
    Eigen::MatrixXd gradEj = Eigen::MatrixXd::Zero(12, 1);

    // Contact force gradient
    Eigen::MatrixXd gradDi = Eigen::MatrixXd::Zero(12, 1);
    Eigen::MatrixXd gradDj = Eigen::MatrixXd::Zero(12, 1);

    // A list of control points position
    Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
      controlPointsI(Q.data() + ii*3);
    Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
      controlPointsJ(Q.data() + jj*3);

    // A list of control points velocity
    Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
      controlPointsVI(dQ.data() + ii*3);
    Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
      controlPointsVJ(dQ.data() + jj*3);

    // A list of sample points position
    Eigen::Matrix3Xd curveI(3, params.contactForceSamples);
    Eigen::Matrix3Xd curveJ(3, params.contactForceSamples);

    // A list of sample points velocity
    Eigen::Matrix3Xd curveVI(3, params.contactForceSamples);
    Eigen::Matrix3Xd curveVJ(3, params.contactForceSamples);

    // Sample the curve
    for (int s = 0; s < params.contactForceSamples; s++) {
      curveI.col(s) = (catmullRomCoefficient.row(s) * controlPointsI).transpose();
      curveJ.col(s) = (catmullRomCoefficient.row(s) * controlPointsJ).transpose();
      curveVI.col(s) = (catmullRomCoefficient.row(s) * controlPointsVI).transpose();
      curveVJ.col(s) = (catmullRomCoefficient.row(s) * controlPointsVJ).transpose();
    }

    // Integrate the force
    const double r = yarns.yarns[0].radius;
    const double thresh2 = 4.0 * r * r; // (2r)^2

    for (int i = 0; i < params.contactForceSamples; i++) {
      auto Pi = curveI.col(i);
      auto Vi = curveVI.col(i);

      for (int j = 0; j < params.contactForceSamples; j++) {
        auto Pj = curveJ.col(j);
        auto Vj = curveVJ.col(j);

        Eigen::Vector3d Pdiff = Pj - Pi;
        double distance = Pdiff.squaredNorm();
        if (distance >= thresh2)
          continue;

        Eigen::Vector3d Vdiff = Vj - Vi;

        double tmp = -2 * (kDt - kDn) / sqrt(distance);
        double coeff = -thresh2 / distance / distance + 1 / thresh2;

        for (int kk = 0; kk < 4; kk++) {
          // Contact Energy
          gradEi.block<3, 1>(3ll * kk, 0)
            += -2 * catmullRomCoefficient(i, kk) * coeff * Pdiff;
          gradEj.block<3, 1>(3ll * kk, 0)
            += 2 * catmullRomCoefficient(j, kk) * coeff * Pdiff;

          // Contact Damping
          gradDi.block<3, 1>(3ll * kk, 0)
            += kDt * (-2) * catmullRomCoefficient(i, kk) * Vdiff
              + tmp * (-catmullRomCoefficient(i, kk)) * Pdiff;
          gradDj.block<3, 1>(3ll * kk, 0)
            += kDt * 2 * catmullRomCoefficient(j, kk) * Vdiff
              + tmp * (catmullRomCoefficient(j, kk) * Pdiff);
        }
      }
    }

    double coeffE = params.kContact * catmullRomLength[ii] * catmullRomLength[jj];
    double coeffD = catmullRomLength[ii] * catmullRomLength[jj];

    F.block<12, 1>(ii * 3, 0) -= coeffE * step * step * gradEi;
    F.block<12, 1>(ii * 3, 0) -= coeffD * step * step * gradDi;

    F.block<12, 1>(jj * 3, 0) -= coeffE * step * step * gradEj;
    F.block<12, 1>(jj * 3, 0) -= coeffD * step * step * gradDj;
  }
  
  void BaseSimulator::applyContactForce(const StateGetter& cancelled) {
    EASY_FUNCTION();

    // Initialize accumulator for each thread
    std::vector<Eigen::MatrixXd> forces;
    for (int i = 0; i < thread_pool.size(); i++) {
      forces.push_back(std::move(Eigen::MatrixXd(Q.rows(), Q.cols())));
      forces[i].setZero();
    }

    // Find all intersecting segments
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

    // Summarize the result of all threads
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

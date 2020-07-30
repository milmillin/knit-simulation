#include "BaseSimulator.h"

#include "./SimulatorParams.h"
#include "./Helper.h"
#include "../file_format/yarnRepr.h"
#include "./threading/threading.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <functional>
#include <thread>

namespace simulator {
  BaseSimulator::BaseSimulator(file_format::YarnRepr _yarns,
    SimulatorParams _params):
    thread_pool(std::thread::hardware_concurrency()),
    yarns(_yarns),
    params(_params),
    m(_yarns.yarns[0].points.rows()),
    numStep(1),
    Q(flatten(_yarns.yarns[0].points)),
    dQ(Eigen::MatrixXf::Zero(3ll * m, 1)),
    F(Eigen::MatrixXf::Zero(3ll * m, 1)),
    constraints(m)
  {
    log() << "Initializing Simulator" << std::endl;
    log() << "> Found " << m << " control points" << std::endl;

    log() << "Calculating Segment Length" << std::endl;
    segmentLength.resize(m - 1ll);
    for (int i = 0; i < m - 1; i++) {
      segmentLength[i] = params.cInit * (pointAt(Q, i) - pointAt(Q, i + 1)).norm();
    }

    int N = m - 3;
    float totalLength = 0;
    catmullRomLength.resize(N);
    for (int i = 0; i < N; i++) {
      int index = i * 3;
      DECLARE_POINTS2(p, Q, index);

      totalLength += catmullRomLength[i] = params.cInit * integrate<float>([&](float s)->float {
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
    M = Eigen::SparseMatrix<float>(3ll * m, 3ll * m);
    M.setIdentity();
    invM = Eigen::SparseMatrix<float>(3ll * m, 3ll * m);
    invM.setIdentity();
  }

  void BaseSimulator::initialize() {
    log() << "Constructing Mass Matrix and Inverse" << std::endl;
    this->constructMassMatrix();
    log() << "Setting-up constraints" << std::endl;
    this->setUpConstraints();
  }

  ///////////////
  // Stepping
#define WRITE_MATRIX(Q) writeMatrix(#Q"-" + std::to_string(numStep) + ".csv", Q);

  void BaseSimulator::step(const StateGetter& cancelled) {
    log() << "Stepping " << params.steps << " step(s) ("
      << numStep << " to " << (numStep + params.steps - 1) << ")" << std::endl;

    for (int i = 0; i < params.steps; i++) {
      IFDEBUG log() << "Step " << numStep << std::endl;
      if (cancelled()) break;
      this->stepImpl(cancelled);
      if (cancelled()) break;
      this->fastProjection(cancelled);
      if (cancelled()) break;
      this->updateCollisionTree(cancelled);
      if (cancelled()) break;
      this->postStep(cancelled);

      IFDEBUG{
        WRITE_MATRIX(Q);
        WRITE_MATRIX(dQ);
        WRITE_MATRIX(F);
      }

      numStep++;
    }
    log() << "> Done " << params.steps << " step(s)" << std::endl;
  }

  void BaseSimulator::updateCollisionTree(const StateGetter& cancelled) {
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
    int nIter = 0;
    Eigen::MatrixXf constraint;
    Eigen::MatrixXf Qj = Q;
    float& h = params.h;
    float cValue;
    while ((cValue = maxCoeff(constraint = constraints.calculate(Qj))) > params.fastProjErrorCutoff
      && nIter < params.fastProjMaxIter && !cancelled()) {
      IFDEBUG log() << "- iter: " << nIter << ", constraint: " << cValue << std::endl;

      Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;

      Eigen::SparseMatrix<float> dConstraint = constraints.getJacobian(Qj);

      solver.compute(dConstraint * invM * dConstraint.transpose());
      if (solver.info() != Eigen::Success) {
        log() << "--- solve failed (1) at step " << numStep 
          << " iter " << nIter << " STOPPING" << std::endl;
        break;
      }

      Eigen::MatrixXf lambda = solver.solve(constraint);
      if (solver.info() != Eigen::Success) {
        log() << "--- solve failed (2) at step " << numStep 
          << " iter " << nIter << " STOPPING" << std::endl;
        break;
      }

      Eigen::MatrixXf dQj = invM * dConstraint.transpose() * lambda;
      Qj -= dQj;

      nIter++;
    }

    dQ += (Qj - Q) / params.h;
    Q = Qj;
  }

  ///////////////
  // Contact Force

  static inline Eigen::Vector3f contactForce(Eigen::Vector3f direction,
      float distance2, float threshold2) {
    return (-threshold2 / distance2 / distance2 + 1 / threshold2) * direction;
  }

  static void contactForceBetweenSegments
      (int thread_id,
      const Eigen::MatrixXf *Q,
      const Eigen::MatrixXf *dQ,
      std::vector<Eigen::MatrixXf> *forces,
      const SimulatorParams *params,
      int ii, int jj,
      float r,
      float length_i, float length_j) {
    Eigen::MatrixXf &F = (*forces)[thread_id];
    int iIndex = ii * 3;
    int jIndex = jj * 3;

    DECLARE_POINTS2(pi, (*Q), iIndex);
    DECLARE_POINTS2(pj, (*Q), jIndex);
    DECLARE_POINTS2(piD, (*dQ), iIndex);
    DECLARE_POINTS2(pjD, (*dQ), jIndex);

    float coeffE = params->kContact * length_i * length_j;
    float coeffD = length_i * length_j;
    float kDt = params->kDt;
    float kDn = params->kDn;
    //dataLock.unlock();

    float thresh2 = 4.0f * r * r;

    float step = 1.f / params->contactForceSamples;
    float halfStep = step / 2;

    Eigen::MatrixXf gradEi = Eigen::MatrixXf::Zero(12, 1);
    Eigen::MatrixXf gradEj = Eigen::MatrixXf::Zero(12, 1);

    Eigen::MatrixXf gradDi = Eigen::MatrixXf::Zero(12, 1);
    Eigen::MatrixXf gradDj = Eigen::MatrixXf::Zero(12, 1);

    for (int i = 0; i < params->contactForceSamples; i++) {
      float si = i * step + halfStep;
      DECLARE_BASIS2(bi, si);
      Eigen::Vector3f Pi = POINT_FROM_BASIS(pi, bi);
      Eigen::Vector3f PiD = POINT_FROM_BASIS(piD, bi);

      for (int j = 0; j < params->contactForceSamples; j++) {
        float sj = j * step + halfStep;
        DECLARE_BASIS2(bj, sj);

        Eigen::Vector3f Pj = POINT_FROM_BASIS(pj, bj);
        Eigen::Vector3f diff = Pj - Pi;

        float norm2 = diff.squaredNorm();
        if (norm2 >= thresh2) continue;

        Eigen::Vector3f PjD = POINT_FROM_BASIS(pjD, bj);
        Eigen::Vector3f diffD = PjD - PiD;

        float tmp = -2 * (kDt - kDn) / sqrt(norm2);

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
    std::vector<Eigen::MatrixXf> forces;
    for (int i = 0; i < thread_pool.size(); i++) {
      forces.push_back(std::move(Eigen::MatrixXf(Q.rows(), Q.cols())));
      forces[i].setZero();
    }

    threading::submitProducerAndWait(thread_pool,
      [this, &forces](int, ctpl::thread_pool *pool){
        int N = m - 3;
        for (int i = 0; i < N; i++) {
          std::vector<unsigned int> intersections = collisionTree.query(i);
          for (int j : intersections) {
            if (j > i + 1) {
              using namespace std::placeholders;
              auto task = std::bind(contactForceBetweenSegments,
                                    _1,
                                    &Q, &dQ, &forces, &params,
                                    i, j, yarns.yarns[0].radius,
                                    catmullRomLength[i], catmullRomLength[j]);
              pool->push(task);
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
    float length = segmentLength[i];

    Constraints::Func f = [=](const Eigen::MatrixXf& q)->float {
      Eigen::Vector3f p0 = pointAt(q, i);
      Eigen::Vector3f p1 = pointAt(q, i + 1);
      float current = (p1 - p0).norm();
      return current / length - 1;
    };

    Constraints::JacobianFunc fD = [=](const Eigen::VectorXf& q, Constraints::Referrer ref) {
      Eigen::Vector3f p0 = pointAt(q, i);
      Eigen::Vector3f p1 = pointAt(q, i + 1);
      Eigen::Vector3f diff = p1 - p0;
      float norm = diff.norm();
      for (int ax = 0; ax < 3; ax++) {
        ref(i, ax) = -diff(ax) / length / norm;
        ref(i + 1, ax) = diff(ax) / length / norm;
      }
    };
    constraints.addConstraint(f, fD);
  }

  void BaseSimulator::addCatmullRomLengthConstraint(int i) {
    int index = i * 3;
    float length = catmullRomLength[i];

    Constraints::Func f = [=](const Eigen::MatrixXf& q)->float {
      DECLARE_POINTS2(p, q, index);
      float currentLength = integrate<float>([&](float s)->float {
        DECLARE_BASIS_D2(bD, s);
        return POINT_FROM_BASIS(p, bD).norm();
        }, 0, 1);
      return 1 - currentLength / length;
    };

    using Vec12 = Eigen::Matrix<float, 12, 1>;

    Constraints::JacobianFunc fD = [=](const Eigen::MatrixXf& q, const Constraints::Referrer& ref) {
      DECLARE_POINTS2(p, q, index);
      Vec12 res = integrate<Vec12>([&](float s)->Vec12 {
        Vec12 ans;
        DECLARE_BASIS_D2(bD, s);
        Eigen::Vector3f P = POINT_FROM_BASIS(p, bD);
        float norm = P.norm();

        for (int kk = 0; kk < 4; kk++) {
          ans.block<3, 1>(kk * 3ll, 0) = (bD[kk] / norm) * P;
        }
        return ans;
        }, 0, 1);

      res *= -1.f / length;

      for (int ii = 0; ii < 12; ii++) {
        ref(i, ii) += res(ii);
      }
    };
    constraints.addConstraint(f, fD);
  }

  void BaseSimulator::addPinConstraint(int i, Eigen::Vector3f point) {
    int index = i * 3;

    for (int ax = 0; ax < 3; ax++) {
      Constraints::Func f = [=](const Eigen::MatrixXf& q)->float {
        return coordAt(q, i, ax) - point(ax);
      };

      Constraints::JacobianFunc fD = [=](const Eigen::MatrixXf& q, const Constraints::Referrer& ref) {
        ref(i, ax) += 1;
      };

      constraints.addConstraint(f, fD);
    }
  }

  const file_format::YarnRepr& BaseSimulator::getYarns() {
    yarns.yarns[0].points = inflate(Q);
    return this->yarns;
  }
} // namespace simulator

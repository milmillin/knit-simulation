#include "BaseSimulator.h"

#include <functional>
#include <thread>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include "easy_profiler_stub.h"
#include "spdlog/spdlog.h"

#include "file_format/yarnRepr.h"
#include "threading/threading.h"
#include "./SimulatorParams.h"
#include "./Helper.h"
#include "./JacobiMethod.h"

namespace simulator {
BaseSimulator::BaseSimulator(file_format::YarnRepr _yarns,
  SimulatorParams _params) :
  thread_pool(std::thread::hardware_concurrency()),
  yarns(reparameterizeYarns(_yarns, _params.reparameterFactor, &segmentLength)),
  params(_params),
  nControlPoints(yarns.vertices.rows()),
  Q(flatten(yarns.vertices)),
  dQ(Eigen::MatrixXd::Zero(3ll * nControlPoints, 1)),
  F(Eigen::MatrixXd::Zero(3ll * nControlPoints, 1)),
  constraints(nControlPoints, &thread_pool),
  contactModelCache(nControlPoints, nControlPoints)
{
  SPDLOG_INFO("Initializing Simulator");
  SPDLOG_INFO("> Found {} yarn(s) with total of {} control points.", yarns.numYarns(), nControlPoints);

  SPDLOG_INFO("Calculating Segment Length");
  segmentLength *= params.cInit;
  double totalLength = 0;
  for (const auto& yarn : yarns.yarns) {
    totalLength += segmentLength * (yarn.size() - 3);
  }
  SPDLOG_INFO("> Total Length: {}", totalLength);


  SPDLOG_INFO("Initializing AABB tree");
  collisionTree = aabb::Tree(3, 0.05, nControlPoints, true);
  std::vector<double> lowerBound;
  std::vector<double> upperBound;
  for (const auto& yarn : yarns.yarns) {
    for (size_t i = yarn.begin; i < yarn.end - 3; i++) {
      catmullRomBoundingBox(Q, i, lowerBound, upperBound, yarn.radius);
      collisionTree.insertParticle((unsigned int)i, lowerBound, upperBound);
    }
  }
}

// Initialize identity mass matrix by default
void BaseSimulator::constructMassMatrix() {
  M = Eigen::SparseMatrix<double>(3ull * nControlPoints, 3ull * nControlPoints);
  M.setIdentity();
  invM = Eigen::SparseMatrix<double>(3ull * nControlPoints, 3ull * nControlPoints);
  invM.setIdentity();
}

void BaseSimulator::initialize() {
  SPDLOG_INFO("Constructing Mass Matrix and Inverse");
  this->constructMassMatrix();
  SPDLOG_INFO("Setting-up constraints");
  this->setUpConstraints();

  initializeContactForceMetaData();
}

void BaseSimulator::setPosition(const file_format::YarnRepr& yarn) {
  const Eigen::MatrixXd& v = yarn.vertices;
  assert(v.cols() == 3 && 3 * v.rows() == dQ.rows());
  Q = flatten(v);
}

void BaseSimulator::setVelocity(const file_format::YarnRepr& yarn) {
  const Eigen::MatrixXd& v = yarn.vertices;
  assert(v.cols() == 3 && 3 * v.rows() == dQ.rows());
  dQ = flatten(v);
}

const file_format::YarnRepr& BaseSimulator::getYarns() {
  yarns.vertices = inflate(Q);
  return yarns;
}

file_format::YarnRepr BaseSimulator::getVelocityYarns() {
  file_format::YarnRepr yarn = yarns.createAlike();
  yarn.vertices = inflate(dQ);
  return yarn;
}

///////////////
// Stepping
#define WRITE_MATRIX(Q) writeMatrix(#Q"-" + std::to_string(numStep) + ".csv", Q);

void BaseSimulator::step(const StateGetter& cancelled) {
  EASY_FUNCTION();
  if (params.statistics) {
    statistics.clear();
  }

  for (int i = 0; i < params.steps; i++) {
    if (params.debug) {
      SPDLOG_INFO("Step {}", i);
    }

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

  }

  if (params.statistics) {
    printStatistics();
  }
}

void BaseSimulator::updateCollisionTree(const StateGetter& cancelled) {
  EIGEN_UNUSED_VARIABLE(cancelled)
  EASY_FUNCTION();

  if (params.debug)
    SPDLOG_INFO("Updating collision tree");

  // Update AABB tree
  std::vector<double> lowerBound;
  std::vector<double> upperBound;
  for (const auto& yarn : yarns.yarns) {
    for (size_t i = yarn.begin; i < yarn.end - 3; i++) {
      catmullRomBoundingBox(Q, i, lowerBound, upperBound, yarn.radius);
      collisionTree.updateParticle((unsigned)i, lowerBound, upperBound);
    }
  }
}

void BaseSimulator::fastProjection(const StateGetter& cancelled) {
  EASY_FUNCTION();

  int nIter = 0;
  Eigen::MatrixXd constraint;
  Eigen::MatrixXd& Qj = Q;
  double cValue;
  while ((cValue = maxCoeff(constraint = constraints.calculate(Qj))) > params.fastProjErrorCutoff
    && nIter < params.fastProjMaxIter && !cancelled()) {
    if (params.debug)
      SPDLOG_INFO("- iter: {}, constraint: {}", nIter, cValue);

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
      SPDLOG_ERROR("--- solve failed (1) at iteration {}, skippking this step", nIter);
      break;
    }

    Eigen::MatrixXd lambda = solver.solve(constraint);
    if (solver.info() != Eigen::Success) {
      SPDLOG_ERROR("--- solve failed (2) at iteration {}, skippking this step", nIter);
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

void BaseSimulator::contactForce(size_t i, size_t j, ControlPoints* forceI, ControlPoints* forceJ) {
  // EASY_FUNCTION();
  forceI->setZero();
  forceJ->setZero();

  const double step = 1.0 / params.contactForceSamples;

  // A list of control points position
  Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
    controlPointsI(Q.data() + i * 3);
  Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
    controlPointsJ(Q.data() + j * 3);

  // A list of sample points position
  Eigen::Matrix3Xd curveI(3, params.contactForceSamples);
  Eigen::Matrix3Xd curveJ(3, params.contactForceSamples);

  // Sample the curve
  for (int s = 0; s < params.contactForceSamples; s++) {
    curveI.col(s) = (catmullRomCoefficient.row(s) * controlPointsI).transpose();
    curveJ.col(s) = (catmullRomCoefficient.row(s) * controlPointsJ).transpose();
  }

  // Integrate the force
  const double r = yarns.yarns[0].radius;
  const double thresh2 = 4.0 * r * r; // (2r)^2

  for (int i = 0; i < params.contactForceSamples; i++) {
    auto Pi = curveI.col(i);

    for (int j = 0; j < params.contactForceSamples; j++) {
      auto Pj = curveJ.col(j);

      Eigen::Vector3d Pdiff = Pj - Pi;
      double distance = Pdiff.squaredNorm();
      if (distance >= thresh2)
        continue;

      double coeff = -thresh2 / distance / distance + 1 / thresh2;

      for (int kk = 0; kk < 4; kk++) {
        // Contact Energy
        forceI->block<3, 1>(3ll * kk, 0)
          += -2 * catmullRomCoefficient(i, kk) * coeff * Pdiff;
        forceJ->block<3, 1>(3ll * kk, 0)
          += 2 * catmullRomCoefficient(j, kk) * coeff * Pdiff;
      }
    }
  }

  double coeffE = params.kContact * segmentLength * segmentLength;

  (*forceI) *= -coeffE * step * step;
  (*forceJ) *= -coeffE * step * step;
}

void BaseSimulator::buildLinearModel(size_t i, size_t j, LinearizedContactModel* model) {
  // EASY_FUNCTION();
  const double step = 1.0 / params.contactForceSamples;

  // Save location
  model->baseQI = Q.block<12, 1>(i * 3ull, 0);
  model->baseQJ = Q.block<12, 1>(j * 3ull, 0);

  // Save offset
  Eigen::Vector3d center = Eigen::Vector3d::Zero();  // Center of control points
  for (int k = 0; k < 4; k++) {
    center += pointAt(model->baseQI, k);
    center += pointAt(model->baseQJ, k);
  }
  center /= 8;

  for (int k = 0; k < 4; k++) {
    pointAt(model->RelativeQI, k) = pointAt(model->baseQI, k) - center;
    pointAt(model->RelativeQJ, k) = pointAt(model->baseQJ, k) - center;
  }

  // Save contact force
  contactForce(i, j, &model->baseForceI, &model->baseForceJ);

  // Sample the curve
  Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
    controlPointsI(Q.data() + i * 3);
  Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>
    controlPointsJ(Q.data() + j * 3);

  Eigen::Matrix3Xd curveI(3, params.contactForceSamples);
  Eigen::Matrix3Xd curveJ(3, params.contactForceSamples);
  for (int s = 0; s < params.contactForceSamples; s++) {
    curveI.col(s) = (catmullRomCoefficient.row(s) * controlPointsI).transpose();
    curveJ.col(s) = (catmullRomCoefficient.row(s) * controlPointsJ).transpose();
  }

  // Integrate the jecobian using samples
  model->dForceII.setZero();
  model->dForceIJ.setZero();
  model->dForceJI.setZero();
  model->dForceJJ.setZero();

  const double r = yarns.yarns[0].radius;
  const double thresh2 = 4.0 * r * r; // (2r)^2

  for (int i = 0; i < params.contactForceSamples; i++) {
    auto Pi = curveI.col(i);

    for (int j = 0; j < params.contactForceSamples; j++) {
      auto Pj = curveJ.col(j);

      Eigen::Vector3d Pdiff = Pj - Pi;
      double distance2 = Pdiff.squaredNorm();
      if (distance2 >= thresh2)
        continue;

      double coeff = -thresh2 / distance2 / distance2 + 1 / thresh2;

      Eigen::Matrix3d term =
        2.0 * (16.0 * r * r / distance2 / distance2 / distance2
          * Pdiff * Pdiff.transpose()
          + coeff * Eigen::Matrix3d::Identity());

      for (size_t kf = 0; kf < 4; kf++) {
        for (size_t kx = 0; kx < 4; kx++) {
          model->dForceII.block<3, 3>(kf * 3ull, kx * 3ull)
            += catmullRomCoefficient(i, kf) * catmullRomCoefficient(i, kx) * term;
          model->dForceIJ.block<3, 3>(kf * 3ull, kx * 3ull)
            += -catmullRomCoefficient(i, kf) * catmullRomCoefficient(j, kx) * term;

          model->dForceJI.block<3, 3>(kf * 3ull, kx * 3ull)
            += -catmullRomCoefficient(j, kf) * catmullRomCoefficient(i, kx) * term;
          model->dForceJJ.block<3, 3>(kf * 3ull, kx * 3ull)
            += catmullRomCoefficient(j, kf) * catmullRomCoefficient(j, kx) * term;
        }
      }
    }
  }

  double coefficient = -params.kContact * segmentLength * segmentLength * step * step;

  (model->dForceII) *= coefficient;
  (model->dForceIJ) *= coefficient;
  (model->dForceJI) *= coefficient;
  (model->dForceJJ) *= coefficient;

  // Initialize jacobi method metadata
  model->jacobiV = Eigen::Matrix3d::Identity();

  // Mark as valid
  model->lastUpdate = 0;
  model->valid = true;

  // Update statistics
  if (params.statistics) {
    statistics.linearizedModelRebuildCount++;
  }
}

bool BaseSimulator::applyApproxContactForce(size_t i, size_t j,
  Eigen::MatrixXd& forces, LinearizedContactModel* model) {
  // EASY_FUNCTION();

  auto QI = Q.block<12, 1>(i * 3ull, 0);
  auto QJ = Q.block<12, 1>(j * 3ull, 0);
  ControlPoints dQI = QI - model->baseQI;
  ControlPoints dQJ = QJ - model->baseQJ;

  // Apply offset to the model to adjust for moving
  Eigen::Vector3d positionOffset = Eigen::Vector3d::Zero();
  for (int k = 0; k < 4; k++) {
    positionOffset += pointAt(dQI, k);
    positionOffset += pointAt(dQJ, k);
  }
  positionOffset /= 8.0;

  checkNaN(positionOffset);

  for (int k = 0; k < 4; k++) {
    pointAt(dQI, k) -= positionOffset;
    pointAt(dQJ, k) -= positionOffset;
  }

  // Apply rotation to the model to adjust for rotation
  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for (int k = 0; k < 4; k++) {
    center += pointAt(QI, k);
    center += pointAt(QJ, k);
  }
  center /= 8.0;
  checkNaN(center);

  Eigen::Matrix3d Apq = Eigen::Matrix3d::Zero();
  for (int k = 0; k < 4; k++) {
    Apq += (pointAt(QI, k) - center) * pointAt(model->RelativeQI, k).transpose();
  }
  for (int k = 0; k < 4; k++) {
    Apq += (pointAt(QJ, k) - center) * pointAt(model->RelativeQJ, k).transpose();
  }

  Eigen::Matrix3d A = Apq.transpose() * Apq;

  // Store the ground truth calculated by Eigen
  // Compare it with Jacobi method later
  Eigen::Matrix3d rotationGroundTruth = Eigen::Matrix3d::Zero();
  if (params.statistics) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(A);
    if (solver.info() != Eigen::Success) {
    }
    else {
      rotationGroundTruth = Apq * solver.operatorInverseSqrt();
    }
  }

  // Warm start Jacobi method
  A = model->jacobiV.transpose() * A * model->jacobiV;

  // Estimate rotation using Jacobi method
  cyclicJacobi(&A, &(model->jacobiV), 1);
  inverseSquareRoot(&A, model->jacobiV);
  Eigen::Matrix3d rotation = Apq * A;

  // The inverse of a rotational matrix is the the transpose
  Eigen::Matrix3d inverseRotation = rotation.transpose();
  checkNaN(inverseRotation);

  // Rotate the points
  for (int k = 0; k < 4; k++) {
    pointAt(dQI, k) = inverseRotation * pointAt(dQI, k);
    pointAt(dQJ, k) = inverseRotation * pointAt(dQJ, k);
  }

  // Check the error by Jacobi method
  if (params.statistics) {
    double error = (rotation - rotationGroundTruth).cwiseAbs().sum() / rotationGroundTruth.cwiseAbs().sum();
    double oldError = statistics.jacobiMethodError;
    while (!statistics.jacobiMethodError.compare_exchange_strong(
      oldError, error + oldError)) {
      oldError = statistics.jacobiMethodError;
    }
    statistics.jacobiMethodErrorCount++;
  }

  // EASY_END_BLOCK; EASY_BLOCK("Deviation");

  // === Check if a model rebuild is needed ===
  double deviation = std::max(dQI.cwiseAbs().maxCoeff(), dQJ.cwiseAbs().maxCoeff());
  const double r = yarns.yarns[0].radius;
  if (deviation >= params.contactModelTolerance * r) {
    // Need a model rebuild
    return false;
  }

  // EASY_END_BLOCK; EASY_BLOCK("Estimation");

  // === Estimate the force ===
  ControlPoints forceI = model->baseForceI + model->dForceII * dQI + model->dForceIJ * dQJ;
  ControlPoints forceJ = model->baseForceJ + model->dForceJI * dQI + model->dForceJJ * dQJ;

  // EASY_END_BLOCK; EASY_BLOCK("Revert Rotation");

  // === Revert the rotation ===
  for (int k = 0; k < 4; k++) {
    pointAt(forceI, k) = rotation * pointAt(forceI, k);
    pointAt(forceJ, k) = rotation * pointAt(forceJ, k);
  }

  // EASY_END_BLOCK; EASY_BLOCK("Apply");

  // === Apply force ===
  forces.block<12, 1>(i * 3, 0) += forceI;
  forces.block<12, 1>(j * 3, 0) += forceJ;

  // EASY_END_BLOCK; EASY_BLOCK("Statistics");

  // === Calculate statistics ===
  if (params.statistics) {
    ControlPoints realForceI;
    ControlPoints realForceJ;
    contactForce(i, j, &realForceI, &realForceJ);

    double error = (realForceI - forceI).cwiseAbs().sum()
      + (realForceI - forceI).cwiseAbs().sum();

    statistics.contactForceErrorDivider += 1;
    double oldError = statistics.contactForceTotalError;
    while (!statistics.contactForceTotalError
      .compare_exchange_weak(oldError, oldError + error)) {
      oldError = statistics.contactForceTotalError;
    }
  }
  model->lastUpdate++;

  // Estimation applied
  return true;
}

void BaseSimulator::applyContactForceBetweenSegments
(int thread_id,
  std::vector<Eigen::MatrixXd>* forces,
  size_t ii, size_t jj) {
  EASY_FUNCTION();

  // Update statistics
  if (params.statistics) {
    statistics.totalContactCount++;
  }

  // Get linearized model
  // We swap `ii` and `jj` so that there's less lock contention
  auto& model = contactModelCache.lock(jj, ii);

  // Try to use approximation
  if (model.valid && model.lastUpdate < params.maxContactModelUpdateInterval
    && applyApproxContactForce(ii, jj, (*forces)[thread_id], &model)) {
    statistics.approximationUsedCount++;
  }
  else {
    // Need to rebuild the linearized model
    buildLinearModel(ii, jj, &model);
    applyApproxContactForce(ii, jj, (*forces)[thread_id], &model);
  }

  // Release the lock
  contactModelCache.unlock(jj, ii);
}

void BaseSimulator::applyContactForce(const StateGetter& cancelled) {
  EIGEN_UNUSED_VARIABLE(cancelled)
  EASY_FUNCTION();

  // Initialize accumulator for each thread
  std::vector<Eigen::MatrixXd> forces;
  for (size_t i = 0; i < thread_pool.size(); i++) {
    forces.push_back(Eigen::MatrixXd::Zero(Q.rows(), Q.cols()));
  }

  // Find all intersecting segments
  threading::submitProducerAndWait(thread_pool,
    [this, &forces](int, ctpl::thread_pool* thread_pool) {

      for (const auto& yarn : yarns.yarns) {
        for (size_t i = yarn.begin; i < yarn.end - 3; i++) {
          for (auto j : collisionTree.query(i)) {
            if (j > i + 1) {
              using namespace std::placeholders;
              auto task = std::bind(&BaseSimulator::applyContactForceBetweenSegments,
                this, _1,
                &forces,
                i, (size_t)j);
              thread_pool->push(task);
            }
          }
        }
      }
    });

  // Summarize the result of all threads
  for (auto force : forces) {
    F += force;
    checkNaN(force);
    checkNaN(F);
  }
}

///////////////////////
// Length spring
void BaseSimulator::applyLengthSpringForce() {
  for (const auto& yarn : yarns.yarns) {
    for (size_t i = yarn.begin; i < yarn.end - 1; i++) {
      Eigen::Vector3d force = pointAt(Q, i + 1) - pointAt(Q, i);
      double distance = force.norm();
      force *= params.kLen * (distance - segmentLength) / distance;
      F.block<3, 1>(i * 3, 0) += force;
      F.block<3, 1>((i + 1) * 3, 0) -= force;
    }
  }
}

///////////////
// Constraints

void BaseSimulator::addSegmentLengthConstraint(size_t i) {
  double length = segmentLength;

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

void BaseSimulator::addCatmullRomLengthConstraint(size_t i) {
  int index = i * 3;
  double length = segmentLength;

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

void BaseSimulator::addPinConstraint(size_t i, Eigen::Vector3d point) {
  for (int ax = 0; ax < 3; ax++) {
    Constraints::Func f = [=](const Eigen::MatrixXd& q)->double {
      return coordAt(q, i, ax) - point(ax);
    };

    Constraints::JacobianFunc fD = [=](const Eigen::MatrixXd& q, const Constraints::Referrer& ref) {
      EIGEN_UNUSED_VARIABLE(q)

      ref(i, ax) += 1;
    };

    constraints.addConstraint(f, fD);
  }
}

} // namespace simulator

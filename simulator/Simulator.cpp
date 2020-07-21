#include "Simulator.h"

#include "./SimulatorParams.h"
#include "./Helper.h"
#include "../file_format/yarnRepr.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <algorithm>

namespace simulator {


//////////////////////////////////////////////
//
// Segment Length
//

void Simulator::calculateSegmentLength() {
	log() << "Calculating Segment Length" << std::endl;
	int N = m - 3;
	segmentLength = std::vector<float>(N);
	
	// Compress yarns
	float coeff = params.cInit;

	float totalLength = 0;
	for (int i = 0; i < N; i++) {
		// Calculate arclength of segment[i].
		// Note that the length is fixed throughout the simulation.

		int index = i * 3;
		DECLARE_POINTS2(p, q, index);

		totalLength += segmentLength[i] = coeff * integrate<float>([&](float s)->float {
			DECLARE_BASIS_D2(bD, s);
			return POINT_FROM_BASIS(p, bD).norm();
			}, 0, 1);
	}
	log() << "Total Length: " << totalLength << std::endl;
}

void Simulator::addSegmentLengthConstraint() {
	int N = m - 3;

	for (int i = 0; i < N; i++) {
		constraints.addLengthConstrain(i, segmentLength[i]);
	}
}

//////////////////////////////////////////////
//
// Mass Matrix
//


// Mass contribution matrix of Catmull-Rom curve with tightness 0.5.
// b0(s) = (1/2)(-s + 2s^2 - s^3)
// b1(s) = (1/2)(2 - 5s^2 + 3s^3)
// b2(s) = (1/2)(s + 4s^2 - 3s^3)
// b3(s) = (1/2)(-s^2 + s^3)
//
// MASS_CONTRIBUTION[i][j] = integrate(bi(s)*bj(s), (s, 0, 1))
// TODO: check correctness
static const float MASS_CONTRIBUTION[4][4] =
{
	{1.f / 420, -47.f / 1680, -1.f / 56, 1.f / 560},
	{-47.f / 1680, 17.f / 42, 307.f / 1680, -1.f / 56},
	{-1.f / 56, 307.f / 1680, 17.f / 42, -47.f / 1680},
	{1.f / 560, -1.f / 56, -47.f / 1680, 1.f / 420}
};

void Simulator::constructMassMatrix() {
	log() << "Constructing Mass Matrix" << std::endl;
	M = Eigen::SparseMatrix<float>(3 * m, 3 * m);

	// Iterate though each segment
	int N = m - 3;

	float mUnit = params.m;
	for (int i = 0; i < N; i++) {
		// Iterate through combination of control points to find mass contribution.
		for (int j = 0; j < 4; j++) {
			for (int k = j; k < 4; k++) {
				// M[(i+j)*3][(i+k)*3] = mUnit * l[i] * integrate(bj(s), bk(s), (s, 0, 1))
				float contribution = mUnit * segmentLength[i] * MASS_CONTRIBUTION[j][k];
				M.coeffRef((i + j) * 3, (i + k) * 3) += contribution;
				M.coeffRef((i + j) * 3 + 1, (i + k) * 3 + 1) += contribution;
				M.coeffRef((i + j) * 3 + 2, (i + k) * 3 + 2) += contribution;

				// transpose
				if (j != k) {
          M.coeffRef((i + k) * 3, (i + j) * 3) += contribution;
          M.coeffRef((i + k) * 3 + 1, (i + j) * 3 + 1) += contribution;
          M.coeffRef((i + k) * 3 + 2, (i + j) * 3 + 2) += contribution;
				}
			}
		}
	}

	log() << "Calculating Mass Matrix Inverse" << std::endl;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
	solver.compute(M);
	if (solver.info() != Eigen::Success) {
		log() << "Decomposition Failed" << std::endl;
		return;
	}

	Eigen::SparseMatrix<float> I(M.rows(), M.cols());
	I.setIdentity();

	MInverse = solver.solve(I);
	if (solver.info() != Eigen::Success) {
		log() << "Solve Failed" << std::endl;
		return;
	}
}

//////////////////////////////////////////////
//
// Energy Gradient
//

void Simulator::calculateGradient(const std::function<bool()> cancelled) {
	log() << "Calculating Gradient" << std::endl;
	int N = m - 3;

	// Energy
	log() << "- Length Energy" << std::endl;
	for (int i = 0; i < N; i++) {
		calculateLengthEnergyGradient(i);
	}

	log() << "- Bending Energy" << std::endl;
	for (int i = 0; i < N; i++) {
		calculateBendingEnergyGradient(i);
	}

	log() << "- Collision Energy" << std::endl;
	{
		ParallelWorker worker(cancelled);

    for (int i = 0; i < N; i++) {
			worker.addWork([this, i, N]() {
				for (int j = i + 2; j < N; j++) {
					calculateCollisionGradient(i, j);
				}
				});
    }

		worker.run();
	}

  // Damping
	log() << "- Global Damping" << std::endl;
	for (int i = 0; i < N; i++) {
		calculateGlobalDampingGradient(i);
	}
}


//////////////////////////////////////////////
//
// Fast Projection
//

float maxCoeff(const Eigen::MatrixXf& m) {
	if (m.rows() == 0 || m.cols() == 0) return 0;
	return std::max(fabs(m.maxCoeff()), fabs(m.minCoeff()));
}

void Simulator::fastProjection() {
	log() << "Fast Projection" << std::endl;
	float h = params.h;			// timestep
	static const float eps = 1e-5;

	// TODO: check sign of f
	// Unconstrained step
	Eigen::MatrixXf qj = q + h * qD - (h * h) * (MInverse * (gradE + gradD - f));

	// Projection
	int iter = 0;
	Eigen::MatrixXf constraint;
	float cValue;
	while ((cValue = maxCoeff(constraint = constraints.calculate(qj))) > eps && iter < 100) {
		log() << "- iter: " << iter << ", constraint: " << cValue << std::endl;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;

		Eigen::SparseMatrix<float> jC = constraints.getJacobian(qj);

		//writeMatrix("jC-" + std::to_string(stepCnt) + "-" + std::to_string(cValue) + ".csv", Eigen::MatrixXf(jC));

		solver.compute((h * h) * jC * MInverse * jC.transpose());
		if (solver.info() != Eigen::Success) {
			log() << "--- solve failed. STOPPING" << std::endl;
			break;
		}

		Eigen::MatrixXf lambda = solver.solve(constraint);
		if (solver.info() != Eigen::Success) {
			log() << "--- solve failed. STOPPING" << std::endl;
			break;
		}

		Eigen::MatrixXf qjD = (-h * h) * MInverse * jC.transpose() * lambda;
		qj += qjD;

		iter++;
	}
  log() << "- iter: " << iter << ", constraint: " << cValue << std::endl;

	// Update velocity and position;
	qD = (qj - q) / h;
	q = qj;
}

//////////////////////////////////////////////
//
// Simulator implementation
//

Simulator::Simulator(file_format::YarnRepr yarns, SimulatorParams params_) : constraints(0), stepCnt(0) {
	log() << "Initializing Simulator" << std::endl;
	this->yarns = yarns;
	this->params = params_;
	if (yarns.yarns.size() > 0) {
		q = yarns.yarns[0].points;
	}
	m = q.rows();

	log() << "Found " << m << " control points" << std::endl;

	// Initialize constraints
	constraints = Constraints(m);

	// Resize to vector 3m x 1
	q = flatten(q);

	// Initialize to zero
	qD = Eigen::MatrixXf::Zero(3 * m, 1);
	gradE = Eigen::MatrixXf::Zero(3 * m, 1);
	gradD = Eigen::MatrixXf::Zero(3 * m, 1);
	f = Eigen::MatrixXf::Zero(3 * m, 1);

	// Calculate Segment Length;
	calculateSegmentLength();
	addSegmentLengthConstraint();

	// Construct Mass Matrix
	constructMassMatrix();

	// Pin first and last control points
	constraints.addPinConstrain(0, q.block<3, 1>(0, 0));
	constraints.addPinConstrain(m - 1, q.block<3, 1>(3 * (m - 1), 0));
	//constraints.addGlueConstrain(0, 1);
	//constraints.addGlueConstrain(m - 2, m - 1);

	writeToFile();

	log() << "Simulator Initialized" << std::endl;
}

Simulator::~Simulator() {
	log() << "Simulator Destroyed" << std::endl;
}

void Simulator::writeToFile() const {
	writeMatrix("q-" + std::to_string(stepCnt) + ".csv", q);
	writeMatrix("qD-" + std::to_string(stepCnt) + ".csv", qD);
	writeMatrix("gradE-" + std::to_string(stepCnt) + ".csv", gradE);
}

void Simulator::step(const std::function<bool()>& cancelled) {
	stepCnt++;
	log() << "Step (" << stepCnt << ")" << std::endl;

	gradE.setZero();
	gradD.setZero();

	// update gradE, gradD, f
	if (!cancelled()) {
    calculateGradient(cancelled);
	}

	if (!cancelled()) {
    fastProjection();
	}

	// save to YarnRepr
	// supports one yarn for now
	yarns.yarns[0].points = inflate(q);

	writeToFile();

	log() << "Done Step (" << stepCnt << ")" << std::endl;
}

};  // namespace simulator

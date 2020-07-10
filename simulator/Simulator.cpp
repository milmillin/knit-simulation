#include "Simulator.h"

#include "./SimulatorParams.h"
#include "./Helper.h"
#include "../file_format/yarnRepr.h"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <fstream>

namespace simulator {


//////////////////////////////////////////////
//
// Helper Functions
//

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
void writeMatrix(std::string filename, const Eigen::MatrixXf& q) {
	std::ofstream f;
	f.open(filename);
	f << q.format(CSVFormat) << "\n";
	f.close();
}

// Reshapes #m x 3 matrix into a vector of (3 * #m) rows.
// Returns a new matrix.
Eigen::MatrixXf flatten(const Eigen::MatrixXf& v) {
	Eigen::MatrixXf res(v.rows() * v.cols(), 1);
	int r = v.rows();
	int c = v.cols();
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			res(i * c + j, 0) = v(i, j);
		}
	}
	return res;
}

// Reshapes a vector of (3 * #m) rows to a #m x 3 matrix.
// Returns a new matrix.
Eigen::MatrixXf inflate(const Eigen::MatrixXf& v, size_t col) {
	assert(v.cols() == 1);
	assert(v.rows() % col == 0);

	int r = v.rows() / col;

	Eigen::MatrixXf res(r, col);
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < col; j++) {
			res(i, j) = v(i * col + j, 0);
		}
	}
	return res;
}

//////////////////////////////////////////////
//
// Segment Length
//

void Simulator::calculateSegmentLength() {
	int N = m - 3;
	segmentLength = std::vector<float>(N);

	float totalLength = 0;
	for (int i = 0; i < N; i++) {
		// Calculate arclength of segment[i].
		// Note that the length is fixed throughout the simulation.
		segmentLength[i] = catmullRomArcLength(q, i * 3);
		totalLength += segmentLength[i];
	}

	std::cout << "Total Length: " << totalLength << std::endl;
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
	std::cout << "Constructing Mass Matrix" << std::endl;
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
			}
		}
	}

	std::cout << "Calculating Mass Matrix Inverse" << std::endl;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> solver;
	solver.compute(M);
	if (solver.info() != Eigen::Success) {
		std::cout << "Decomposition Failed" << std::endl;
		return;
	}

	Eigen::SparseMatrix<float> I(M.rows(), M.cols());
	I.setIdentity();

	MInverse = solver.solve(I);
	if (solver.info() != Eigen::Success) {
		std::cout << "Solve Failed" << std::endl;
		return;
	}
}

//////////////////////////////////////////////
//
// Energy Gradient
//

void Simulator::calculateGradient() {
	int N = m - 3;
	for (int i = 0; i < N; i++) {
		// Energy
		//calculateBendingEnergyGradient(i);
		calculateLengthEnergyGradient(i);

		// Damping
		//calculateGlobalDampingGradient(i);
	}
}


//////////////////////////////////////////////
//
// Fast Projection
//

void Simulator::fastProjection() {
	float h = params.h;			// timestep
	static const float eps = 1e-8;

	// TODO: check sign of f

	// Unconstrained step
	Eigen::MatrixXf qj = q + h * qD - (h * h) * (MInverse * (gradE + gradD - f));

	// Projection
	while (constraints.calculateMax(qj) > eps) {
		Eigen::MatrixXf jC = constraints.getJacobian(qj);
		Eigen::MatrixXf tmp = jC * MInverse * jC.transpose();

		Eigen::MatrixXf lambda = tmp.inverse() * (constraints.calculate(qj) / (h * h));
		Eigen::MatrixXf qjD = (-h * h) * MInverse * jC.transpose() * lambda;
		qj += qjD;
	}

	// Update velocity and position;
	qD = (qj - q) / h;
	q = qj;
}

//////////////////////////////////////////////
//
// Simulator implementation
//

Simulator::Simulator(file_format::YarnRepr yarns, SimulatorParams params_) : params(params_), constraints(0), stepCnt(0) {
	this->yarns = yarns;
	if (yarns.yarns.size() > 0) {
		q = yarns.yarns[0].points;
	}
	m = q.rows();

	// Initialize constraints
	constraints = Constraints(m);

	// Resize to vector 3m x 1
	q = flatten(q);
	// q.resize(3 * m, 1); won't work because Column-major


	// Initialize to zero
	qD = Eigen::MatrixXf::Zero(3 * m, 1);
	gradE = Eigen::MatrixXf::Zero(3 * m, 1);
	gradD = Eigen::MatrixXf::Zero(3 * m, 1);
	f = Eigen::MatrixXf::Zero(3 * m, 1);

	// Calculate Segment Length;
	calculateSegmentLength();

	// Construct Mass Matrix
	constructMassMatrix();

	writeMatrix("mass.csv", Eigen::MatrixXf(M));
	writeMatrix("massInverse.csv", Eigen::MatrixXf(MInverse));

	writeToFile();

	std::cout << "Simulator Initialized" << std::endl;
}


void Simulator::writeToFile() const {
	writeMatrix("q-" + std::to_string(stepCnt) + ".csv", q);
	writeMatrix("qD-" + std::to_string(stepCnt) + ".csv", qD);
	writeMatrix("gradE-" + std::to_string(stepCnt) + ".csv", gradE);
}

void Simulator::step() {
	std::cout << "Step" << std::endl;

	gradE.setZero();
	gradD.setZero();

	// update gradE, gradD, f
	calculateGradient();

	fastProjection();

	// save to YarnRepr
	// supports one yarn for now
	yarns.yarns[0].points = inflate(q, 3);

	stepCnt++;
	writeToFile();

	std::cout << ">> Done Step" << std::endl;
}

};  // namespace simulator

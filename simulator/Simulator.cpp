#include "Simulator.h"

#include "SimulatorParams.h"
#include "Helper.h"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace simulator {


//////////////////////////////////////////////
//
// Helper Functions
//

// Reshapes #m x 3 matrix into a vector of (3 * #m) rows.
// Returns a new matrix.
Eigen::MatrixXf flatten(const Eigen::MatrixXf &v) {
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
	
	float totalLength = 0;
	float mUnit = params.m;
	for (int i = 0; i < N; i++) {
		// Calculate arclength of segment[i].
		// Note that the length is fixed throughout the simulation.
		float length = catmullRomArcLength(q, i * 3);
		totalLength += length;

		// Iterate through combination of control points to find mass contribution.
		for (int j = 0; j < 4; j++) {
			for (int k = j; k < 4; k++) {
				// M[(i+j)*3][(i+k)*3] = mUnit * l[i] * integrate(bj(s), bk(s), (s, 0, 1))
				float contribution = mUnit * length * MASS_CONTRIBUTION[j][k];
				M.coeffRef((i + j) * 3, (i + k) * 3) += contribution;
				M.coeffRef((i + j) * 3 + 1, (i + k) * 3 + 1) += contribution;
				M.coeffRef((i + j) * 3 + 2, (i + k) * 3 + 2) += contribution;
			}
		}
	}

	std::cout << "Total Length: " << totalLength << std::endl;

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
// Fast Projection
//

float Simulator::constraint(const Eigen::MatrixXf& q)
{
	return 0.0f;
}

void Simulator::fastProjection() {
	int j = 0;							// current projection iteration
	float h = params.h;			// timestep
	static const float eps = 1e-8;

	// TODO: check sign of f
	Eigen::MatrixXf qj = q + h * qD - (h * h) * (MInverse * (gradE + gradD - f));

	while (constraint(qj) > eps) {
		// TODO:
	}

	// Update velocity and position;
	qD = (qj - q) / h;
	q = qj;
}

//////////////////////////////////////////////
//
// Simulator implementation
//

Simulator::Simulator(Eigen::MatrixXf q_, SimulatorParams params_) : params(params_) {
	assert(q_.cols() == 3);

	m = q_.rows();

	// Resize to vector 3m x 1
	q = flatten(q_);
	// q.resize(3 * m, 1); won't work because Column-major


	// Initialize to zero
	qD = Eigen::MatrixXf::Zero(3 * m, 1);
	gradE = Eigen::MatrixXf::Zero(3 * m, 1);
	gradD = Eigen::MatrixXf::Zero(3 * m, 1);
	f = Eigen::MatrixXf::Zero(3 * m, 1);

	// Construct Mass Matrix
	constructMassMatrix();

	std::cout << "Simulator Initialized" << std::endl;
}

Eigen::MatrixXf Simulator::getControlPoints() const {
	return inflate(q, 3);
}

void Simulator::step() {
	//TODO:
	// update gradE, gradD, f
	fastProjection();
}

};  // namespace simulator

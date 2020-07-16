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
	log() << "Calculating Segment Length" << std::endl;
	int N = m - 3;
	segmentLength = std::vector<float>(N);
	
	// Compress yarns
	float coeff = params.cInit;

	for (int i = 0; i < N; i++) {
		// Calculate arclength of segment[i].
		// Note that the length is fixed throughout the simulation.

		int index = i * 3;
		DECLARE_POINTS(p, index)

		segmentLength[i] = coeff * integrate([&](float s)->float {
			DECLARE_BASIS_D(b, s);
			return sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2) +
				pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2) +
				pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2));
			}, 0, 1);
	}
}

void Simulator::addSegmentLengthConstraint() {
	int N = m - 3;

	for (int i = 0; i < N; i++) {
		int index = i * 3;
		float curSegmentLength = segmentLength[i];

    auto constraint = [=](const Eigen::MatrixXf& q)->float {
      DECLARE_POINTS(p, index)
      float currentLength = integrate([&](float s)->float {
        DECLARE_BASIS_D(b, s);
        return sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2) +
          pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2) +
          pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2));
        }, 0, 1);
      return 1 - currentLength / curSegmentLength;
    };

		std::vector<Constraints::Entry> constraintD{
      Constraints::Entry {    // px1
        index,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD1 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // py1
        index + 1,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD1 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // pz1
        index + 2,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD1 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // px2
        index + 3,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD2 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // py2
        index + 4,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD2 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // pz2
        index + 5,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD2 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // px3
        index + 6,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD3 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // py3
        index + 7,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD3 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // pz3
        index + 8,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD3 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // px4
        index + 9,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
						DECLARE_BASIS_D(b, s);
            return bD4 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // py4
        index + 10,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
            DECLARE_BASIS_D(b, s);
            return bD4 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      },
      Constraints::Entry {    // pz4
        index + 11,
        [=](const Eigen::MatrixXf& q)->float {
          DECLARE_POINTS(p, index)
          float lengthD = integrate([&](float s)->float {
						DECLARE_BASIS_D(b, s);
            return bD4 * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4,2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4,2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4,2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4);
          }, 0, 1);
          return -lengthD / curSegmentLength;
        }
      }
    };

		constraints.addConstraint(constraint, constraintD);
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

void Simulator::calculateGradient() {
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

	for (int i = 0; i < N; i++) {
		if (i % 10 == 0) {
      log() << "- Collision Energy (" << i << "/" << N << ")" << std::endl;
		}
		for (int j = i + 2; j < N; j++) {
			calculateCollisionEnergyGradient(i, j);
		}
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
	while ((cValue = maxCoeff(constraint = constraints.calculate(qj))) > eps && iter < 10) {
		log() << "- iter: " << iter << ", constraint: " << cValue << std::endl;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;

		Eigen::SparseMatrix<float> jC = constraints.getJacobian(qj);

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

Simulator::Simulator(file_format::YarnRepr yarns, SimulatorParams params_) : params(params_), constraints(0), stepCnt(0) {
	log() << "Initializing Simulator" << std::endl;
	this->yarns = yarns;
	if (yarns.yarns.size() > 0) {
		q = yarns.yarns[0].points;
	}
	m = q.rows();
	history.push_back(yarns);

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
	addSegmentLengthConstraint();

	// Construct Mass Matrix
	constructMassMatrix();

	// Pin first and last control points
	constraints.addPinConstrain(0, q(0), q(1), q(2));
	constraints.addPinConstrain(m - 1, q((m - 1) * 3), q((m - 1) * 3 + 1), q((m - 1) * 3 + 2));

	writeToFile();

	log() << "Simulator Initialized" << std::endl;
}


void Simulator::writeToFile() const {
	writeMatrix("q-" + std::to_string(stepCnt) + ".csv", q);
	writeMatrix("qD-" + std::to_string(stepCnt) + ".csv", qD);
	writeMatrix("gradE-" + std::to_string(stepCnt) + ".csv", gradE);
	//writeMatrix("contactE-" + std::to_string(stepCnt) + ".csv", contactE);
}

void Simulator::step() {
	stepCnt++;
	log() << "Step (" << stepCnt << ")" << std::endl;

	gradE.setZero();
	gradD.setZero();

	// update gradE, gradD, f
	calculateGradient();

	fastProjection();

	// save to YarnRepr
	// supports one yarn for now
	// yarns.yarns[0].points = inflate(q, 3);
	history.push_back(history.back().createAlike());
	history.back().yarns[0].points = inflate(q, 3);

	writeToFile();

	log() << "Done Step (" << stepCnt << ")" << std::endl;
}

std::ostream& Simulator::log() const {
	char buf[20];
	time_t rawTime;
	struct tm* timeInfo;

	time(&rawTime);
	timeInfo = localtime(&rawTime);

	strftime(buf, 20, "%D %T", timeInfo);
	return std::cout << "[" << buf << "] ";
}

};  // namespace simulator

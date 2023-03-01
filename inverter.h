#pragma once
#include <iostream>
#include <cassert>
#include <iomanip>
#include "matrix.h"
#include "utils.h"

#define assertm(exp, msg) assert(((void)msg, exp));

using namespace utils;

class Inverter {
protected:
	const int n;
	Matrix invertedMatrix;
	Matrix invertingMatrix;

	void stepReverseGauss(int k) {
		double m_kk = 1. / invertingMatrix[k * (n + 1)];
		for (int i = 0; i < k; ++i) {
			double coef = -(invertingMatrix[i * n + k] * m_kk);
			invertedMatrix.addRow(i, k, coef);
		}
		invertedMatrix.multiplicateRow(k, m_kk);
	}

	int chooseMaxRow(int k) {
		int kNew = invertingMatrix.argmaxRow(k);
		if (kNew < 0) return -1;

		invertingMatrix.swapRows(k, kNew, k);
		invertedMatrix.swapRows(k, kNew);

		return 0;
	}

	int chooseNonzeroRow(int k) {
		if (!isEqual(invertingMatrix[k * n + k], 0.)) return 0;
		for (int i = k + 1; i < n; ++i) {
			if (abs(invertingMatrix[i * n + k]) > EPS) {
				invertingMatrix.swapRows(k, i);
				invertedMatrix.swapRows(k, i);
				return 0;
			}
		}
		return -1;
	}

	virtual int stepDirectAlgorithm(int k) = 0;

public:
	Inverter(std::shared_ptr<Matrix> argument) noexcept : n((*argument.get()).n), invertedMatrix(Matrix(n, n, &e)), invertingMatrix(*argument.get()) {
		assertm((*argument.get()).n == (*argument.get()).m, "Only square matrices are invertible!");
	}

	std::unique_ptr<Matrix> returnResult() {
		return std::make_unique<Matrix>(invertedMatrix);
	}
		
	void reverseGauss() {
		for (int k = n - 1; k > -1; --k) {
			stepReverseGauss(k);
		}
	}

	void checkRes(const Matrix& copied) {
		std::cout << std::endl << "Inverse * Inverting:" << std::endl;
		Matrix product = (invertedMatrix * copied);
		product.printTrunkated(10, 1. / 100);

		double norm = 0;
		for (int i = 0; i < n; ++i) {
			double temp = 0;
			for (int j = 0; j < n; ++j) {
				int isDiag = (int)(i == j);
				temp += abs(product[i * n + j] - isDiag);
			}
			if (temp > norm) {
				norm = temp;
			}
		}

		std::cout << "Residual of inverse matrix = " << norm << std::endl;;
	}

	void directAlgorithm() {
		for (int k = 0; k < n; ++k) {
			if (stepDirectAlgorithm(k) < 0) {
				std::cout << "Matrix is irrevertible!";
				break;
			}
		}

		std::cout << std::endl << "Triagonalized matrix: " << std::endl;
		invertingMatrix.printTrunkated(10, 1./100);
	}
};

class GaussianInverter : public virtual Inverter {
protected:
	int stepDirectAlgorithm(int k) {
		int irreversibleError = chooseMaxRow(k);
		if (irreversibleError < 0) return irreversibleError;

		double m_kk = 1. / invertingMatrix[k * n + k];

		for (int i = k + 1; i < n; ++i) {
			double coef = -(invertingMatrix[i * n + k] * m_kk);
			invertingMatrix.addRow(i, k, coef, k);
			invertedMatrix.addRow(i, k, coef);
		}

		return 0;
	}

public:
	GaussianInverter(std::shared_ptr<Matrix> argument) noexcept : Inverter(argument) {
		std::cout << "Gaussian inverter created" << std::endl;
		std::cout << "Matrix to be inverted" << std::endl;
		invertingMatrix.printTrunkated(10);
	}
};

class RotationInverter : public virtual Inverter {
protected:
	int stepDirectAlgorithm(int k) {
		int irreversibleError = chooseNonzeroRow(k);
		if (irreversibleError < 0) return irreversibleError;

		for (int i = k + 1; i < n; ++i) {
			double x = invertingMatrix[k * n + k];
			double y = invertingMatrix[i * n + k];
			if (isEqual(y, 0.)) continue;
			double r = Q_rsqrt(x * x + y * y);
			double cos = x * r, sin = - y * r;
			invertingMatrix.rotateRow(k, i, cos, sin, k);
			invertedMatrix.rotateRow(k, i, cos, sin);
		}

		return 0;
	}

public:
	RotationInverter(std::shared_ptr<Matrix> argument) noexcept : Inverter(argument) {
		std::cout << "Rotation inverter created" << std::endl;
		std::cout << "Matrix to be inverted" << std::endl;
		invertingMatrix.printTrunkated(10);
	}
};

class ReflectionInverter : public virtual Inverter {
protected:
	std::shared_ptr<double[]> supportingVector = std::make_shared<double[]>(n);

	inline void refillSupportingVector(int k, double aNorm, double xNorm) noexcept {
		supportingVector[k] = (invertingMatrix[k * n + k] - aNorm) * xNorm;
		for (int i = k + 1; i < n; ++i) {
			supportingVector[i] = invertingMatrix[i * n + k] * xNorm;
		}
	}

	int stepDirectAlgorithm(int k) {
		int irreversibleError = chooseNonzeroRow(k);
		if (irreversibleError < 0) return irreversibleError;

		double trunkatedNorm = 0;
		double aNorm, xNorm;
		double a_kk = invertingMatrix[k * n + k];

		for (int i = k+1; i < n; ++i) {
			trunkatedNorm += invertingMatrix[i * n + k] * invertingMatrix[i * n + k];
		}

		aNorm = std::sqrt(trunkatedNorm + a_kk * a_kk);
		xNorm = 1. / std::sqrt(trunkatedNorm + (a_kk - aNorm) * (a_kk - aNorm));

		refillSupportingVector(k, aNorm, xNorm);

		invertedMatrix.reflectSubmatrix(k, supportingVector);
		invertingMatrix.reflectSubmatrix(k, supportingVector);
	}

public:
	ReflectionInverter(std::shared_ptr<Matrix> argument) noexcept : Inverter(argument) {
		std::cout << "Reflection inverter created" << std::endl;
		std::cout << "Matrix to be inverted" << std::endl;
		invertingMatrix.printTrunkated(10);
	}
};
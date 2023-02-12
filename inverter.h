#pragma once
#include <iostream>
#include <cassert>
#include <iomanip>
#include <vector>
#include "matrix.h"
#include "utils.h"

#define assertm(exp, msg) assert(((void)msg, exp));
#define invertingMatrix (*invertingPtr.get())

using namespace utils;

class Inverter {
protected:
	const int n;
	bool notInvertible = false;
	Matrix invertedMatrix;
	std::shared_ptr<Matrix> invertingPtr;

public:
	Inverter(std::shared_ptr<Matrix> argument) noexcept : n((*argument.get()).n), invertedMatrix(Matrix(n, n, &e)), invertingPtr(argument) {
		assertm((*argument.get()).n == (*argument.get()).m, "Only square matrices are invertible!");
	}

	std::unique_ptr<Matrix> returnResult() {
		return std::move(std::make_unique<Matrix>(invertedMatrix));
	}
		
	void reverseGauss() {
		for (int k = n - 1; k > 0; --k) {
			stepReverseGauss(k);
		}
	}

	void checkRes(const Matrix& copied) {
		std::cout << std::endl << "Inverse * Inverting:" << std::endl;
		(invertedMatrix * copied).printTrunkated(10);
	}

	void stepReverseGauss(int k) {
		double m_kk = 1. / invertingMatrix[k * (n + 1)];
		for (int i = 0; i < k; ++i) {
			double coef = - (invertingMatrix[i * n + k] * m_kk);
			invertedMatrix.addRow(i, k, coef);
		}
		invertedMatrix.multiplicateRow(k, m_kk);
	}

	int chooseMaxRow(int k) {
		int kNew = invertingMatrix.argmaxRow(k);
		if (kNew < 0) return -1;

		invertingMatrix.swapRows(k, kNew);
		invertedMatrix.swapRows(k, kNew);

		return 0;
	}

	void directAlgorithm() {
		for (int k = 0; k < n; ++k) {
			if (stepDirectAlgorithm(k) < 0) {
				std::cout << "Matrix is irrevertible!";
				break;
			}
		}
	}

	virtual int stepDirectAlgorithm(int k) = 0;
};

class GaussianInverter : public virtual Inverter {
public:
	GaussianInverter(std::shared_ptr<Matrix> argument) noexcept : Inverter(argument) {
		assertm((*argument.get()).n == (*argument.get()).m, "Only square matrices are invertible!");

		std::cout << "Gaussian inverter created" << std::endl;
		std::cout << "Matrix to be inverted" << std::endl;
		invertingMatrix.printTrunkated(10);
	}

	int stepDirectAlgorithm(int k) {
		int errorCode = chooseMaxRow(k);
		if (errorCode < 0) return errorCode;

		double m_kk = 1. / invertingMatrix[k * n + k];

		for (int i = k + 1; i < n; ++i) {
			double coef = - (invertingMatrix[i * n + k] * m_kk);
			invertingMatrix.addRow(i, k, coef, k);
			invertedMatrix.addRow(i, k, coef);
		}

		return 0;
	}
};

class RotationInverter : public virtual Inverter {
public:
	RotationInverter(std::shared_ptr<Matrix> argument) noexcept : Inverter(argument) {
		assertm((*argument.get()).n == (*argument.get()).m, "Only square matrices are invertible!");

		std::cout << "Rotation inverter created" << std::endl;
		std::cout << "Matrix to be inverted" << std::endl;
		invertingMatrix.printTrunkated(10);
	}

	int stepDirectAlgorithm(int k) {
		int kNew = invertingMatrix.argmaxRow(k);
		if (kNew < 0) return -1;

		invertingMatrix.swapRows(k, kNew);
		invertedMatrix.swapRows(k, kNew);

		double x = invertingMatrix[k * n + k];
		for (int i = k; i < n; ++i) {
			double y = invertingMatrix[i * n + k];
			double r = std::sqrt(x * x + y * y);
			double cos = x / r, sin = - y / r;
			invertingMatrix.rotate(k, i, cos, sin, k);
			invertedMatrix.rotate(k, i, cos, sin);
		}

		return 0;
	}
};
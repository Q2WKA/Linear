#pragma once
#include <iostream>
#include <iomanip>
#include <memory>
#include "inverter.h"
#include "utils.h"

#define assertm(exp, msg) assert(((void)msg, exp));

class Matrix {
protected:
	const int n, m;
	std::unique_ptr<double[]> data = std::make_unique<double[]>(n * m);
	inline void set(int linearizedIndex, double new_value) {
		data[linearizedIndex] = new_value;
	};

	friend class Inverter;
	friend class RotationInverter;
	friend class GaussianInverter;

public:
	Matrix(int n, int m) : n(n), m(m) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data[i * m + j] = 0;
			}
		}
	}

	Matrix(int n, int m, double (*f)(int, int)) : n(n), m(m) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data[i * m + j] = f(i, j);
			}
		}
	}

	Matrix(const Matrix& other): n(other.n), m(other.m) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data[i * m + j] = other.data[i * m + j];
			}
		}
	}

	Matrix operator=(const Matrix& other) {
		Matrix newMatrix(other.n, other.m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				newMatrix.data[i * m + j] = other.data[i * m + j];
			}
		}
		return newMatrix;
	}

	Matrix operator*(const Matrix& other) const {
		assertm(m == other.n, "Matrices must have fitting sizes!");

		Matrix newMatrix(n, other.m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < other.m; ++j) {
				double m_ij = 0;
				for (int k = 0; k < m; ++k) {
					m_ij += data[i * m + k] * other.data[k * other.m + j];
				}
				newMatrix.data[i * other.m + j] = m_ij;
			}
		}
		return std::move(newMatrix);
	}

	inline double operator[](int linearizedIndex) const {
		return data[linearizedIndex];
	}
		
	inline void addRow(int i_1, int i_2, double coef) noexcept {
		for (int j = 0; j < m; ++j) {
			data[i_1 * m + j] += data[i_2 * m + j] * coef;
		}
	}

	inline void addRow(int i_1, int i_2, double coef, int k) noexcept {
		for (int j = k; j < m; ++j) {
			data[i_1 * m + j] += data[i_2 * m + j] * coef;
		}
	}

	inline void multiplicateRow(int i, double coef) noexcept {
		for (int j = 0; j < m; ++j) {
			data[i * m + j] *= coef;
		}
	}

	inline void swapRows(int i_1, int i_2) {
		double temp;
		for (int j = 0; j < m; ++j) {
			temp = data[i_1 * m + j];
			data[i_1 * m + j] = data[i_2 * m + j];
			data[i_2 * m + j] = temp;
		}
	}

	void print() const noexcept {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				if (utils::isEqual(data[i * m + j], 0.)) std::cout << std::setw(7) << std::setprecision(3) << 0 << "  ";
				else std::cout << std::setw(7) << std::setprecision(3) << data[i * m + j] << "  ";
			}
			std::cout << std::endl;
		}
	}

	void printTrunkated(int threshold) const noexcept {
		int n_ = std::min<int>(threshold, n), m_ = std::min<int>(threshold, m);

		for (int i = 0; i < n_; ++i) {
			for (int j = 0; j < m_; ++j) {
				if (utils::isEqual(data[i * m + j], 0.)) std::cout << std::setw(7) << std::setprecision(3) << 0 << "  ";
				else std::cout << std::setw(7) << std::setprecision(3) << data[i * m + j] << "  ";
			}
			std::cout << std::endl;
		}
	}
};
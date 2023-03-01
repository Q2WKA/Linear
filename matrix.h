#pragma once
#include <iostream>
#include <iomanip>
#include <memory>
#include "utils.h"

#define assertm(exp, msg) assert(((void)msg, exp));

class Matrix {
public:
	const int n, m;

protected:
	std::unique_ptr<double[]> data = std::make_unique<double[]>(n * m);
	inline void set(int linearizedIndex, double new_value) {
		data[linearizedIndex] = new_value;
	};

	friend class Inverter;
	friend class RotationInverter;

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

	Matrix(Matrix&& other) noexcept : n(other.n), m(other.m) {
		data = std::move(other.data);
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
		return newMatrix;
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

	inline void swapRows(int i_1, int i_2, int k) {
		double temp;
		for (int j = k; j < m; ++j) {
			temp = data[i_1 * m + j];
			data[i_1 * m + j] = data[i_2 * m + j];
			data[i_2 * m + j] = temp;
		}
	}

	inline void rotateRow(int i_1, int i_2, double cos, double sin) {
		for (int j = 0; j < m; ++j) {
			double x_i1 = data[i_1 * m + j];
			double x_i2 = data[i_2 * m + j];
			data[i_1 * m + j] = x_i1 * cos - x_i2 * sin;
			data[i_2 * m + j] = x_i1 * sin + x_i2 * cos;
		}
	}

	inline void rotateRow(int i_1, int i_2, double cos, double sin, int k) {
		for (int j = k; j < m; ++j) {
			double x_i1 = data[i_1 * m + j];
			double x_i2 = data[i_2 * m + j];
			data[i_1 * m + j] = x_i1 * cos - x_i2 * sin;
			data[i_2 * m + j] = x_i1 * sin + x_i2 * cos;
		}
	}

	void reflectSubmatrix(int k, std::shared_ptr<double[]> x) {
		for (int j = k; j < m; ++j) {
			double scalarProduct = 0;
			for (int i = k; i < n; ++i) {
				scalarProduct += data[i * m + j] * x[i];
			}
			scalarProduct *= 2;
			for (int i = k; i < n; ++i) {
				data[i * m + j] -= x[i] * scalarProduct;
			}
		}
	}

	int argmaxRow(int k) const {
		double maxValue = data[k * n + k];
		int argmax = k;
		for (int i = k + 1; i < n; ++i) {
			if (abs(data[i * n + k]) > abs(maxValue)) {
				maxValue = data[i * n + k];
				argmax = i;
			}
		}

		if (utils::isEqual(maxValue, 0))
			return -1;

		return argmax;
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

	void print(double zeroPrecision) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				if (utils::isEqual(data[i * m + j] * zeroPrecision, 0.)) std::cout << std::setw(7) << std::setprecision(3) << 0 << "  ";
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

	void printTrunkated(int threshold, double zeroPrecision) const noexcept {
		int n_ = std::min<int>(threshold, n), m_ = std::min<int>(threshold, m);

		for (int i = 0; i < n_; ++i) {
			for (int j = 0; j < m_; ++j) {
				if (utils::isEqual(data[i * m + j] * zeroPrecision, 0.)) std::cout << std::setw(7) << std::setprecision(3) << 0 << "  ";
				else std::cout << std::setw(7) << std::setprecision(3) << data[i * m + j] << "  ";
			}
			std::cout << std::endl;
		}
	}
};
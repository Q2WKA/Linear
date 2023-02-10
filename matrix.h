#pragma once
#include <memory>
#include <iostream>
#include <iomanip>
#include "linear_vector.h"

#define assertm(exp, msg) assert(((void)msg, exp));

class matrix {
protected:
	const int n, m;
	std::unique_ptr<double[]> data = std::make_unique<double[]>(n * m);

public:
	matrix(int n, int m) : n(n), m(m) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data[i * m + j] = 0;
			}
		}
	}

	matrix(int n, int m, double (*f)(int, int)) : n(n), m(m) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data[i * m + j] = f(i, j);
			}
		}
	}

	void print() const noexcept {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				std::cout << std::setw(6) << std::setprecision(3) << data[i * m + j] << "  ";
			}
			std::cout << std::endl;
		}
	}

	inline double operator()(int i, int j) const noexcept {
		return data[i * m + j];
	}

	linear_vector operator*(const linear_vector& vector) const noexcept {
		assertm(m == vector.dim, "Vectors must be of same dim!");

		linear_vector newVector(n);
		for (int i = 0; i < n; ++i) {
			int tempRes = 0;
			for (int j = 0; j < m; ++j) {
				tempRes += vector[j] * data[i * m + j];
			}
			newVector.set(i, tempRes);
		}
		return newVector;
	}
};
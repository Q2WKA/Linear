#include <iostream>
#include <cassert>
#include <memory>
#include <cmath>
#include "matrix.h"

#define assertm(exp, msg) assert(((void)msg, exp));

class linear_vector {
protected:
	const int dim;
	std::unique_ptr<double[]> data = std::make_unique<double[]>(dim);

	double operator[](int i) noexcept {
		assertm(i < dim, "Asked index is out of bonds!");
		return data[i];
	}

	inline void set(int i, double value) {
		data[i] = value;
	}

	friend class matrix;

public:
	linear_vector(int dim) noexcept : dim(dim) {
		for (int i = 0; i < dim; ++i) {
			data[i] = 0;
		}
	}

	linear_vector(const linear_vector& other) noexcept : dim(other.dim) {
		for (int i = 0; i < dim; ++i) {
			data[i] = other.data[i];
		}
	}

	linear_vector& operator=(const linear_vector& other) noexcept {
		assertm(dim == other.dim, "Vectors must be of same dim!");

		for (int i = 0; i < dim; ++i) {
			data[i] = other.data[i];
		}

		return *this;
	}

	linear_vector(linear_vector&& other) noexcept : dim(other.dim) {
		std::swap(data, other.data);
	}

	linear_vector& operator=(linear_vector&& other) noexcept {
		std::swap(data, other.data);

		return *this;
	}

	void fill() noexcept {
		for (int i = 0; i < dim; ++i) {
			std::cout << "Enter " << i + 1 << "-th vector coordinate" << std::endl;
			std::cin >> data[i];
		}
	}

	void print() noexcept {
		for (int i = 0; i < dim; ++i) {
			std::cout << data[i] << "  ";
		}
		std::cout << std::endl;
	}

	bool operator==(const linear_vector& other) noexcept {
		assertm(dim == other.dim, "Vectors must be of same dim!");
		for (int i = 0; i < dim; ++i) {
			if (data[i] != other.data[i]) return false;
		}
		return true;
	}

	linear_vector& operator+=(const linear_vector& other) noexcept {
		assertm(dim == other.dim, "Vectors must be of same dim!");
		for (int i = 0; i < dim; ++i) {
			data[i] += other.data[i];
		}
		return *this;
	}

	inline linear_vector operator+(const linear_vector& other) noexcept {
		assertm(dim == other.dim, "Vectors must be of same dim!");
		linear_vector newVector(dim);
		for (int i = 0; i < dim; ++i) {
			newVector.data[i] = data[i] + other.data[i];
		}
		return newVector;
	}

	linear_vector& operator*=(int lambda) noexcept {
		for (int i = 0; i < dim; ++i) {
			data[i] *= lambda;
		}
		return *this;
	}

	inline linear_vector operator*(int lambda) noexcept {
		linear_vector newVector(dim);
		for (int i = 0; i < dim; ++i) {
			newVector.data[i] = data[i] * lambda;
		}
		return newVector;
	}

	friend double dot(const linear_vector& lhs, const linear_vector& rhs) noexcept {
		assertm(lhs.dim == rhs.dim, "Vectors must be of same dim!");
		double res = 0;
		for (int i = 0; i < lhs.dim; ++i) {
			res += lhs.data[i] * rhs.data[i];
		}
		return res;
	}

	double l_norm(double l) const noexcept {
		double res = 0;
		for (int i = 0; i < dim; ++i) {
			res += pow(data[i], l);
		}
		return res;
	}

	double operator[](int i) const noexcept {
		assertm(i < dim, "Asked index is out of bonds!");
		return data[i];
	}
};
#include <chrono>
#include <iostream>
#include "inverter.h"
#include "matrix.h"
#include "utils.h"

using namespace utils;

int main(int argc, char* argv[]) {
	std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

	Matrix to_invert1 = Matrix(1000, 1000, &revAbsDiff);

	GaussianInverter gaussianInverter(std::make_shared<Matrix>(to_invert1));
	gaussianInverter.directAlgorithm();
	gaussianInverter.reverseGauss();
	// gaussianInverter.checkRes(to_invert);

	std::unique_ptr<Matrix> res = gaussianInverter.returnResult();
	std::cout << std::endl << "Inverted matrix: " << std::endl;
	res.get()->printTrunkated(10, 1. / 100);

	std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[micro seconds]" << std::endl;




	std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

	RotationInverter rotationInverter(std::make_shared<Matrix>(to_invert1));
	rotationInverter.directAlgorithm();
	rotationInverter.reverseGauss();
	// gaussianInverter.checkRes(to_invert);

	res = rotationInverter.returnResult();
	std::cout << std::endl << "Inverted matrix: " << std::endl;
	res.get()->printTrunkated(10, 1. / 100);

	std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[micro seconds]" << std::endl;

    return 0;
}

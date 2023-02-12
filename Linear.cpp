#include <chrono>
#include <iostream>
#include "inverter.h"
#include "matrix.h"
#include "utils.h"

using namespace utils;

int main() {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	Matrix to_invert = Matrix(100, 100, &revAbsDiff);
	Matrix copied = to_invert;

	GaussianInverter gaussianInverter(std::make_shared<Matrix>(to_invert));
	gaussianInverter.directAlgorithm();
	gaussianInverter.reverseGauss();
	gaussianInverter.checkRes(copied);

	std::unique_ptr<Matrix> res = gaussianInverter.returnResult();
	std::cout << std::endl << "Inverted matrix: " << std::endl;
	res.get()->printTrunkated(10);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro seconds]" << std::endl;

    return 0;
}

#include <chrono>
#include <iostream>
#include "inverter.h"
#include "matrix.h"
#include "utils.h"

using namespace utils;

int main(int argc, char* argv[]) {
	Matrix to_invert = Matrix(1000, 1000, &revAbsDiff);



	std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
	GaussianInverter gaussianInverter(std::make_shared<Matrix>(to_invert));
	gaussianInverter.directAlgorithm();
	gaussianInverter.reverseGauss();
	gaussianInverter.checkRes(to_invert);
	std::unique_ptr<Matrix> res1 = gaussianInverter.returnResult();
	std::cout << std::endl << "Inverted matrix: " << std::endl;
	res1.get()->printTrunkated(10, 1. / 100);
	std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[micro seconds]" << std::endl << std::endl << std::endl;




	std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
	RotationInverter rotationInverter(std::make_shared<Matrix>(to_invert));
	rotationInverter.directAlgorithm();
	rotationInverter.reverseGauss();
	rotationInverter.checkRes(to_invert);
	std::unique_ptr<Matrix> res2 = rotationInverter.returnResult();
	res2 = rotationInverter.returnResult();
	std::cout << std::endl << "Inverted matrix: " << std::endl;
	res2.get()->printTrunkated(10, 1. / 100);
	std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[micro seconds]" << std::endl << std::endl << std::endl;




	std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();
	ReflectionInverter reflectionInverter(std::make_shared<Matrix>(to_invert));
	reflectionInverter.directAlgorithm();
	reflectionInverter.reverseGauss();
	reflectionInverter.checkRes(to_invert);
	std::unique_ptr<Matrix> res3 = rotationInverter.returnResult();
	res3 = rotationInverter.returnResult();
	std::cout << std::endl << "Inverted matrix: " << std::endl;
	res3.get()->printTrunkated(10, 1. / 100);
	std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end3 - begin3).count() << "[micro seconds]" << std::endl;

    return 0;
}

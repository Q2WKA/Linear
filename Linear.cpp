#include <iostream>
#include "matrix.h"
#include "utils.h"

using namespace utils;

int main() {
	matrix m(3, 5, &abs_diff);

	linear_vector v(5);

	v.fill();

	linear_vector v1 = m*v;

	v1.print();

	m.print();

    return 0;
}

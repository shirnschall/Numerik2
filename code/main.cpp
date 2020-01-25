#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>

#include "Linag/Eigen/Dense"

int main() {

    srand(5);

    linag::Vector<double> test{1,2,3,4,5,6};

    std::cout << test << std::endl << test.l2norm();


    return 0;
}

#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>

#include "Linag/Eigen/Dense"

int main() {

    srand(5);

    linag::DenseMatrix<double> test(200,200);
    linag::DenseMatrix<double> pinv(200,200);
    linag::Vector<double> b(200);
    b.rand();
    test.randSPD(9);
    pinv.diag(test);
    pinv = pinv.inverse();

    linag::SparseMatrix<double> testS(test);
    linag::SparseMatrix<double> pinvS(pinv);

    int count = 0;


    std::cout << test.conjugateGradientSolver(b,10e-10, &count, nullptr) <<
              std::endl << count << std::endl;

    std::cout << testS.preCondConjugateGradientSolver(pinvS,b,10e-10, &count, nullptr) <<
                std::endl << count << std::endl;



    return 0;
}

#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <ctime>


int main() {

    srand(5);

    std::cout << "n\tt" << std::endl;
    for (int i = 1; i<100000; i*=2) {
        linag::DenseMatrix<double> test{i,i};

        test.rand();
        std::clock_t begin = std::clock();
        test = test*test.transpose();
        std::clock_t end = std::clock();

        std::cout << "" << i << "\t" << double(end - begin) / CLOCKS_PER_SEC << std::endl;

    }

    //test.randSPD(3);
    //linag::SparseMatrix<double> tests(test);
    //std::cout << b <<std::endl;
    //std::cout << tests*tests.conjugateGradientSolver(b,0.000001) << std::endl;
    //std::cout << test.conjugateGradientSolver(b,0.000001) << std::endl;
    //std::cout << b.l2norm() << test << std::endl;
    //linag::SparseMatrix<double> a;
    //linag::DenseMatrix<double> d(a);


    return 0;
}

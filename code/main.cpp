#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>


int main() {

    srand(5);

    //measure time with <chrono>:
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //...

    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    const int notZero = 4;
        linag::DenseMatrix<double> test(10, 10);
        linag::Vector<double> b(test.dim().rows);

        test.randSPD(4);
        b.rand();




    //test.randSPD(3);
    linag::SparseMatrix<double> tests(test);
    std::cout << b <<std::endl;
    std::cout << tests*tests.conjugateGradientSolver(b,0.000001) << std::endl;
    std::cout << test.conjugateGradientSolver(b,0.000001) << std::endl;
    //std::cout << b.l2norm() << test << std::endl;
    //linag::SparseMatrix<double> a;
    //linag::DenseMatrix<double> d(a);


    return 0;
}

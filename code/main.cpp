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

        for (int i = 10; i < 10000; i+=100) {
            linag::DenseMatrix<double> test(i, i);
            linag::Vector<double> b(test.dim().rows);

            test.randSPD(notZero>i?i: notZero);
            b.rand();
            linag::SparseMatrix<double> testS(test);

            std::cout << "n = " << i << std::endl;


            std::cout << "DenseMatrix: " << std::endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                //test.conjugateGradientSolver(b,10e-10);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
            std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
            std::cout << "SparseMatrix: " << std::endl;
            begin = std::chrono::steady_clock::now();
                testS.conjugateGradientSolver(b,10e-10);
            end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
            std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl << std::endl;

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

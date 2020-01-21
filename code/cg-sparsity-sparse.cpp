#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>


int main() {
    //seed std::rand()
    srand(5);

    const int notZero = 5;
    //csv header
    std::cout << "n\tt\ttype\tdensity\n";

        for (int n = 10; n < 1000; n+=50) {
            for (double j = 0.00; j < 1; j+=0.2) { // density in %

            //create random symmetric positive definite matrix A and vector b
            linag::DenseMatrix<double> A(n, n);
            linag::Vector<double> b(A.dim().rows);

            A.randSPD((int)std::floor(n*j)%2?std::floor(n*j):std::floor(n*j)+1);
            b.rand();

            linag::SparseMatrix<double> ASparse(A);


            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            ASparse.conjugateGradientSolver(b,10e-10);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << n << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\tsparse\t" << j << std::endl;

            }

        }

    return 0;
}

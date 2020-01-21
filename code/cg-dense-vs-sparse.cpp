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
    std::cout << "n\tt\ttype\n";

    for (int n = 10; n < 1000; n+=50) {
        //create random symmetric positive definite matrix A and vector b
        linag::DenseMatrix<double> A(n, n);
        linag::Vector<double> b(A.dim().rows);

        A.randSPD(notZero>n?n:notZero);
        b.rand();

        linag::SparseMatrix<double> ASparse(A);


        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        A.conjugateGradientSolver(b,10e-10);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << n << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\tdense\n";

        begin = std::chrono::steady_clock::now();
        ASparse.conjugateGradientSolver(b,10e-10);
        end = std::chrono::steady_clock::now();
        std::cout << n << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "\tsparse\n";


    }
    return 0;
}

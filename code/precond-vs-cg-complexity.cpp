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

    for (int n = 10; n < 2000; n+=50) {
        //create random symmetric positive definite matrix A and vector b
        linag::DenseMatrix<double> A(n, n);
        linag::Vector<double> b(A.dim().rows);
        linag::DenseMatrix<double> pinv(n,n);

        A.randSPD(notZero>n?n:notZero);
        b.rand();
        for (int i = 0; i < pinv.dim().rows; ++i) {
            for (int j = 0; j < pinv.dim().cols; ++j) {
                if(i==j)
                    pinv.at(i,j) = 1/A.at(i,j);
                else
                    pinv.at(i,j) = 0;
            }
        }
        linag::SparseMatrix<double> pinvS(pinv);

        linag::SparseMatrix<double> ASparse(A);

        int count = 0;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        ASparse.conjugateGradientSolver(b,10e-10,&count);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << n << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/count << "\tcg\n";

        begin = std::chrono::steady_clock::now();
        ASparse.preCondConjugateGradientSolver(pinvS,b,10e-10,&count);
        end = std::chrono::steady_clock::now();
        std::cout << n << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/count << "\tprecond-cg\n";


    }
    return 0;
}

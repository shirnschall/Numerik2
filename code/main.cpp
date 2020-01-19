#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"

int main() {

    srand(5);

    linag::DenseMatrix<double> test(4,4);
    linag::Vector<double> b{1,2,3,4};
    linag::Vector<double> z{1e-10,1e-11,1e-12,0.00001};
    //b.rand();
    //test.id();
    //test.at(0,3)=0.2;
    //test.at(0,2) = 17.2;
    //test.randLT();
    test.randSPD(4);
    std::cout << test <<std::endl;
    std::cout << test*test.conjugateGradientSolver(b,0.000001) << std::endl;
    std::cout << test.conjugateGradientSolver(b,0.000001) << std::endl;
    //std::cout << b.l2norm() << test << std::endl;
    //linag::SparseMatrix<double> a;

    //linag::DenseMatrix<double> d(a);
    return 0;
}

#include <iostream>
#include <complex>
#include "densematrix.h"
//#include "vector.h"

int main() {

    srand(5);

    linag::DenseMatrix<double> test {{1,1,1,1}, {1,2,2,2}, {2,3,4,5}, {6,3,5,4}};
    test = test*test.transpose();
    linag::Vector<double> b{1,2,3,4};
    linag::Vector<double> z{1e-10,1e-11,1e-12,0.00001};
    //b.rand();
    //test.id();
    //test.at(0,3)=0.2;
    //test.at(0,2) = 17.2;
    test.conjugateGradientSolver(b,0.000001);
    //std::cout << b.l2norm() << test << std::endl;
    return 0;
}

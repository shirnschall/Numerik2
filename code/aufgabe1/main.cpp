#include <iostream>
#include <complex>
#include "densematrix.h"
int main() {

    linag::DenseMatrix<double> test {{1,0,17.2,3}, {0,1.2,5,4}, {0,0,4,0}, {0,0,0,1}};
    //test.at(0,2) = 17.2;
    std::cout << test << std::endl;
    //test.inverse();
    std::cout << std::endl << test*test.inverse() <<std::endl;
    return 0;
}

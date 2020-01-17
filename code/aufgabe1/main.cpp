#include <iostream>
#include <complex>
#include "densematrix.h"
int main() {

    linag::DenseMatrix<std::complex<double>> test(2,2);



    std::cout << test << std::endl;


    return 0;
}

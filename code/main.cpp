#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>

#include "Linag/Eigen/Dense"

int main() {

    srand(5);

    linag::DenseMatrix<double> test(501,501);
    test.randSPD(1);
    linag::DenseMatrix<double> test2(501,501);
    test2.randSPD(501);
    linag::DenseMatrix<double> test3(501,501);
    test3.randSPD(21);
    linag::Vector<double> b(test.dim().rows);
    b.rand();


    std::cout << test.cond() << std::endl << std::endl;

    linag::Vector<double> rs(test.dim().rows);
    linag::Vector<double> rs2(test.dim().rows);
    linag::Vector<double> rs3(test.dim().rows);
    test.conjugateGradientSolver(b,10e-10, nullptr,&rs);
    test2.conjugateGradientSolver(b,10e-10, nullptr,&rs2);
    test3.conjugateGradientSolver(b,10e-10, nullptr,&rs3);

    std::cout << "solved\n";

    double testcond = test.cond();
    double test2cond = test2.cond();
    double test3cond = test3.cond();

    std::cout << "k\tr\ti\n";
    for (int i = 0; i < rs.length(); ++i) {
        std::cout << testcond << '\t' << rs.at(i) << '\t' << i << '\n';
        std::cout << test2cond << '\t' << rs2.at(i) << '\t' << i << '\n';
        std::cout << test3cond << '\t' << rs3.at(i) << '\t' << i << '\n';
    }


    return 0;
}

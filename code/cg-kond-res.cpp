#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>

#include "Linag/Eigen/Dense"

int main() {

    srand(time(NULL));

    linag::DenseMatrix<double> test(501,501);
    test.randSPD(21);
    linag::DenseMatrix<double> test2(501,501);
    test2.randSPD(501);
    linag::DenseMatrix<double> test3(501,501);
    test3.randSPD(21);
    linag::Vector<double> b(test.dim().rows);
    b.rand();


    //std::cout << test.cond() << std::endl << std::endl;

    linag::Vector<linag::Vector<double>*> xs(test.dim().rows);

    linag::Vector<double> res (test.conjugateGradientSolver(b,10e-10, nullptr,&xs));


    //std::cout << "solved\n";

    double testcond = test.cond();

    std::cout << "k\tr\ti\n";
    for (int i = 0; i < xs.length(); ++i) {
        if(!xs.at(i))
            break;
        std::cout << testcond << '\t' << ((*xs.at(i))-res).Anorm(test) << '\t' << i << '\n';
        delete(xs.at(i));

    }


    std::cout << std::endl;
    return 0;
}

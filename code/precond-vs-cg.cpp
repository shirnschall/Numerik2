#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>

#include "Linag/Eigen/Dense"

int main() {

    srand(5);

    linag::DenseMatrix<double> test(201,201);
    linag::DenseMatrix<double> pinv(201,201);
    linag::Vector<double> b(201);
    b.rand();
    test.randSPD(9);
    pinv.diag(test);
    pinv = pinv.inverse();

    linag::SparseMatrix<double> testS(test);
    linag::SparseMatrix<double> pinvS(pinv);

    int count = 0;

    linag::Vector<double> rs(b.length());
    linag::Vector<linag::Vector<double>*> xs(b.length());
    linag::Vector<double> res(b.length());

    //std::cout << test.cond() << std::endl;

    std::cout << "r\tt\ttype" << std::endl;


    res = testS.conjugateGradientSolver(b,10e-16, &count, &xs,&rs);
    for (int i = 0; i < rs.length(); ++i) {
        if(rs.at(i))
            std::cout << (*xs.at(i) - res).Anorm(test) << '\t' << i << "\tcg" << std::endl;
    }

    res = testS.preCondConjugateGradientSolver(pinvS, b,10e-16, &count, &xs,&rs);
    for (int i = 0; i < rs.length(); ++i) {
        if(rs.at(i))
            std::cout << (*xs.at(i) - res).Anorm(test) << '\t' << i << "\tprecg" << std::endl;
    }


    return 0;
}

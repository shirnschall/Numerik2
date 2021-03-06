#include <iostream>
#include "Linag/sparsematrix.h"
#include "Linag/densematrix.h"
#include "Linag/vector.h"
#include <chrono>

#include "Linag/Eigen/Dense"

int main() {

    srand(5);


    std::cout << "i\tn\ttype" << std::endl;
    for (int i = 10; i < 1000; i+=10) {
    linag::DenseMatrix<double> test(i,i);
    linag::DenseMatrix<double> pinv(test.dim());
    linag::DenseMatrix<double> id(test.dim());
    id.id();
    linag::Vector<double> b(test.dim().rows);
    b.rand();

    int notZeroPerLine = 9;

    //create matrix
    test.zeros();
    test.randLT();
    test=test+id;

    //test = test + test.transpose() + c * diagM;
    test = test * test.transpose();
    test.randSPD(9);

    //std::cout << test.toEigen().eigenvalues() << std::endl;

    pinv.diag(test);
    pinv = pinv.inverse();

    linag::SparseMatrix<double> testS(test);
    linag::SparseMatrix<double> pinvS(pinv);

    int count = 0;

    linag::Vector<double> rs(b.length());
    linag::Vector<linag::Vector<double>*> xs(b.length());
    linag::Vector<double> res(b.length());

    //std::cout << test.cond() << std::endl;

    res = test.conjugateGradientSolver(b,10e-16);
    res = testS.conjugateGradientSolver(b,10e-16, &count, &xs,&rs);
    std::cout << count << '\t' << i << "\tcg" << std::endl;

    res = testS.preCondConjugateGradientSolver(pinvS, b,10e-16, &count, &xs,&rs);
            std::cout << count<< '\t' << i << "\tprecg" << std::endl;


    }

    return 0;
}

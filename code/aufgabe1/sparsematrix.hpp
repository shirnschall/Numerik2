//
// Created by Sebastian Hirnschall on 30.12.19.
//

#ifndef AUFGABE1_SPARSEMATRIX_HPP
#define AUFGABE1_SPARSEMATRIX_HPP

#include <stdlib.h>

namespace linag {
    template<typename T>
    struct vector {
        int length;
        T *data;
    };

    double **generateLSData(int);

    vector<double> conjugateGradientSolver(double ** A, int n, double * b, double * x, double tau);


    class SparseMatrix {
    private:
        linag::vector<double> v;
        linag::vector<int> J;
        linag::vector<int> I;

        int sum(int *, int);

    public:
        SparseMatrix(int, int *, int);

        ~SparseMatrix();

        SparseMatrix &operator=(const SparseMatrix &);

        SparseMatrix(const SparseMatrix &);

        SparseMatrix(const double **, int, int);


        const linag::vector<double> &getv() const { return v; }

        const linag::vector<int> &getI() const { return I; }

        const linag::vector<int> &getJ() const { return J; };

        linag::vector<double> &setv() { return v; }

        linag::vector<int> &setI() { return I; }

        linag::vector<int> &setJ() { return J; };

        SparseMatrix* transpose();

        int size() const;

    };

}
#endif //AUFGABE1_SPARSEMATRIX_HPP

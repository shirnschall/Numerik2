//
// Created by Sebastian Hirnschall on 30.12.19.
//

#ifndef AUFGABE1_SPARSEMATRIX_HPP
#define AUFGABE1_SPARSEMATRIX_HPP

#include <stdlib.h>

template <typename T>
struct array {
    int length;
    T* data;
};

class SparseMatrix {
private:
    array<double> v;
    array<int> J;
    array<int> I;

    int sum(int*,int);

public:
    SparseMatrix(int,int*,int);
    ~SparseMatrix();

    SparseMatrix& operator=(const SparseMatrix&);
    SparseMatrix(const SparseMatrix&);
    SparseMatrix(const double**,int,int);


    const array<double>& getv() const { return v;}
    const array<int>& getI() const { return I;}
    const array<int>& getJ() const { return J;};
    int size() const;

};


#endif //AUFGABE1_SPARSEMATRIX_HPP

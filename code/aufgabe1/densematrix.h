//
// Created by Sebastian Hirnschall on 17.01.20.
//

#ifndef AUFGABE1_DENSEMATRIX_H
#define AUFGABE1_DENSEMATRIX_H

#include "linag.h"
//eigen lib
#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Eigenvalues>
#include <iostream>

namespace linag {

template<typename T>
class DenseMatrix{
private:
    size dimension;
    T* data;
public:
    DenseMatrix(int rows,int cols);
    DenseMatrix(linag::size dimension);
    ~DenseMatrix();

    DenseMatrix(const DenseMatrix<T> &);
    DenseMatrix<T> &operator=(const DenseMatrix<T> &rhs);

    operator Eigen::Matrix() const;

    const DenseMatrix<T> operator-() const;

    const DenseMatrix<T> inverse();
    const DenseMatrix<T> transpose();


    T& at(int,int);
    const T& at(int,int) const;
    Vector<T> colToVector(int col);
    Vector<T> rowToVector(int col);
    size dim();
};

template<typename T>
const DenseMatrix<T> operator+(const DenseMatrix<T>& x,const DenseMatrix<T>& y);
template<typename T>
const DenseMatrix<T> operator-(const DenseMatrix<T>& x,const DenseMatrix<T>& y);
template<typename T>
const DenseMatrix<T> operator*(const DenseMatrix<T>& x,const DenseMatrix<T>& y);
template<typename T>
const DenseMatrix<T> operator*(const DenseMatrix<T>& x,const T y);
template<typename T>
const DenseMatrix<T> operator*(const T x,const DenseMatrix<T>& y);
template<typename T>
const Vector<T> operator*(const Vector<T>& x,const DenseMatrix<T>& y);
template<typename T>
const Vector<T> operator*(const DenseMatrix<T>& x,const Vector<T>& y);

template<typename T>
std::ostream& operator<<(std::ostream& output,const DenseMatrix<T>& x);


}

#endif //AUFGABE1_DENSEMATRIX_H

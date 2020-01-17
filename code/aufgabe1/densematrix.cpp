//
// Created by Sebastian Hirnschall on 17.01.20.
//

#include "densematrix.h"

template<typename T>
std::ostream& linag::operator<<(std::ostream& output,const linag::DenseMatrix<T>& x){
    for (int i = 0; i < x.dim().rows; ++i) {
        for (int j = 0; j < x.dim().cols; ++j) {
            output << x.at(i,j);
        }
        output << '\n';
    }
    return output;
}

template<typename T>
const linag::Vector<T> linag::operator*(const linag::DenseMatrix<T>& x,const linag::Vector<T>& y){
    assert(y.dim() == x.dim().cols);

    linag::Vector<T> res(x.dim().rows);

    for (int i = 0; i < res.dim(); ++i) {
        res.at(i) = 0;
        for (int j = 0; j < y.dim(); ++j) {
            res.at(i) += x.at(i,j) * y.at(j);
        }
    }
    return res;
}

template<typename T>
const linag::Vector<T> operator*(const Vector<T>& x,const linag::DenseMatrix<T>& y){
    assert(y.dim() == x.dim().rows);

    linag::Vector<T> res(x.dim().cols);

    for (int i = 0; i < res.dim(); ++i) {
        res.at(i) = 0;
        for (int j = 0; j < y.dim(); ++j) {
            res.at(i) += x.at(j,i) * y.at(j);
        }
    }
    return res;
}

template<typename T>
const linag::DenseMatrix<T> operator*(const T x,const linag::DenseMatrix<T>& y){
    linag::DenseMatrix<T> res(y.dim());
    for (int i = 0; i < y.dim().rows; ++i) {
        for (int j = 0; j < y.dim().cols; ++j) {
            res.at(i,j) = x*y.at(i,j);
        }
    }
    return res;
}

template<typename T>
const linag::DenseMatrix<T> operator*(const linag::DenseMatrix<T>& x,const T y){
    return y*x;
}

template<typename T>
const linag::DenseMatrix<T> operator*(const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y){
    assert(x.dim().cols==y.dim().rows);

    linag::DenseMatrix<T> res(x.dim().rows,y.dim().cols);

    for (int i = 0; i < res.dim().rows; ++i) {
        for (int j = 0; j < res.dim().cols; ++j) {
            res.at(i,j)=0;
            for (int k = 0; k < x.dim().cols; ++k) {
                res.at(i,j) += x.at(i,k) * y.at(k,j);
            }
        }
    }
    return res;
}


template<typename T>
const linag::DenseMatrix<T> operator-(const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y){
    assert(x.dim().cols==y.dim().cols && x.dim().rows==y.dim().rows);

    linag::DenseMatrix<T> res(x.dim().rows,x.dim().cols);

    for (int i = 0; i < x.dim().rows; ++i) {
        for (int j = 0; j < x.dim().cols; ++j) {
            res.at(i,j) = x.at(i,j)-y.at(i,j);
        }
    }
    return res;
}

template<typename T>
const linag::DenseMatrix<T> operator+(const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y){
    return x-(-y);
}

template<typename T>
linag::size linag::DenseMatrix<T>::dim(){
    return dimension;
}


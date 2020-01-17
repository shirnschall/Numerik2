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

    DenseMatrix(const DenseMatrix<T> &rhs);
    DenseMatrix<T> &operator=(const DenseMatrix<T> &rhs);

    operator Eigen::MatrixXd() const;

    const DenseMatrix<T> operator-() const;

    const DenseMatrix<T> inverse() const;
    const DenseMatrix<T> transpose() const;


    T& at(int row,int col);
    const T& at(int row,int col) const;
    Vector<T> colToVector(int col);
    Vector<T> rowToVector(int row);
    const size dim() const;
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

//
// Created by Sebastian Hirnschall on 17.01.20.
//

#include "densematrix.h"
#include <string.h>

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

    auto res = linag::Vector<T>(x.dim().rows);

    for (int i = 0; i < res.dim(); ++i) {
        res.at(i) = 0;
        for (int j = 0; j < y.dim(); ++j) {
            res.at(i) += x.at(i,j) * y.at(j);
        }
    }
    return res;
}

template<typename T>
const linag::Vector<T> operator*(const linag::Vector<T>& x,const linag::DenseMatrix<T>& y){
    assert(y.dim() == x.dim().rows);

    auto res = linag::Vector<T>(x.dim().cols);

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
    auto res = linag::DenseMatrix<T>(y.dim());
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

    auto res = linag::DenseMatrix<T>(x.dim().rows,y.dim().cols);

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

    auto res = linag::DenseMatrix<T>(x.dim().rows,x.dim().cols);

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
const linag::size linag::DenseMatrix<T>::dim() const{
    return dimension;
}

template<typename T>
linag::Vector<T> linag::DenseMatrix<T>::rowToVector(int row){
    assert(dim().row >= 0 && dim().cols < dim().rows);

    auto res = linag::Vector<T>(dim().cols);

    for (int i = 0; i < res.dim(); ++i) {
        res.at(i) = at(row,i);
    }
    return res;
}

template<typename T>
linag::Vector<T> linag::DenseMatrix<T>::colToVector(int col){
    assert(col >= 0 && col < dim().cols);

    auto res = linag::Vector<T>(dim().rows);

    for (int i = 0; i < res.dim(); ++i) {
        res.at(i) = at(i,col);
    }
    return res;
}

template<typename T>
const T& linag::DenseMatrix<T>::at(int row,int col) const{
    assert(row >= 0 && col >= 0 && row < dim().rows && col < dim().cols);

    return data[row + col*dim().rows];
}

template<typename T>
T& linag::DenseMatrix<T>::at(int row,int col){
    assert(row >= 0 && col >= 0 && row < dim().rows && col < dim().cols);

    return data[row + col*dim().rows];
}

template<typename T>
const linag::DenseMatrix<T> linag::DenseMatrix<T>::transpose() const{
    auto res = linag::DenseMatrix<T>(dim().cols,dim().cols);

    for (int i = 0; i < dim().rows; ++i) {
        for (int j = 0; j < dim().cols; ++j) {
            res.at(j,i) = at(i,j);
        }
    }
    return res;
}

template<typename T>
const linag::DenseMatrix<T> linag::DenseMatrix<T>::inverse() const{
    assert(dim().rows == dim().cols);

    linag::DenseMatrix<T> cpy(*this);
    auto res = linag::DenseMatrix<T>(dim());

    for (int i = 0; i < dim().rows; ++i) {
        for (int j = 0; j < dim().cols; ++j) {
            if(i==j)
                res.at(i,j)=1;
            else
                res.at(i,j)=0;
        }
    }

    //gauss-jordan
    for (int k = 0; k < dim().cols; ++k) {
        T diagValue = cpy.at(k,k);
        for (int i = 0; i < dim().rows; ++i) {
            cpy.at(i,k) /= diagValue;
            res.at(i,k) /= diagValue;
        }
        for (int i = 0; i < dim().rows; ++i) {
            if(i==k)
                continue;
            T rowMult = cpy.at(i,k);
            for (int j = i; j <dim().cols; ++j) {
                cpy.at(i,j) -= rowMult * cpy.at(k,j);
                res.at(i,j) -= rowMult * cpy.at(k,j);
            }
        }
    }

    return res;
}

template<typename T>
const linag::DenseMatrix<T> linag::DenseMatrix<T>::operator-() const{
    return -1* (*this);
}

template <typename T>
linag::DenseMatrix<T>::operator Eigen::MatrixXd() const{
    Eigen::MatrixXd res = Eigen::MatrixXcd::Zero(dim().rows,dim().cols);

    for (int i = 0; i < dim().rows; ++i) {
        for (int j = 0; j < dim().cols; ++j) {
            res(i,j)=at(i,j);
        }
    }
    return res;
}

template <typename T>
linag::DenseMatrix<T> & linag::DenseMatrix<T>::operator=(const linag::DenseMatrix<T> &rhs){
    if(this != &rhs){
        if(dimension!=rhs.dim()) {
            dimension = rhs.dim();
            if(dim().rows*dim().cols > 0)
            {
                free(data);
                data = (T*) realloc (data, rhs.dim().rows * rhs.dim().cols * sizeof(T));
            }
            else
                data = (T*) nullptr;
        }
        //memcpy is a "dumb" function that only copies bytes
        std::memcpy(data,rhs.data,rhs.dim().rows * rhs.dim().cols * sizeof(T));
    }
    return *this;
}

template <typename T>
linag::DenseMatrix<T>::DenseMatrix(const DenseMatrix<T> &rhs){
    dimension = rhs.dim();
    if(dimension.rows*dimension.cols > 0)
    {
        data = (T*) malloc(rhs.dim().rows * rhs.dim().cols * sizeof(T));
        assert(data != nullptr);
        //memcpy is a "dumb" function that only copies bytes
        std::memcpy(data,rhs.data,rhs.dim().rows * rhs.dim().cols * sizeof(T));
    }
    else
        data = (T*) nullptr;
}

template <typename T>
linag::DenseMatrix<T>::~DenseMatrix(){
    if(data!= nullptr)
        free(data);
}

template <typename T>
linag::DenseMatrix<T>::DenseMatrix(linag::size dimension):dimension(dimension){
    if(dim().rows*dim().cols > 0)
    {
        data = (T*) malloc(dim().rows * dim().cols * sizeof(T));
        assert(data != nullptr);
    }
    else
        data = (T*) nullptr;
}

template <typename T>
linag::DenseMatrix<T>::DenseMatrix(int rows,int cols){
    dimension.rows=rows;
    dimension.cols=cols;
    if(rows*cols > 0)
    {
        data = (T*) malloc(rows * cols * sizeof(T));
        assert(data != nullptr);
    }
    else
        data = (T*) nullptr;
}



#endif //AUFGABE1_DENSEMATRIX_H

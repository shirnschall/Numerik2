//
// Created by Sebastian Hirnschall on 17.01.20.
//

#ifndef AUFGABE1_VECTOR_H
#define AUFGABE1_VECTOR_H

#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <cassert>
#include <cstring>

#include "size.h"

namespace linag {
    template<typename T>
    class Vector{
    private:
        Size dimension;
        T *data;
    public:
        ~Vector();

        Vector(const Vector<T> &rhs);
        Vector<T> &operator=(const Vector<T> &);
        explicit Vector(int rows,int cols = 1);
        explicit Vector(const Size dimension);

        Vector(std::initializer_list<T> init);

        const Vector<T> operator-() const;

        T& at(int index);
        const T &at(int index) const;

        const Size dim() const;
        unsigned long length() const;

        double l2norm();

        void zeros();
        void ones();
        void rand();

    };

    template<typename T>
    const Vector<T> operator+(const Vector<T>& x,const Vector<T>& y);
    template<typename T>
    const Vector<T> operator-(const Vector<T>& x,const Vector<T>& y);
    template<typename T>
    T operator*(const Vector<T>& x,const Vector<T>& y);
    template<typename T>
    const Vector<T> operator*(const Vector<T>& x,const T y);
    template<typename T>
    const Vector<T> operator/(const Vector<T>& x,const T y);
    template<typename T>
    const Vector<T> operator*(const T x,const Vector<T>& y);

    template<typename T>
    std::ostream& operator<<(std::ostream& output,const Vector<T>& x);
}

template<typename T>
std::ostream& linag::operator<<(std::ostream& output,const linag::Vector<T>& x){
    for (int i = 0; i < x.dim().rows; ++i) {
        for (int j = 0; j < x.dim().cols; ++j) {
            output << x.at(i+j) << ",\t";
        }
        output << '\n';
    }
    return output;
}

template<typename T>
const linag::Vector<T> linag::operator*(const T x,const linag::Vector<T>& y){
    linag::Vector<T> res(y.dim());
    for (int i = 0; i < y.dim().rows * y.dim().cols; ++i) {
        res.at(i) = x*y.at(i);
    }
    return res;
}


template<typename T>
const linag::Vector<T> linag::operator*(const linag::Vector<T>& x,const T y){
    linag::Vector<T> res(x.dim());
    for (int i = 0; i < x.dim().rows * x.dim().cols; ++i) {
        res.at(i) = x.at(i) * y;
    }
    return res;
}


template<typename T>
T linag::operator*(const linag::Vector<T>& x,const linag::Vector<T>& y){
    assert(x.dim().rows*x.dim().cols == y.dim().rows*y.dim().cols);
    T res = (T)0;
    for (int i = 0; i < x.dim().rows*x.dim().cols; ++i) {
        res+= y.at(i) *x.at(i);
    }
    return res;
}


template<typename T>
const linag::Vector<T> operator/(const linag::Vector<T>& x,const T y){
    linag::Vector<T> res(x);
    for (int i = 0; i < res.dim().cols*res.dim().rows; ++i) {
        res.at(i)/=y;
    }
    return res;
}


template<typename T>
const linag::Vector<T> linag::operator-(const linag::Vector<T>& x,const linag::Vector<T>& y){
    return x + (-y);
}

template<typename T>
const linag::Vector<T> linag::operator+(const linag::Vector<T>& x,const linag::Vector<T>& y){
    assert(x.dim()== y.dim());
    linag::Vector<T> res(x.dim());
    for (int i = 0; i < x.dim().rows*x.dim().cols; ++i) {
            res.at(i) = x.at(i) + y.at(i);
    }
    return res;
}


template <typename T>
void linag::Vector<T>::zeros(){
    for (int i = 0; i < length(); ++i) {
        at(i) = 0;
    }
}
template <typename T>
void linag::Vector<T>::ones(){
    for (int i = 0; i < length(); ++i) {
        at(i) = 1;
    }
}

template <typename T>
void linag::Vector<T>::rand(){
    for (int i = 0; i < dim().rows*dim().cols; ++i) {
        at(i) = (T)std::rand()/RAND_MAX;
    }
}

template<typename T>
double linag::Vector<T>::l2norm(){
    double sum = 0;
    for (int i = 0; i < dim().rows; ++i) {
        for (int j = 0; j < dim().cols; ++j) {
            sum += std::fabs(double((at(i+j)*at(i+j))));
        }
    }
    return std::sqrt(sum);
}

template <typename T>
unsigned long linag::Vector<T>::length() const{
    return dim().rows*dim().cols;
}

template <typename T>
const linag::Size linag::Vector<T>::dim() const{
    return dimension;
}

template <typename T>
const T &linag::Vector<T>::at(int index) const{
    assert(index>=0 && index <dim().rows*dim().cols);
    return data[index];
}

template <typename T>
T &linag::Vector<T>::at(int index){
    assert(index>=0 && index <dim().rows*dim().cols);
    return data[index];
}

template <typename T>
const linag::Vector<T> linag::Vector<T>::operator-() const{
    return (T)-1* (*this);
}

template <typename T>
linag::Vector<T>::Vector(std::initializer_list<T> init){
    dimension.rows = init.size();
    dimension.cols = 1;
    if(dim().rows*dim().cols > 0)
    {
        data = (T*) malloc(dim().rows * dim().cols * sizeof(T));
        assert(data != nullptr);
        //copy
        int i=0;
        for(auto item:init) {
                at(i) = item;
                ++i;
        }
    }
    else
        data = (T*) nullptr;
}

template <typename T>
linag::Vector<T>::Vector(const linag::Size dimension):dimension(dimension){
    assert(dimension.rows == 1 || dimension.cols == 1);
    if(dim().rows*dim().cols > 0)
    {
        data = (T*) malloc(dim().rows * dim().cols * sizeof(T));
        assert(data != nullptr);
    }
    else
        data = (T*) nullptr;
}

template <typename T>
linag::Vector<T>::Vector(int rows,int cols){
    assert(rows == 1 || cols == 1);
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

template <typename T>
linag::Vector<T> &linag::Vector<T>::operator=(const linag::Vector<T> &rhs){
    if(this != &rhs){
        if(dimension!=rhs.dim()) {
            dimension = rhs.dim();
            if(dim().rows*dim().cols > 0)
            {
                if(data == nullptr){
                    data = (T *) malloc(rhs.dim().rows * rhs.dim().cols * sizeof(T));
                    assert(data != nullptr);
                }
                else {
                    data = (T *) realloc(data, rhs.dim().rows * rhs.dim().cols * sizeof(T));
                    assert(data != nullptr);
                }
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
linag::Vector<T>::Vector(const linag::Vector<T> &rhs){
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
linag::Vector<T>::~Vector(){
    if(data!= nullptr)
        free(data);
}


#endif //AUFGABE1_VECTOR_H

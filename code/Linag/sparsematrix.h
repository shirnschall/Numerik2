//
// Created by Sebastian Hirnschall on 17.01.20.
//

#ifndef AUFGABE1_SPARSEMATRIX_H
#define AUFGABE1_SPARSEMATRIX_H

#include "vector.h"
#include "cmath"
#include "iostream"

namespace linag {
    template <typename T>
    class SparseMatrix {
    private:
        Vector<T> v;
        Vector<int> I;
        Vector<int> J;

    public:
        SparseMatrix();
        ~SparseMatrix();
    };


    template<typename T>
    const SparseMatrix<T> operator*(const SparseMatrix<T>& x,const T y);
    template<typename T>
    const SparseMatrix<T> operator*(const T x,const SparseMatrix<T>& y);
    template<typename T>
    const Vector<T> operator*(const Vector<T>& x,const SparseMatrix<T>& y);
    template<typename T>
    const Vector<T> operator*(const SparseMatrix<T>& x,const Vector<T>& y);

    template<typename T>
    std::ostream& operator<<(std::ostream& output,const SparseMatrix<T>& x);
}

template<typename T>
const linag::SparseMatrix<T> linag::operator*(const SparseMatrix<T>& x,const T y){
    linag::SparseMatrix<T> res(x);
    for (int i = 0; i < res.v.length(); ++i) {
        res.v.at(i) *= y;
    }
    return res;
}

template<typename T>
const linag::SparseMatrix<T> linag::operator*(const T x,const SparseMatrix<T>& y){
    return y*x;
}

template<typename T>
const linag::Vector<T> linag::operator*(const Vector<T>& x,const SparseMatrix<T>& y){

}

template <typename T>
linag::SparseMatrix<T>::SparseMatrix():v(2),I(2),J(2){

    std::cout << "constructor"<<std::endl;
}
template <typename T>
linag::SparseMatrix<T>::~SparseMatrix(){
    std::cout << "destructor"<<std::endl;
}

#endif //AUFGABE1_SPARSEMATRIX_H

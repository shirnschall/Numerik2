//
// Created by Sebastian Hirnschall on 17.01.20.
//

#ifndef AUFGABE1_SPARSEMATRIX_H
#define AUFGABE1_SPARSEMATRIX_H

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
    std::ostream& operator<<(std::ostream& output,const DenseMatrix<T>& x);
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

#endif //AUFGABE1_SPARSEMATRIX_H

#ifndef AUFGABE1_SPARSEMATRIX_H
#define AUFGABE1_SPARSEMATRIX_H

#include "cmath"
#include "iostream"
#include "vector.h"

#define THREAD_COUNT 8


namespace linag {
    //say class exists without defining it
    template <typename T> class DenseMatrix;


    template <typename T>
    class SparseMatrix {
    private:
        Vector<T> v;
        Vector<int> I;
        Vector<int> J;

        Size dimension;

    public:
        SparseMatrix()= default;
        ~SparseMatrix() = default;

        explicit SparseMatrix(const DenseMatrix<T>& rhs);
        SparseMatrix<T> & operator=(const DenseMatrix<T>& rhs);

        SparseMatrix(const SparseMatrix<T>& rhs);
        SparseMatrix<T> &operator=(const SparseMatrix<T> &rhs);

        const SparseMatrix<T> operator-() const;

        const Size dim() const;

        char isSymmetric() const;

        const Vector<T>& getV() const{ return v;};
        const Vector<int>& getI() const{ return I;};
        const Vector<int>& getJ() const{ return J;};

        Vector<T> conjugateGradientSolver(linag::Vector<T> b, double tau);
    };


    template<typename T>
    const SparseMatrix<T> operator*(const SparseMatrix<T>& x,const T y);
    template<typename T>
    const SparseMatrix<T> operator*(const T x,const SparseMatrix<T>& y);
    template<typename T>
    const Vector<T> operator*(const Vector<T>& x,const SparseMatrix<T>& y);
    template<typename T>
    const Vector<T> operator*(const SparseMatrix<T>& x,const Vector<T>& y);
    //template <typename T>
    //void mult(linag::Vector<T> &res, const linag::SparseMatrix<T>& x,const linag::Vector<T>& y,int idThread,int numThreads);

}

template<typename T>
const linag::SparseMatrix<T> linag::operator*(const linag::SparseMatrix<T>& x,const T y){
    linag::SparseMatrix<T> res(x);
    for (int i = 0; i < res.v.length(); ++i) {
        res.v.at(i) *= y;
    }
    return res;
}

template<typename T>
const linag::SparseMatrix<T> linag::operator*(const T x,const linag::SparseMatrix<T>& y){
    return y*x;
}

template<typename T>
const linag::Vector<T> linag::operator*(const linag::Vector<T>& x,const linag::SparseMatrix<T>& y){
    std::cout << "undefined" << std::endl;
}

template<typename T>
const linag::Vector<T> linag::operator*(const linag::SparseMatrix<T>& x,const linag::Vector<T>& y){
    assert(x.dim().cols == y.length());
    linag::Vector<T> res(y.length());
    res.zeros();
    for (int i = 0; i < res.length(); ++i) {
        for (int j = x.getI().at(i); j < x.getI().at(i + 1); ++j) {
            res.at(i) += y.at(x.getJ().at(j)) * x.getV().at(j);
        }
    }

    return res;
}

//template <typename T>
//void linag::mult(linag::Vector<T> &res, const linag::SparseMatrix<T>& x,const linag::Vector<T>& y,int idThread,int numThreads){

//}

template <typename T>
linag::Vector<T> linag::SparseMatrix<T>::conjugateGradientSolver(linag::Vector<T> b, double tau){
    assert(tau>0 && dim().rows == b.length());

    linag::Vector<T> r1(dim().rows);
    linag::Vector<T> r2(dim().rows);
    linag::Vector<T> d(dim().rows);
    linag::Vector<T> x(dim().rows);
    linag::Vector<T> z(dim().rows);
    x.rand();
    T alpha;
    T betta;
    unsigned long t = 0;
    r1 = b - (*this)*x;
    d = r1;

    do{
        z = (*this)*d;
        alpha = (r1*r1)/(d*z);
        x = x + alpha*d;
        r2 = r1 - alpha*z;
        betta = (r2*r2)/(r1*r1);
        d = r2 + betta*d;

        r1=r2;
    }while (r2.l2norm()>tau);

    return x;
}

template <typename T>
char linag::SparseMatrix<T>::isSymmetric() const{
    return dim().cols == dim().rows?1:0;
}

template <typename T>
const linag::SparseMatrix<T> linag::SparseMatrix<T>::operator-() const{
    linag::SparseMatrix<T> res(*this);
    res.v = res.v *(-1);
    return res;
}

template<typename T>
const linag::Size linag::SparseMatrix<T>::dim() const{
    return dimension;
}

template <typename T>
linag::SparseMatrix<T>::SparseMatrix(const linag::SparseMatrix<T>& rhs){
    dimension = rhs.dim();
    v = rhs.v;
    I = rhs.I;
    J = rhs.J;
}

template <typename T>
linag::SparseMatrix<T> &linag::SparseMatrix<T>::operator=(const linag::SparseMatrix<T> &rhs){
    if(this != &rhs) {
        dimension = rhs.dim();
        v = rhs.v;
        I = rhs.I;
        J = rhs.J;
    }
    return *this;
}

template <typename T>
linag::SparseMatrix<T>::SparseMatrix(const linag::DenseMatrix<T>& rhs):
I(0),J(0),v(0),dimension(rhs.dim()){
    //calculate array size
    int vc = 0;
    int Ic = rhs.dim().rows+1;
    for (int i = 0; i < rhs.dim().rows; ++i) {//rows
        for (int j = 0; j < rhs.dim().cols; ++j) {//cols
            if(std::fabs(rhs.at(i,j)) > 10e-10){
                ++vc;
            }
        }
    }
    //set array size
    I = linag::Vector<int>(Ic);
    J = linag::Vector<int>(vc);
    v = linag::Vector<T>(vc);

    //convert dense matrix to sparse matrix
    vc=0;
    Ic=-1;
    int Jc=0;
    for (int i = 0; i < rhs.dim().rows; ++i) {//rows
        for (int j = 0; j < rhs.dim().cols; ++j) {//cols
            if(std::fabs(rhs.at(i,j)) > 10e-10){
                if(Ic != i){
                    I.at(++Ic) = vc;
                }
                v.at(vc++) = rhs.at(i,j);
                J.at(Jc++) = j;
            }
        }
    }
    I.at(++Ic) = vc;
}

template <typename T>
linag::SparseMatrix<T> & linag::SparseMatrix<T>::operator=(const linag::DenseMatrix<T>& rhs){
    //calculate array size
    int vc=0;
    int Ic =rhs.dim().rows+1;
    for (int i = 0; i < rhs.dim().rows; ++i) {//rows
        for (int j = 0; j < rhs.dim().cols; ++j) {//cols
            if(std::fabs(rhs.at(i,j)) > 10e-10){
                ++vc;
            }
        }
    }
    //rows/cols
    dimension.rows=rhs.dim().rows;
    dimension.cols=rhs.dim().cols;
    //set array size
    I = linag::Vector<int>(Ic);
    J = linag::Vector<int>(vc);
    v = linag::Vector<T>(vc);

    //convert dense matrix to sparse matrix
    vc=0;
    Ic=0;
    int Jc=0;
    for (int i = 0; i < rhs.dim().rows; ++i) {//rows
        for (int j = 0; j < rhs.dim().cols; ++j) {//cols
            if(std::fabs(rhs.at(i,j)) > 10e-10){
                if(!Ic || Ic != i){
                    I.at(Ic++) = vc;
                }
                v.at(vc++) = rhs.at(i,j);
                J.at(Jc++) = j;
            }
        }
    }
}


#endif //AUFGABE1_SPARSEMATRIX_H

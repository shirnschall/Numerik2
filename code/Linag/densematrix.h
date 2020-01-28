#ifndef AUFGABE1_DENSEMATRIX_H
#define AUFGABE1_DENSEMATRIX_H

//eigen lib
#include "Eigen/Dense"
#include <iostream>
#include <cstring>
#include "size.h"
#include <thread>
#define THREAD_COUNT 8


namespace linag {
    template <typename T> class SparseMatrix;
    template <typename T> class Vector;

template<typename T>
class DenseMatrix{
private:
    Size dimension;
    T* data;
public:
    DenseMatrix(int rows = 0,int cols = 0);
    explicit DenseMatrix(linag::Size dimension);
    ~DenseMatrix();

    DenseMatrix(std::initializer_list<std::initializer_list<T>> init);


    DenseMatrix(const DenseMatrix<T> &rhs);
    DenseMatrix<T> &operator=(const DenseMatrix<T> &rhs);


    explicit DenseMatrix(const SparseMatrix<T>& rhs);
    DenseMatrix<T> &operator=(const SparseMatrix<T> &rhs);



    Eigen::MatrixXd toEigen() const;

    const DenseMatrix<T> operator-() const;

    const DenseMatrix<T> inverse() const;
    const DenseMatrix<T> transpose() const;


    T& at(int row,int col);
    const T& at(int row,int col) const;
    Vector<T> colToVector(int col);
    Vector<T> rowToVector(int row);
    const Size dim() const;

    void zeros();
    void id();
    void diag(T value);
    void diag(const DenseMatrix<T>& rhs);
    void rand();
    //upper tirangular matrix
    void randLT();

    void randDiag();

    //rand sym,pos def
    void randSPD(int notZeroPerLine);

    char isSymmetric() const;

    Vector<T> conjugateGradientSolver(linag::Vector<T> b, double tau, int* count = nullptr,Vector<linag::Vector<double>*>* xs = nullptr, linag::Vector<double>* rs = nullptr);

    double cond();

};

template<typename T>
const DenseMatrix<T> operator+(const DenseMatrix<T>& x,const DenseMatrix<T>& y);
template<typename T>
const DenseMatrix<T> operator-(const DenseMatrix<T>& x,const DenseMatrix<T>& y);
template<typename T>
const DenseMatrix<T> operator*(const DenseMatrix<T>& x,const DenseMatrix<T>& y);
template <typename T>
void mult_DenseMatrix_DenseMatrix(linag::DenseMatrix<T>& res,const DenseMatrix<T>& x,const DenseMatrix<T>& y,int idThread,int numThreads);



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

//template spezialication

    template <>
    void linag::DenseMatrix<double >::rand(){
        for (int i = 0; i < dim().rows; ++i) {
            for (int j = 0; j < dim().cols; ++j) {
                at(i,j) = (double)std::rand()/RAND_MAX;
            }
        }
    }

    template<>
    void linag::DenseMatrix<double>::randLT() {
        for (int i = 0; i < dim().cols; ++i) {
            for (int j = 0; j < i+1; ++j) {
                at(i,j) = (double)std::rand()/RAND_MAX;
            }
            for (int j = i+1; j < dim().rows; ++j) {
                at(i,j) = 0;
            }
        }
    }

    template <>
    void linag::DenseMatrix<double>::randDiag(){
        int n = dim().cols<dim().rows?dim().cols:dim().rows;
        for (int i = 0; i < dim().rows; ++i) {
            for (int j = 0; j < dim().cols; ++j) {
                if(i == j)
                    at(i,j) = (double)std::rand()/RAND_MAX;
                else
                    at(i,j) = 0;
            }
        }
    }

    template<>
    void linag::DenseMatrix<double>::randSPD(int notZeroPerLine) {
        assert(isSymmetric() && notZeroPerLine <= dim().cols && notZeroPerLine%2);

        randLT();
        double c = 50;
        linag::DenseMatrix<double> diagM(dim());
        diagM.randDiag();

        (*this) = (*this) + transpose() + c * diagM;

        for (int i = 0; i < dim().rows; ++i) {
            for (int j = i+std::ceil((double)notZeroPerLine/2); j < dim().cols; ++j) {
                at(i,j) = 0;
            }
            for (int j = 0; j <= i-std::ceil((double)notZeroPerLine/2); ++j) {
                at(i,j) = 0;
            }
        }
    }
}



template<typename T>
std::ostream& linag::operator<<(std::ostream& output,const linag::DenseMatrix<T>& x){
    for (int i = 0; i < x.dim().rows; ++i) {
        for (int j = 0; j < x.dim().cols; ++j) {
            output << x.at(i,j) << ",\t";
        }
        output << '\n';
    }
    return output;
}

template<typename T>
const linag::Vector<T> linag::operator*(const linag::DenseMatrix<T>& x,const linag::Vector<T>& y){
    assert(x.dim().cols == y.dim().cols*y.dim().rows);

    linag::Vector<T> res(x.dim().rows);

    for (int i = 0; i < res.dim().rows*res.dim().cols; ++i) {
        res.at(i) = 0;
        for (int j = 0; j < x.dim().cols; ++j) {
            res.at(i) += x.at(i,j) * y.at(j);
        }
    }
    return res;
}

template<typename T>
const linag::Vector<T> linag::operator*(const linag::Vector<T>& x,const linag::DenseMatrix<T>& y){
    assert(y.dim().rows == x.dim().cols*x.dim().rows);

    linag::Vector<T> res(y.dim().cols);

    for (int i = 0; i < res.dim().rows*res.dim().cols; ++i) {
        res.at(i) = 0;
        for (int j = 0; j < y.dim().rows; ++j) {
            res.at(i) += x.at(j) * y.at(i,j);
        }
    }
    return res;
}

template<typename T>
const linag::DenseMatrix<T> linag::operator*(const T x,const linag::DenseMatrix<T>& y){
    linag::DenseMatrix<T> res(y.dim());
    for (int i = 0; i < y.dim().rows; ++i) {
        for (int j = 0; j < y.dim().cols; ++j) {
            res.at(i,j) = x*y.at(i,j);
        }
    }
    return res;
}

template<typename T>
const linag::DenseMatrix<T> linag::operator*(const linag::DenseMatrix<T>& x,const T y){
    return y*x;
}

template<typename T>
const linag::DenseMatrix<T> linag::operator*(const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y){
    assert(x.dim().cols==y.dim().rows);

    linag::DenseMatrix<T> res(x.dim().rows,y.dim().cols);

    //multithreading
    linag::Vector<std::thread*> threads(THREAD_COUNT>res.dim().cols?res.dim().cols:THREAD_COUNT);
    for (int l = 0; l < threads.length(); ++l) {
        //create threads
        threads.at(l) = new std::thread(linag::mult_DenseMatrix_DenseMatrix<T>,std::ref(res),std::ref(x),std::ref(y),l,threads.length());
    }
    for (int l = 0; l < threads.length(); ++l) {
        threads.at(l)->join();
    }
    //std::thread test(std::thread(linag::mult<T>,std::ref(res),std::ref(x),std::ref(y),0,1));
    //test.join();
    for (int l = 0; l < threads.length(); ++l) {
        delete threads.at(l);
    }

    return res;
}

template <typename T>
void linag::mult_DenseMatrix_DenseMatrix(linag::DenseMatrix<T>& res,const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y,int idThread,int numThreads){
    for (int i = idThread; i < res.dim().cols; i+=numThreads) {
        for (int j = 0; j < x.dim().rows; ++j) {
            res.at(j,i)=0;
            for (int k = 0; k < x.dim().cols; ++k) {
                res.at(j,i) += x.at(j,k) * y.at(k,i);
            }
        }
    }
}


template<typename T>
const linag::DenseMatrix<T> linag::operator-(const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y){
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
const linag::DenseMatrix<T> linag::operator+(const linag::DenseMatrix<T>& x,const linag::DenseMatrix<T>& y){
    return x-(-y);
}

template<typename T>
const linag::Size linag::DenseMatrix<T>::dim() const{
    return dimension;
}

template<typename T>
linag::Vector<T> linag::DenseMatrix<T>::rowToVector(int row){
    assert(dim().row >= 0 && dim().cols < dim().rows);

    linag::Vector<T> res(dim().cols);

    for (int i = 0; i < res.dim(); ++i) {
        res.at(i) = at(row,i);
    }
    return res;
}

template<typename T>
linag::Vector<T> linag::DenseMatrix<T>::colToVector(int col){
    assert(col >= 0 && col < dim().cols);

    linag::Vector<T> res(dim().rows);

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
    linag::DenseMatrix<T> res(dim().cols,dim().cols);

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
    linag::DenseMatrix<T> res(dim());

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
    //int k = 2;
        T diagValue = cpy.at(k,k);
        for (int i = 0; i < dim().cols; ++i) {
            cpy.at(k,i) /= diagValue;
            res.at(k,i) /= diagValue;
        }
        for (int i = 0; i < dim().rows; ++i) {
            if(i==k)
                continue;
            T rowMult = cpy.at(i,k);
            for (int j = 0; j < dim().cols; ++j) {
                cpy.at(i,j) -= rowMult * cpy.at(k,j);
                res.at(i,j) -= rowMult * res.at(k,j);
            }
        }
    }
    //std::cout << cpy << std::endl << res << std::endl;

    return res;
}

template<typename T>
const linag::DenseMatrix<T> linag::DenseMatrix<T>::operator-() const{
    return (T)-1* (*this);
}

template <>
Eigen::MatrixXd linag::DenseMatrix<double>::toEigen () const{
    Eigen::MatrixXd res = Eigen::MatrixXd(dim().rows,dim().cols);

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
                if(!data)
                {
                    data = (T*) malloc (rhs.dim().rows * rhs.dim().cols * sizeof(T));
                    assert(data != nullptr);
                }else {
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
linag::DenseMatrix<T>::DenseMatrix(const DenseMatrix<T> &rhs):dimension(rhs.dim()){
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
linag::DenseMatrix<T>::DenseMatrix(std::initializer_list<std::initializer_list<T>> init){
    dimension.rows = init.size();
    dimension.cols = init.begin()->size();
    //check if all rows have same length
    for(auto row : init){
        assert(dim().cols == row.size());
    }

    if(dim().rows*dim().cols > 0)
    {
        data = (T*) malloc(dim().rows * dim().cols * sizeof(T));
        assert(data != nullptr);
        //copy
        int i=0,j;
        for(auto row:init) {
            j=0;
            for (auto item:row) {
                at(i,j) = item;
                ++j;
            }
            ++i;
        }
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
linag::DenseMatrix<T>::DenseMatrix(linag::Size dimension):dimension(dimension){
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
    dimension.rows = rows;
    dimension.cols = cols;
    if(rows*cols > 0)
    {
        data = (T*) malloc(rows * cols * sizeof(T));
        assert(data != nullptr);
    }
    else
        data = (T*) nullptr;
}

template <typename T>
void linag::DenseMatrix<T>::zeros(){
    for (int i = 0; i < dim().rows; ++i) {
        for (int j = 0; j < dim().cols; ++j) {
            at(i,j) = (T)0;
        }
    }
}

template <typename T>
void linag::DenseMatrix<T>::id(){
    for (int i = 0; i < dim().rows; ++i) {
        for (int j = 0; j < dim().cols; ++j) {
            if(i==j)
                at(i,j) = 1;
            else
                at(i,j) = (T)0;
        }
    }
}

template <typename T>
linag::Vector<T> linag::DenseMatrix<T>::conjugateGradientSolver(linag::Vector<T> b, double tau, int* count,linag::Vector<linag::Vector<double>*>* xs, linag::Vector<double>* rs){
    assert(tau>0 && dim().rows == b.length());
    if(xs)
        assert(xs->length() == dim().rows);//exact result after n iterations
    if(rs)
        assert(rs->length() == dim().rows);//exact result after n iterations

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
    if(xs) {
        for (int i = 1; i < xs->length(); ++i) {
            xs->at(i) = nullptr;
        }
        xs->at(0) = new linag::Vector<double>(x);
    }
    if(rs) {
        for (int i = 1; i < rs->length(); ++i) {
            rs->at(i) = 0;
        }
        rs->at(0) = r1.l2norm();
    }
    do{
        z = (*this)*d;
        alpha = (r1*r1)/(d*z);
        x = x + alpha*d;
        r2 = r1 - alpha*z;
        betta = (r2*r2)/(r1*r1);
        d = r2 + betta*d;

        r1=r2;
        if(xs && t < xs->length())
            xs->at(t) = new linag::Vector<double>(x);
        if(rs && t < rs->length())
            rs->at(t) = r2.l2norm();
        ++t;
    }while (r2.l2norm()>tau);
    if(count)
        *count = t;
    return x;
}

template <typename T>
char linag::DenseMatrix<T>::isSymmetric() const{
    return dim().cols == dim().rows?1:0;
}

template <typename T>
void linag::DenseMatrix<T>::diag(T value){
    int n = dim().cols<dim().rows?dim().cols:dim().rows;
    for (int i = 0; i < n; ++i) {
        at(i,i) = value;
    }
}

template <typename T>
void linag::DenseMatrix<T>::diag(const linag::DenseMatrix<T>& rhs){
    assert(dim() == rhs.dim());
    zeros();
    for (int i = 0; i < rhs.dim().rows; ++i) {
        for (int j = 0; j < rhs.dim().cols; ++j) {
            if(i == j)
                at(i,j) = rhs.at(i,j);
        }
    }
}


template <typename T>
linag::DenseMatrix<T>::DenseMatrix(const linag::SparseMatrix<T>& rhs):dimension(rhs.dim()){
    if(dim().rows*dim().cols > 0)
   {
        data = (T*) malloc(dim().rows * dim().cols * sizeof(T));
        assert(data != nullptr);

        zeros();
        for (int i = 0; i < dim().rows; ++i) {
            for (int j = rhs.getI().at(i); j < rhs.getI().at(i+1); ++j) {
                at(i,rhs.getJ().at(j))=rhs.getV().at(j);
            }
        }
    }
    else
        data = (T*) nullptr;
}

template <typename T>
linag::DenseMatrix<T> &linag::DenseMatrix<T>::operator=(const linag::SparseMatrix<T> &rhs){
            dimension = rhs.dim();
            if (rhs.dim().rows * rhs.dim().cols > 0) {
                if (!data) {
                    data = (T *) malloc(rhs.dim().rows * rhs.dim().cols * sizeof(T));
                    assert(data != nullptr);
                } else {
                    data = (T *) realloc(data, rhs.dim().rows * rhs.dim().cols * sizeof(T));
                    assert(data != nullptr);
                }
            }else
            data = (T *) nullptr;


        zeros();
        for (int i = 0; i < dim().rows; ++i) {
            for (int j = rhs.getI().at(i); j < rhs.getI().at(i + 1); ++j) {
                at(i, rhs.getJ().at(j)) = rhs.getV().at(j);
            }
        }
    return *this;
}

template <typename T>
double linag::DenseMatrix<T>::cond(){
    Eigen::VectorXcd eigenvalues = toEigen().eigenvalues();

    double min=std::fabs(eigenvalues(0).real()),max = std::fabs(eigenvalues(0).real());
    for (int i = 1; i < eigenvalues.size(); ++i) {
        if(std::fabs(eigenvalues(i).real()) > max)
            max = std::fabs(eigenvalues(i).real());

        if(std::fabs(eigenvalues(i).real()) < min)
            min = std::fabs(eigenvalues(i).real());
    }

    return max/min;
}



#endif //AUFGABE1_DENSEMATRIX_H

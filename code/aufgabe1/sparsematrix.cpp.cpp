//
// Created by Sebastian Hirnschall on 30.12.19.
//

#include "sparsematrix.hpp"

#include <string.h>
#include <stdlib.h>
#include <time.h>       /* time */

template <typename T>
linag::SparseMatrix<T>::SparseMatrix(int n,linag::vector<int> notzero){
    //generate
    //rows/cols
    rows=n;
    cols=n;
    v.length = n * sum(notzero);
    v.data = (double*)malloc(v.length * sizeof(double));
    I.length = n+1;
    I.data = (int*)malloc(I.length * sizeof(int));
    J.length = v.length;
    J.data = (int*)malloc(J.length * sizeof(int));
    //generate data



}
template <typename T>
linag::SparseMatrix<T>::~SparseMatrix() {
    free(v.data);
    free(I.data);
    free(J.data);
}
//copy sparsematrix
template <typename T>
linag::SparseMatrix<T>::SparseMatrix(const SparseMatrix &other) {
    if(this!=&other) { //reference cannot be null
        v.length = other.getv().length;
        I.length = other.getI().length;
        J.length = other.getJ().length;
        rows=other.rows;
        cols=other.cols;
        //deep copy
        memcpy(v.data,other.getv().data,v.length * sizeof(double));
        memcpy(I.data,other.getI().data,I.length * sizeof(int));
        memcpy(J.data,other.getJ().data,J.length * sizeof(int));
    }
}
//= sparsematrix
template <typename T>
linag::SparseMatrix<T>& linag::SparseMatrix<T>::operator=(const SparseMatrix<T> &other) {
    if(this!=&other) { //reference cannot be null
        v.length = other.getv().length;
        I.length = other.getI().length;
        J.length = other.getJ().length;
        rows=other.rows;
        cols=other.cols;
        //deep copy
        memcpy(v.data,other.getv().data,v.length * sizeof(double));
        memcpy(I.data,other.getI().data,I.length * sizeof(int));
        memcpy(J.data,other.getJ().data,J.length * sizeof(int));
    }
    return *this;
}
//copy matrix
template <typename T>
linag::SparseMatrix<T>::SparseMatrix(const matrix<T> &other ) {
        //calculate array size
        v.length=0;
        I.length=other.rows+1;
        for (int i = 0; i < other.rows; ++i) {//rows
            for (int j = 0; j < other.cols; ++j) {//cols
                if(other.data[i][j] > 10e-10 || other.data[i][j] < -10e-10){
                    ++v.length;
                }
            }
        }
        //rows/cols
        rows=other.rows;
        cols=other.cols;
        //set array size
        J.length=v.length;
        v.data = (double*)malloc(v.length * sizeof(double));
        I.data = (int*)malloc(I.length * sizeof(int));
        J.data = (int*)malloc(J.length * sizeof(int));

        //convert dense matrix to sparse matrix
        int vc=0;
        int Ic=0;
        int Jc=0;
        for (int i = 0; i < other.rows; ++i) {//rows
            for (int j = 0; j < other.cols; ++j) {//cols
                if(other.data[i][j] > 10e-10 || other.data[i][j] < -10e-10){
                    if(!Ic || Ic != i){
                        I.data[Ic++] = vc;
                    }
                    v.data[vc++] = other.data[i][j];
                    J.data[Jc++] = j;
                }
            }
        }
}

template <typename T>
int linag::SparseMatrix<T>::sum(linag::vector<int> a){
    if(!a.data)
        return 0;
    int sum = a.data[0];
    for(int i =1;i<a.length;++i) {
        sum += a.data[i];
    }
    return sum;
}


template <typename T>
linag::Vector<T> linag::conjugateGradientSolver(linag::DenseMatrix<T> A,linag::Vector<T> b,linag::Vector<T> x,double tau){
    linag::Vector<T> r1(A.dim().rows);
    linag::Vector<T> r2(A.dim().rows);
    linag::Vector<T> d(A.dim().rows);
    linag::Vector<T> res(A.dim().rows);
    T alpha;
    T betta;
    unsigned long t = 0;
    r1 = b - prod(A,x);
    d = r1;
    do{
        alpha = linag::prod(r1,r1)/linag::prod(linag::prod(d,A),d);
        res = res + linag::prod(linag::prod(alpha,A),d);
        r2 = r1 - linag::prod(linag::prod(alpha,A),d);
        betta = linag::prod(r2,r2)/linag::prod(r1,r1);
        d = r2 + linag::prod(betta,d);

        r1=r2;
    }while (r2.norm()>tau);


    return res;
}


template <typename T>
linag::matrix<T> linag::SparseMatrix<T>::todense() {
    linag::matrix<T> M={};
    M.rows = rows;
    M.cols = cols;
    M.data = (T**)malloc(rows* sizeof(T*));
    for (int i = 0; i < rows; ++i) {
        //calloc allocates the memory and sets all values to 0
        M.data[i] = (double*)calloc(cols, sizeof(double));
        for (int j = J.data[i]; j < J.data[i+1]-J.data[i]; ++j) {
            M.data[i][I.data[j]]=v.data[j];
        }
    }
    return M;//use freeMatrix to free M.data
}

template<typename T>
linag::DenseMatrix<T>& linag::DenseMatrix<T>::operator=(const linag::DenseMatrix<T> &){

}

template <typename T>
linag::DenseMatrix<T> linag::genRandomMatrix(int n,int numberOfZerosPerLine){
    assert(numberOfZerosPerLine<n);

    srand (time(NULL));

    linag::DenseMatrix<T> randMatrix(n,n);
    linag::DenseMatrix<T> res(n,n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i+1; ++j) {
            randMatrix->at(i,j) = rand();
        }
        for (int j = i+1; j < n; ++j) {
            randMatrix->at(i,j) = 0;
        }
    }

    randMatrix = linag::prod(randMatrix,linag::trans(randMatrix));
    int index;
    for (int i = 0; i < n; ++i) {   //rows
        int zerosInThisRow = 0;
        for (int l = 0; l <i; ++l) {
            if(std::abs(res.at(i,l))<10e-10)
                ++zerosInThisRow;
        }
        for (int k = 0; k < numberOfZerosPerLine-zerosInThisRow; ++k) {
            do{
                index = (int)floor((i)+(rand()/RAND_MAX)*(n-(i+1)));
            }while(std::abs(res.at(i,index))<10e-10);
            res.at(i,index) = 0;
        }
        for (int j = i+1; j < n; ++j) {
            if(abs(res.at(i,j))<10e-10)
                res.at(i,j) = 0;
        }
    }

    return res;
}

template <typename T>
linag::DenseMatrix<T> linag::trans(linag::DenseMatrix<T> other){
    linag::DenseMatrix<T> res(other.cols,other.rows);
    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            res.at(i,j) = res.at(j,i);
        }
    }
    return res;
}

template <typename T>
double linag::Vector<T>::norm(){
    double sum = 0;
    for(int j = 0; j<l;++j) {
        sum += std::abs(data[j]) * std::abs(data[j]);
    }
    return std::sqrt(sum);
}

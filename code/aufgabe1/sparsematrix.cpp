//
// Created by Sebastian Hirnschall on 30.12.19.
//

#include "sparsematrix.hpp"

#include <string.h>
#include <stdlib.h>


linag::SparseMatrix::SparseMatrix(int n,int* ungleichNull, int ungleichNullc){
    v.length = n * sum(ungleichNull,ungleichNullc);
    v.data = (double*)malloc(v.length * sizeof(double));
    I.length = n+1;
    I.data = (int*)malloc(I.length * sizeof(int));
    J.length = v.length;
    J.data = (int*)malloc(J.length * sizeof(int));
}
linag::SparseMatrix::~SparseMatrix() {
    free(v.data);
    free(I.data);
    free(J.data);
}
//copy sparsematrix
linag::SparseMatrix::SparseMatrix(const SparseMatrix &other) {
    if(this!=&other) { //reference cannot be null
        v.length = other.getv().length;
        I.length = other.getI().length;
        J.length = other.getJ().length;
        //deep copy
        memcpy(v.data,other.getv().data,v.length * sizeof(double));
        memcpy(I.data,other.getI().data,I.length * sizeof(int));
        memcpy(J.data,other.getJ().data,J.length * sizeof(int));
    }
}
//= sparsematrix
linag::SparseMatrix& linag::SparseMatrix::operator=(const SparseMatrix &other) {
    if(this!=&other) { //reference cannot be null
        v.length = other.getv().length;
        I.length = other.getI().length;
        J.length = other.getJ().length;
        //deep copy
        memcpy(v.data,other.getv().data,v.length * sizeof(double));
        memcpy(I.data,other.getI().data,I.length * sizeof(int));
        memcpy(J.data,other.getJ().data,J.length * sizeof(int));
    }
    return *this;
}
//copy matrix
linag::SparseMatrix::SparseMatrix(const double** other,int n,int m) {
    if(other) {//check if other is nullptr
        //calculate array size
        v.length=0;
        I.length=n+1;
        for (int i = 0; i < n; ++i) {//rows
            for (int j = 0; j < m; ++j) {//cols
                if(other[n][m] > 10e-10 || other[n][m] < -10e-10){
                    ++v.length;
                }
            }
        }
        //set array size
        J.length=v.length;
        v.data = (double*)malloc(v.length * sizeof(double));
        I.data = (int*)malloc(I.length * sizeof(int));
        J.data = (int*)malloc(J.length * sizeof(int));

        //convert matrix to sparse matrix
        int vc=0;
        int Ic=0;
        int Jc=0;
        for (int i = 0; i < n; ++i) {//rows
            for (int j = 0; j < m; ++j) {//cols
                if(other[n][m] > 10e-10 || other[n][m] < -10e-10){
                    if(!Ic || Ic != i){
                        I.data[Ic++] = vc;
                    }
                    v.data[vc++] = other[n][m];
                    J.data[Jc++] = j;
                }
            }
        }
    }
}

int linag::SparseMatrix::sum(int* a,int ac){
    if(!a)
        return 0;
    int sum = a[0];
    for(int i =1;i<ac;++i) {
        sum += a[i];
    }
    return sum;
}

double** linag::generateLSData(int n) {
    double** A = (double**)malloc(n* sizeof(double*));
    for(int i=0;i<n;++i){
        A[i] = (double*)malloc(n* sizeof(double));
    }
    //generate data


    return A;
}

linag::vector<double> conjugateGradientSolver(double** data,int n,double* b,double* x,double tau){
    double* r0 = (double*)malloc(n* sizeof(double));
    double* r1 = (double*)malloc(n* sizeof(double));
    double alpha0;
    double alpha1;

}

linag::SparseMatrix* linag::SparseMatrix::transpose() {
    linag::SparseMatrix* trans = new linag::SparseMatrix(*this); //copy this matrix
    //transpose:
    int vc = 0;
    int Ic = 0;
    int Jc=0;
    for (int i = 0; i < trans->size()[0]; ++i) {    //cols
        for (int j = 0; j < trans->getv().length; ++j) {
              if(v.data[j] == i){
                  if(!Ic || Ic != i){
                      trans->setI().data[Ic++] = vc;
                  }
                  trans->setv().data[vc++] = v.data[j];
              }
        }
    }


    return trans;//trans is  pointer. use delete();
}


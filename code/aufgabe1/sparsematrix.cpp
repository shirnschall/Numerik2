//
// Created by Sebastian Hirnschall on 30.12.19.
//

#include "sparsematrix.hpp"

#include <string.h>
#include <stdlib.h>


linag::SparseMatrix::SparseMatrix(int n,linag::vector<int> notzero){
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
        rows=other.rows;
        cols=other.cols;
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
linag::SparseMatrix::SparseMatrix(const matrix<double> &other ) {
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

int linag::SparseMatrix::sum(linag::vector<int> a){
    if(!a.data)
        return 0;
    int sum = a.data[0];
    for(int i =1;i<a.length;++i) {
        sum += a.data[i];
    }
    return sum;
}

double** linag::generateLSData(int n,linag::vector<int> notzero) {
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
    for (int i = 0; i < trans->getI().length; ++i) {
        trans->setI().data[i] = 0;
    }
    for (int i = 0; i < trans->size()[0]; ++i) {    //cols
        for (int j = 0; j < trans->getv().length; ++j) {
              if(v.data[j] == i){
                  //increment number of elements in this col (=rowptr after trans)
                  ++trans->setI().data[i];
                  //calculate new col index



                  //change order of elements in v
                  trans->setv().data[vc++] = v.data[j];
              }
        }
    }


    return trans;//trans is  pointer. use delete();
}

linag::matrix<double> linag::SparseMatrix::todense() {
    linag::matrix<double> M={};
    M.rows = rows;
    M.cols = cols;
    M.data = (double**)malloc(rows* sizeof(double*));
    for (int i = 0; i < rows; ++i) {
        //calloc allocates the memory and sets all values to 0
        M.data[i] = (double*)calloc(cols, sizeof(double));
        for (int j = J.data[i]; j < J.data[i+1]-J.data[i]; ++j) {
            M.data[i][I.data[j]]=v.data[j];
        }
    }
    return M;//use freeMatrix to free M.data
}


void linag::freeMatrix(matrix<double> &M){
    for (int i = 0; i < M.rows; ++i) {
        free(M.data[i]);
    }
    free(M.data);
}

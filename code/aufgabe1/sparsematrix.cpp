//
// Created by Sebastian Hirnschall on 30.12.19.
//

#include "sparsematrix.hpp"

#include <string.h>
#include <stdlib.h>


SparseMatrix::SparseMatrix(int n,int* ungleichNull, int ungleichNullc){
    v.length = n * sum(ungleichNull,ungleichNullc);
    v.data = (double*)malloc(v.length * sizeof(double));
    I.length = n+1;
    I.data = (int*)malloc(I.length * sizeof(int));
    J.length = v.length;
    J.data = (int*)malloc(J.length * sizeof(int));
}
SparseMatrix::~SparseMatrix() {
    free(v.data);
    free(I.data);
    free(J.data);
}
//copy sparsematrix
SparseMatrix::SparseMatrix(const SparseMatrix &other) {
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
SparseMatrix& SparseMatrix::operator=(const SparseMatrix &other) {
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
SparseMatrix::SparseMatrix(const double** other,int n,int m) {
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

int SparseMatrix::sum(int* a,int ac){
    if(!a)
        return 0;
    int sum = a[0];
    for(int i =1;i<ac;++i) {
        sum += a[i];
    }
    return sum;
}
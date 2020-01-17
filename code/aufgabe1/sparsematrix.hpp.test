//
// Created by Sebastian Hirnschall on 30.12.19.
//

#ifndef AUFGABE1_SPARSEMATRIX_HPP
#define AUFGABE1_SPARSEMATRIX_HPP

#include <stdlib.h>
//eigen lib
#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Eigenvalues>



namespace linag {
    template<typename T>
    struct vector {
        int length;
        T *data;
    };
    template<typename T>
    struct matrix {
        int rows,cols;
        T **data;
    };

    struct size {
        int rows,cols;
    };


    template <typename T>
    class Vector{
    private:
        int l;
        T* data;
    public:
        int length(){return l;};
        T& at(int,int);
        const T& at(int,int) const;

        double norm();
    };

    template<typename T>
    class DenseMatrix{
    private:
        int rows,cols;
        T** data;
    public:
        DenseMatrix(int,int);
        ~DenseMatrix();

        DenseMatrix(const DenseMatrix<T> &);
        DenseMatrix &operator=(const DenseMatrix<T> &);
        DenseMatrix(const matrix<T>&);
        DenseMatrix &operator=(const matrix<T> &);

        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> toEigenMatrix();

        T& at(int,int);
        const T& at(int,int) const;
        linag::Vector<T> colToVector(int col);
        linag::Vector<T> rowToVector(int col);
        linag::size dim();
    };


    template<typename T>
    class SparseMatrix {
    private:
        linag::vector<T> v;
        linag::vector<int> J;
        linag::vector<int> I;
        int rows,cols;

        int sum(linag::vector<int>);

    public:
        SparseMatrix(int,linag::vector<int>);

        ~SparseMatrix();

        SparseMatrix &operator=(const SparseMatrix &);

        SparseMatrix(const SparseMatrix &);

        SparseMatrix(const matrix<T>&);

        SparseMatrix &operator=(const matrix<T> &);

        const linag::vector<T> &getv() const { return v; }

        const linag::vector<int> &getI() const { return I; }

        const linag::vector<int> &getJ() const { return J; };

        linag::vector<T> &setv() { return v; }

        linag::vector<int> &setI() { return I; }

        linag::vector<int> &setJ() { return J; };

        linag::matrix<T> todense();

        T& at(int,int);

        int size() const;

    };


    template <typename T>
    linag::DenseMatrix<T> genRandomMatrix(int,int);

    template <typename T>
    Vector<T> conjugateGradientSolver(linag::DenseMatrix<T> A,linag::Vector<T> b,linag::Vector<T> x,double tau);
    template <typename T>
    linag::DenseMatrix<T> prod(linag::DenseMatrix<T>,linag::DenseMatrix<T>);
    template <typename T>
    linag::Vector<T> prod(linag::DenseMatrix<T>,linag::Vector<T>);
    template <typename T>
    linag::Vector<T> prod(linag::Vector<T>,linag::DenseMatrix<T>);
    template <typename T>
    T prod(linag::Vector<T>,linag::Vector<T>);
    template <typename T>
    T prod(T,linag::Vector<T>);
    template <typename T>
    T prod(T,linag::DenseMatrix<T>);

    template <typename T>
    linag::DenseMatrix<T> trans(linag::DenseMatrix<T>);

}
#endif //AUFGABE1_SPARSEMATRIX_HPP

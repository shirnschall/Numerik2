#ifndef RANDOMMATRIX_HPP_INCLUDED
#define RANDOMMATRIX_HPP_INCLUDED

#include <stdlib.h>
#include <cassert>
#include "Vecotr.hpp"

class RandomMatrix
{
    private:
        int n; //dimension
        int k; //nonzero elements
        double * coeff;
    public:
        RandomMatrix(int, int);  //constructor
        ~RandomMatrix();    //destructor
        RandomMatrix(const RandomMatrix&);  //copyconstructor
        RandomMatrix& operator=(const RandomMatrix&);


        int getDimension() const;
        double getCoefficient(int, int) const;
        void printRandomMatrix();

};


#endif // RANDOMMATRIX_HPP_INCLUDED

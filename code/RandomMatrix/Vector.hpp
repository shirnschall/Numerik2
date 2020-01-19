#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>

class Vector
{
	private:
		// dimension of the vector
		int dim;
		// dynamic coefficient vector
		double* coeff;

	public:
		// constructors and destructor
		Vector(int dim, double init = 0);
		//Dreierregel
		Vector(const Vector& rhs);
		~Vector();
		Vector& operator=(const Vector& rhs);

		// return vector dimension
		int size() const;

		// read and write vector coefficients
		void set(int k, double value);
		double get(int k) const;

		//Vector ausgeben
		void printVector() const;
};


#endif // VECTOR_HPP_INCLUDED

#include "Vector.hpp"

Vector::Vector(int dim, double init)
{
	assert(dim > 0);
	this->dim = dim;
	coeff = new double[dim];
	assert(coeff != (double*) 0);

	for (int j = 0; j < dim; j++)
    {
		coeff[j] = init;
	}
}

Vector::Vector(const Vector& rhs)
{
	dim = rhs.size();

	if(dim == 0)
    {
		coeff = (double*) 0;
	}
	else
	{
		coeff = new double[dim];
		for(int i = 0; i < dim; i++)
		{
			coeff[i] = rhs.get(i);
		}
	}
}

Vector::~Vector()
{
	if (dim > 0)
    {
		delete coeff;
	}
}

Vector& Vector::operator=(const Vector& rhs)
{
	if (this != &rhs)
    {
		if (dim != rhs.dim)
		{
			if (dim > 0)
			{
				delete coeff;
			}

			dim = rhs.dim;

			if (dim > 0)
            {
				coeff = new double[dim];
			}
			else
			{
				coeff = (double*) 0;
			}
		}

		for (int j=0; j<dim; ++j)
		{
			coeff[j] = rhs.get(j);
		}
	}
	return *this;
}

int Vector::size() const
{
	return dim;
}

void Vector::set(int k, double value)
{
	assert(k >= 0 && k < dim);
	coeff[k] = value;
}

double Vector::get(int k) const
{
	assert(k >= 0 && k < dim);
	return coeff[k];
}

void Vector::printVector() const
{
	for(int i = 0; i < dim; i++)
    {
		std::cout << coeff[i] << std::endl;
	}
}

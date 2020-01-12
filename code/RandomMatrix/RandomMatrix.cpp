#include "RandomMatrix.hpp"

RandomMatrix::RandomMatrix(int n, int n)
{
    assert(n > 0);
    assert(k <= n && k > 0);
	this->n = n;
	this->k = k
	coeff = new double[n*n];
	assert(coeff != (double*) 0);


}

RandomMatrix::~RandomMatrix()
{
    if(n > 0)
		delete coeff;
}

RandomMatrix::RandomMatrix(const RandomMatrix& A)
{
	n = A.getDimension();

	if(n > 0)
    {
        coeff = new double[n*n];
        assert(coeff != (double*) 0);

		for(int i = 0; i < length; i++)
		{
			coeff[i] = A.coeff[i];
		}
	}
}

RandomMatrix& RandomMatrix::operator=(const RandomMatrix& rhs)
{
	if (this != &rhs)
    {
		if (n != rhs.getDimension())
		{
			if (n > 0)
			{
				delete coeff;
			}
			n = rhs.getDimension();

			if (n > 0)
			{
                coeff = new double[n*n];
                assert(coeff != (double*) 0);

				for (int i = 0; i < length; i++)
				{
					coeff[i] = rhs.coeff[i];
				}
			}
			else
			{
				coeff = (double*) 0;
			}
		}
	}
	return *this;
}

int RandomMatrix::getDimension() const
{
    return n;
}

double getCoefficient(int j, int k) const
{
    return coeff[j*n + k];
}
void RandomMatrix::printRandomMatrix()
{
    	for(int i = 0; i < n; i++)
    {
		for(int j = 0; j < n; j++)
		{
			double a = getCoefficient(i,j);
			if(a >= 0)
				std::cout << " ";

			std::cout << (int)(getCoefficient(i, j)*10000)/10000.0 << "\t";
		}
		std::cout << std::endl;
	}
}

#include <iostream>
#include <iomanip>


void print (
	unsigned nRows,
	unsigned nCols,
	double* flatMatrix
) noexcept
{
	int index=0;
	for( unsigned i=0; i<nRows; ++i)
	{
		for(unsigned j=0; j<nCols; ++j, ++index)
			std::cout << std::setprecision(3) << flatMatrix[index] << " ";
		std::cout<<std::endl;
	}
}


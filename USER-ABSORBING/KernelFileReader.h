#ifndef KERNELFILEREADER_H_INCLUDED
#define KERNELFILEREADER_H_INCLUDED

#include "KernelRelationComponent.h"
#include <fstream>
#include <map>
#include <memory>


class KernelFileReaderImp {
	std::ifstream fin;

public:
	KernelFileReaderImp (
		const char* filename
	) noexcept;


	void
	Read (	double tCut,
		double timestep,
		KernelFunctionsMap&
	) noexcept;

private:
	void
	NumericalIntegralScaling (
		double kernelFunction[],
		double scalingFactor,
		unsigned nData
	) const noexcept;


	void
	Interpolate (
		const double function[],
		unsigned nData,
		double dt,
		double interpolateFunction[],
		unsigned nUsed,
		double timeStep
	) const noexcept;


	double
	Interpolate (
		const double function[],
		unsigned nData,
		double dx,
		double x
	) const noexcept;
};

//---------------------------------------------------------------------------//

template < typename LATTICE >
class KernelFileReader {
public:
	void
	ReadKernelFile (
		double tCut,
		double timestep,
		KernelFunctionsMap& kernelFunction
	) const noexcept
	{
		KernelFileReaderImp _Imp( LATTICE::kernelFile );
		_Imp.Read( tCut, timestep, kernelFunction );
	}
};

#endif // KERNELFILEREADER_H_INCLUDED

#include "PrekernelCalculator.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <omp.h>

/*
  output all the kernel function θ = Zk'

  output file name as :
     [unitC]_[iAtomC]_[coordC]-[unitA]_[jAtomA]_[coordA]
*/

template < typename TOPOLOGY,
	typename STIFF >
struct KernelCalculator {
	MKL_INT dof;
	PrekernelCalculator& prekernels;
	const TOPOLOGY& topology;
	const STIFF& stiff;
#ifdef PARALLEL
	MKL_INT begin;
	MKL_INT end;
	std::vector<int> procSize;
	std::vector<int> shift;
#endif


	KernelCalculator (
		MKL_INT dof,
		PrekernelCalculator& prekernels,
		const TOPOLOGY& topology,
		const STIFF& stiff
	) noexcept :
		dof( dof ),
		prekernels( prekernels ),
		topology( topology ),
		stiff( stiff )
	{
#ifdef PARALLEL
		end = 0;
#endif
	}


	void
	calculate (
		MKL_INT aWidth = 0,
		MKL_INT cWidth = 0
	) noexcept
	{
#ifdef PARALLEL
		MKL_INT myPID;
		MKL_INT nProc;
		blacs_pinfo_(&myPID, &nProc);
#endif
		MKL_INT dofPerUnit = STIFF::DOF;
		auto aRange = topology.range(aWidth, dofPerUnit);
		auto cRange = topology.range(cWidth, dofPerUnit);
		for(auto i=cRange.begin(); i!=cRange.end(); ++i)
		{
			for(auto j=aRange.begin(); j<aRange.end(); ++j)
			{
				auto kernel = localCalculation(*i, *j);
#ifdef PARALLEL
				if(myPID==0)
#endif
					dumpKernel(*i, *j, kernel);
#ifdef PARALLEL
				char all = 'A';
				blacs_barrier_(&context, &all);
#endif
			}
			;
		}
	}


#ifndef PARALLEL
	std::vector<double>
	localCalculation (
		MKL_INT idegree,
		MKL_INT jdegree
	) noexcept
	{
		std::vector<double> kernel;
		MKL_INT dofPerLayer = STIFF::DOF * topology.nUnits;
		for(MKL_INT k=0; k<dofPerLayer; ++k)
		{
			MKL_INT kk = k + dofPerLayer;
			std::cout<<"kk: "<<kk<<"\n";
			MKL_INT jj = jdegree;
			auto func = topology(kk, jj, stiff);
			double factor = (stiff.*func)(kk, jj);
			std::cout<<"factor: "<<factor<<"\n";
			if( fabs(factor) < 1.e-5 ) continue;
			axpy(kernel, prekernels[k+dof*idegree], -factor);
		}
		return kernel;
	}


	void
	axpy (	std::vector<double>& y,
		const std::vector<double>& x,
		double factor
	) const noexcept
	{
		if ( y.size() == 0 )
			y.resize( x.size() ); // nts

		#pragma omp parallel for simd
		for(unsigned i=0; i<y.size(); ++i)
			y[i] += factor * x[i];
	}
#else
	std::vector<double>
	localCalculation (
		MKL_INT idegree,
		MKL_INT jdegree
	) noexcept
	{
		std::vector<double> localKernel;
		MKL_INT dofPerLayer = STIFF::DOF * topology.nUnits;
		for(MKL_INT k=0; k<dofPerLayer; ++k)
		{
			MKL_INT kk = k+dofPerLayer;
			MKL_INT jj = jdegree;
			auto func = topology(kk, jj, stiff);
			double factor = (stiff.*func)(kk, jj);	
			std::cout<<"kk:"<<kk<<"kk size:"<<kk.size()<<"\n";
			if( fabs(factor) < 1.e-5 ) continue;
			localAxpy(localKernel, prekernels[k+dof*idegree], -factor);
		}
		std::vector<double> kernel(shift.back());
		MPI_Gatherv(localKernel.data(), end-begin, MPI_DOUBLE,
			kernel.data(), procSize.data(), shift.data(),
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
		return kernel;
	}


	void
	localAxpy (
		std::vector<double>& y,
		const std::vector<double>& x,
		double factor
	) noexcept
	{
		if ( end == 0 ) mpiSizeInit(prekernels.nts);
		MKL_INT dim = end-begin;
		if ( y.size() == 0 )
			y.resize( dim );
		for(MKL_INT i=0; i<dim; ++i)
			y[i] += factor * x[i];
	}


	void
	mpiSizeInit (
		MKL_INT size
	) noexcept
	{
		MKL_INT myPID;
		MKL_INT nProc;
		blacs_pinfo_(&myPID, &nProc);
		for(MKL_INT i=0; i<nProc; ++i)
		{
			MKL_INT begin = i*size/nProc;
			MKL_INT end = (i+1)*size/nProc;
			procSize.push_back(end-begin);
			shift.push_back(begin);
			if ( i==myPID )
			{
				this->begin = begin;
				this->end = end;
			}
		}
		shift.push_back(size);
	}
#endif

	void
	dumpKernel (
		MKL_INT i,
		MKL_INT j,
		const std::vector<double>& kernel
	) const noexcept
	{
		char filename[30];
		prepareFilename(filename, i, j);
		std::ofstream fout(filename);
		for( const auto& value : kernel )
			fout << std::setprecision(18) << value << std::endl;
	}


	void
	prepareFilename (
		char* filename,
		MKL_INT i,
		MKL_INT j
	) const noexcept
	{
		MKL_INT dofPerUnit = STIFF::DOF;

		MKL_INT iUnit = i/dofPerUnit;
		i %= dofPerUnit;
		char ix = alphabet(i%STIFF::DIM);
		i /= STIFF::DIM;

		MKL_INT jUnit = j/dofPerUnit;
		j %= dofPerUnit;
		char jx = alphabet(j%STIFF::DIM);
		j /= STIFF::DIM;

		sprintf(filename, "%d_%d%c-%d_%d%c",
			(int)iUnit, (int)i, ix,
			(int)jUnit, (int)j, jx
		);
	}


	char
	alphabet (
		MKL_INT i
	) const noexcept
	{
		switch(i)
		{
		case 0: return 'x';
			break;
		case 1: return 'y';
			break;
		default:
			return 'z';
		};
	}
};


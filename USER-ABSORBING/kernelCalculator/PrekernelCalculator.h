#include <cmath>
#include <vector>
#include <map>
#include <omp.h>
#include <cstring>
#include <cstddef>
#ifdef PARALLEL
	#include "mpi.h"
#endif
#include "mkl_lapack.h"

#ifndef PARALLEL
void
eigensolver (
	MKL_INT dim,
	std::vector<double>& matrix,
	std::vector<double>& eigenValues
) noexcept
{
	char type = 'V';
	char uplo = 'U';
	MKL_INT lwork = 1+6*dim+2*dim*dim;
	std::vector<double> work(lwork);
	MKL_INT liwork = 3+5*dim;
	std::vector<MKL_INT> iwork(liwork);
	MKL_INT info;
	dsyevd_(&type, &uplo, &dim, matrix.data(), &dim, eigenValues.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
}

#else

void
eigensolver (
	ParallelMatrix_t& parallelMatrix,
	std::vector<double>& matrix,
	std::vector<double>& eigenValues
) noexcept
{
	MKL_INT nprow;
	MKL_INT npcol;
	MKL_INT myrow;
	MKL_INT mycol;
	blacs_gridinfo_(&context, &nprow, &npcol, &myrow, &mycol);

	MKL_INT dim = parallelMatrix.descM[2];
	auto myMatrixSize = matrix.size();

	char type = 'V';
	char uplo = 'U';

	MKL_INT one = 1;
	std::vector<double> eigenVectors(myMatrixSize);

	MKL_INT lwork = -1;
	std::vector<double> work(1);
	MKL_INT liwork = -1;
	std::vector<MKL_INT> iwork(1);
	MKL_INT info;
	pdsyevd_(&type, &uplo, &dim, matrix.data(), &one, &one, parallelMatrix.descM,
		eigenValues.data(), eigenVectors.data(), &one, &one, parallelMatrix.descM,
		work.data(), &lwork, iwork.data(), &liwork, &info);
	lwork = work[0];
	work.resize(lwork);
	liwork = iwork[0];
	iwork.resize(liwork);
	pdsyevd_(&type, &uplo, &dim, matrix.data(), &one, &one, parallelMatrix.descM,
		eigenValues.data(), eigenVectors.data(), &one, &one, parallelMatrix.descM,
		work.data(), &lwork, iwork.data(), &liwork, &info);

	matrix = eigenVectors;
}
#endif

/*
Z = X * sin(wt)/w * X^T
*/

struct PrekernelCalculator {
	MKL_INT begin;
	MKL_INT end;
	std::vector<double>& evectors;
	std::vector<double>& freqs;
	MKL_INT* descM;
	MKL_INT nts;
	std::vector<std::vector<double>> sinF;
	std::map<MKL_INT, std::vector<double>> savePrekernels;

	PrekernelCalculator (
		std::vector<double>& evectors,
		std::vector<double>& freqs,
		double dt,
		double tmax
	) noexcept :
		evectors( evectors ),
		freqs( freqs ),
		descM( nullptr ),
		nts( tmax/dt+0.001 )
	{
#ifndef PARALLEL 
		begin = 0;
		end = nts;
#else
		MKL_INT myPID;
		MKL_INT nProc;
		blacs_pinfo_(&myPID, &nProc);
		begin = myPID*nts/nProc;
		end = (myPID+1)*nts/nProc;
#endif
		MKL_INT dof = freqs.size();
		sinF.resize( nts );
		for ( MKL_INT t=0; t<nts; ++t )
		{
			sinF[t].resize(dof);
			auto time = (t+1)*dt;
			for ( MKL_INT k=0; k<dof; ++k )
				sinF[t][k] = sin( freqs[k]*time ) / freqs[k];
		}
	}

#ifdef PARALLEL
	PrekernelCalculator (
		std::vector<double>& evectors,
		std::vector<double>& freqs,
		ParallelMatrix_t* parallelMatrix,
		double dt,
		double tmax
	) noexcept :
		PrekernelCalculator(evectors, freqs, dt, tmax)
	{
		descM = parallelMatrix->descM;
	}
#endif


	const std::vector<double>&
	operator[] (
		MKL_INT id
	) noexcept
	{
		if ( savePrekernels.find(id) == savePrekernels.end() )
			CalculatePrekernels( id );
		return savePrekernels.at(id);
	}


	void
	CalculatePrekernels (
		MKL_INT id
	) noexcept
	{
		MKL_INT dof = freqs.size();
		MKL_INT i = id / dof;
		MKL_INT j = id % dof;

		std::vector<double> prekernel;
		MKL_INT t = 0;
		do {
			double theta = thetaF(i,j,t);
			if ( t>=begin && t<end )
				prekernel.push_back(theta);
			++t;
		} while( t<nts );
		savePrekernels[id] = std::move(prekernel);
	}

#ifndef PARALLEL
	double thetaF (
		MKL_INT i, //i=idegree
		MKL_INT j, //j=k
		MKL_INT t
	) noexcept
	{
		MKL_INT dof = freqs.size();
		double theta = 0.;
		for(MKL_INT k=0; k<dof; ++k)
			theta += evectors[i+dof*k] * evectors[j+dof*k] * sinF[t][k];
		return theta;
	}
#else
	double thetaF (
		MKL_INT i,
		MKL_INT j,
		MKL_INT t
	) noexcept
	{
		MKL_INT nprow;
		MKL_INT npcol;
		MKL_INT myrow;
		MKL_INT mycol;
		blacs_gridinfo_(&context, &nprow, &npcol, &myrow, &mycol);

		MKL_INT dof = freqs.size();
		MKL_INT block = descM[4];
		MKL_INT myHeight = descM[8];
		MKL_INT myWidth = evectors.size() / myHeight;
		std::vector<double> copy(myWidth, 0.);
		MKL_INT descCopy[9] = {
			descM[0],
			descM[1],
			nprow,
			descM[3],
			1,
			descM[5],
			0,
			0,
			1
		};
		MKL_INT jCopy = (j/block) % nprow;
		if ( jCopy == myrow )
		{
			MKL_INT jj = j/(block*nprow)*block + (j%block);
			for ( MKL_INT kk=0; kk<myWidth; ++kk )
			{
				MKL_INT index = jj+kk*myHeight;
				MKL_INT k = (kk/block)*(block*npcol) + (kk%block) + mycol*block;
				copy[kk] = evectors[index] * sinF[t][k];
			}
		}

		double localTheta = 0.;
		double maxTheta;
		double minTheta;
		MKL_INT one = 1;
		++i;
		++jCopy;
		char all = 'A';
		blacs_barrier_(&context, &all);
		pddot_(&dof, &localTheta, evectors.data(), &i, &one, descM, &dof,
			copy.data(), &jCopy, &one, descCopy, &nprow);
		MPI_Allreduce(&localTheta, &maxTheta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // pddot output is a loccal value, only and only if processor join this calculation has the exact value
		MPI_Allreduce(&localTheta, &minTheta, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		return fabs(maxTheta) > fabs(minTheta) ? maxTheta : minTheta;
	}
#endif
};


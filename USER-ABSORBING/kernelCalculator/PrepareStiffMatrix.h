#include <vector>
#define MAXBLOCK 100

struct PrepareStiffMatrix {
#ifndef PARALLEL
	template < typename TOPOLOGY,
		typename STIFF >
	std::vector<double>
	operator() (
		const TOPOLOGY& topology,
		MKL_INT nLayers,
		const STIFF& stiff,
		double mass = 9648.50353525
	) const noexcept
	{
		MKL_INT nUnits = topology.nUnits;
		constexpr MKL_INT dofPerUnit = STIFF::DOF;
		MKL_INT dof = dofPerUnit * nUnits * nLayers;
		std::vector<double> flatMatrix(dof * dof);

		double** matrix = new double*[dof];
		for(MKL_INT i=0; i<dof; ++i)
			matrix[i] = &flatMatrix[i*dof];

		//-------------------------------------------------//

		double rmass = 9648.50353525/mass;//at.wt. to ps^2
		for(MKL_INT i=0; i<dof; ++i)
			for(MKL_INT j=0; j<dof; ++j)
			{
				MKL_INT ii = i;
				MKL_INT jj = j;
				auto func = topology(ii, jj, stiff);
				matrix[i][j] = (stiff.*func)(ii, jj) * rmass;
			}

		delete [] matrix;
		return flatMatrix;
	}
#else

	template < typename TOPOLOGY,
		typename STIFF >
	std::vector<double>
	operator() (
		const TOPOLOGY& topology,
		MKL_INT nLayers,
		const STIFF& stiff,
		ParallelMatrix_t& parallelMatrix,
		double mass = 9648.50353525
	) const noexcept
	{
		MKL_INT nprow;
		MKL_INT npcol;
		MKL_INT myrow;
		MKL_INT mycol;
		blacs_gridinfo_(&context, &nprow, &npcol, &myrow, &mycol);

		MKL_INT nUnits = topology.nUnits;
		constexpr MKL_INT dofPerUnit = STIFF::DOF;
		MKL_INT dof = dofPerUnit * nUnits * nLayers;
		MKL_INT block = dofPerUnit * nUnits;
		if ( block > MAXBLOCK ) block = MAXBLOCK;
		MKL_INT nBlock = dof / block;
		MKL_INT nRemain = dof % block;
		MKL_INT myWidth = nBlock/npcol*block + (nBlock%npcol>mycol ? block : 0) + (nBlock%npcol==mycol ? nRemain : 0);
		MKL_INT myHeight = nBlock/nprow*block + (nBlock%nprow>myrow ? block : 0) + (nBlock%nprow==myrow ? nRemain : 0);
		std::vector<double> flatMatrix(myWidth * myHeight);
std::cout<<myrow<<" "<<mycol<<" "<<myWidth<<" "<<myHeight<<std::endl;

		//-------------------------------------------------//

		double rmass = 9648.50353525/mass;//at.wt. to ps^2
		for(MKL_INT j=0; j<myWidth; ++j)
		{
			MKL_INT jjBak = (j/block)*(block*npcol) + (j%block) + mycol*block;
			for(MKL_INT i=0; i<myHeight; ++i)
			{
				MKL_INT jj = jjBak;
				MKL_INT ii = (i/block)*(block*nprow) + (i%block) + myrow*block;
				auto func = topology(ii, jj, stiff);
				flatMatrix[j*myHeight+i] = (stiff.*func)(ii, jj) * rmass; // in fortran
			}
		}

		parallelMatrix.descM[0] = 1;
		parallelMatrix.descM[1] = context;
		parallelMatrix.descM[2] = dof;
		parallelMatrix.descM[3] = dof;
		parallelMatrix.descM[4] = block;
		parallelMatrix.descM[5] = block;
		parallelMatrix.descM[6] = 0;
		parallelMatrix.descM[7] = 0;
		parallelMatrix.descM[8] = myHeight;
		return flatMatrix;
	}
#endif
};


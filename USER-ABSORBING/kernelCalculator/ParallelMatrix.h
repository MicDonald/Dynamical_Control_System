#ifdef PARALLEL
#include "mkl_pblas.h"
#include "mkl_blacs.h"
#include "mkl_scalapack.h"
#include <cmath>

MKL_INT context;

struct ParallelMatrix_t {
	MKL_INT descM[9];
};


void
ParallelInit (
	MKL_INT& myPID,
	MKL_INT& nProc
)
{
	blacs_pinfo_(&myPID, &nProc);
	MKL_INT zero = 0;
	MKL_INT nprow = sqrt(double(nProc));
	MKL_INT npcol = nProc/nprow;
	char row = 'R';
	blacs_get_(&zero, &zero, &context);
	blacs_gridinit_(&context, &row, &nprow, &npcol);
}


void
ParallelFinalize()
{
	MKL_INT zero = 0;
	blacs_gridexit_(&context);
	blacs_exit_(&zero);
}

#endif

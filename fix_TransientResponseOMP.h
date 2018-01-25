#ifdef FIX_CLASS

FixStyle(TransientResponseOMP, FixTransientResponseOMP)

#else

#ifndef FIX_TRANS_RESPONSEOMP
#define FIX_TRANS_RESPONSEOMP

#include "omp.h"
#include "fix_TransientResponse.h"
namespace LAMMPS_NS {

class FixTransientResponseOMP: public FixTransientResponse {
public:
	FixTransientResponseOMP(class LAMMPS *, int, char **);
	void init();
	void initial_integrate(int);
	void final_integrate();
	//void computeForce();
	//void pre_force(int);
	//void post_force(int);

protected:
	int Nthreads;
};

}
#endif
#endif

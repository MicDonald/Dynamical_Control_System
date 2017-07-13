#ifdef FIX_CLASS

FixStyle(VRTransitionOMP, FixVRTransitionOMP)

#else

#ifndef FIX_VIRTUALTRANSITIONOMP
#define FIX_VIRTUALTRANSITIONOMP

#include "omp.h"
#include "fix_VRTransition.h"
namespace LAMMPS_NS {

class FixVRTransitionOMP: public FixVRTransition{
public:	
	FixVRTransitionOMP(class LAMMPS *, int, char **);	
	void init();
	void initial_integrate(int);
	void final_integrate();
	//void computeForce();
	//void pre_force(int);
	//void post_force(int);

protected:
	void recount_topology();
	int Nthreads;
	int mode;
	std::vector<Eigen::SparseMatrix<double>> KM;
	int ifNVE;
};

}
#endif
#endif

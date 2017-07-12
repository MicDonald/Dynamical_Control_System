#ifdef FIX_CLASS

FixStyle(VRTransition, FixVRTransition)

#else

#ifndef FIX_VIRTUALTRANSITION
#define FIX_VIRTUALTRANSITION

#include "fix.h"
#include "fix_nve.h"
#include "atom.h"
#include "force.h"
#include <vector>
#include "KernelMatrix.h"
#include "VRAtomInfo.h"
namespace LAMMPS_NS {

class FixVRTransition: public FixNVE{
public:
	FixVRTransition(class LAMMPS *, int, char **);
	int setmask();
	// virtual void init();
	// virtual void initial_integrate(int);
	// virtual void final_integrate();
	// virtual void computeForce();
	//void pre_force(int);
	//void post_force(int);		
	
protected:
	double test;
	void recount_topology();
	bool m_initialized;
	double minX,maxX,minY,maxY,minZ,maxZ;
	double k_mass=1.;
	KernelMatrix K;
	double t=0,tc=10;
	Eigen::MatrixXd pr;
	Eigen::MatrixXd pv;
	std::vector<Eigen::MatrixXd> ur;
	std::vector<bool> conv;
};

}
#endif
#endif

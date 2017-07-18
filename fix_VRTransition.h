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
#include <cstdlib>
#include "lattice.h"
#include "neighbor.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "bond.h"
#include "atom.h"
#include "molecule.h"
#include "memory.h"
#include "verlet.h"
#include "group.h"
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <memory>
#include <ctime>


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
	void recount_topology();
	int iarg=7;
	double test;
	bool equi_initialized=false,zero_initialized=false;
	double minX,maxX,minY,maxY,minZ,maxZ;
	int Vgroupbit;
	double k_mass=1.;
	KernelMatrix K;
	double t=0,tc=10;
	std::string VRtype;
	Eigen::MatrixXd pr;
	Eigen::MatrixXd pv;
	Eigen::MatrixXd pv_all;
	std::vector<Eigen::MatrixXd> ur;
	std::vector<Eigen::SparseMatrix<double>> KM;
	int mode;
	clock_t t_all=0,t_conv=0,t_nve=0,t_KF=0,t_add=0;
};

}
#endif
#endif

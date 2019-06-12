#ifdef FIX_CLASS

FixStyle(TransientResponse, FixTransientResponse)

#else

#ifndef FIX_TRANS_RESPONSE
#define FIX_TRANS_RESPONSE

#include "fix.h"
#include "fix_nve.h"
#include "atom.h"
#include "force.h"
#include "KernelMatrix.h"
#include "AtomInfo.h"
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
#include <vector>
#include <deque>
#include <cstdlib>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <memory>
#include <ctime>
#include <limits>

namespace LAMMPS_NS {

class FixTransientResponse: public FixNVE {
public:
	FixTransientResponse(class LAMMPS *, int, char **);
	int setmask();
	// virtual void init();
	// virtual void initial_integrate(int);
	// virtual void final_integrate();
	// virtual void computeForce();
	//void pre_force(int);
	//void post_force(int);

protected:
	void recount_topology();
	double test;
	double angle_constant, angle;
	bool equi_initialized = false, zero_initialized = false;
	int Vgroupbit;
	double k_mass = 1.;
	KernelMatrix K;
	double t = 0, tc = 0;
	Eigen::MatrixXd pr;
	Eigen::MatrixXd pv;
	Eigen::MatrixXd pv_all;
	std::deque<Eigen::MatrixXd> ur, ar;
	std::vector<Eigen::MatrixXd> KM;
	char mode;
	clock_t t_all = 0, t_conv = 0, t_nve = 0, t_KF = 0;
	int dtau;
};

}
#endif
#endif

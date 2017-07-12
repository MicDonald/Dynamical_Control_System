#include "fix_VRTransition.h"
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
#include <iterator>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <memory>

using namespace Eigen;
namespace LAMMPS_NS{

//       0        1          2         3  4  5  6  7  8           9                10         11
//fix $fix_ID $group_id vrtransition   x  X  y  Y  z  Z  spring_constant/mass    t cutoff   testMode

FixVRTransition::FixVRTransition (
        class LAMMPS *lmp,
        int narg,
        char **arg
) : FixNVE(lmp, narg, arg),
	minX(atof(arg[3])),maxX(atof(arg[4])),minY(atof(arg[5])),maxY(atof(arg[6])),minZ(atof(arg[7])),maxZ(atof(arg[8])),m_initialized( false ),
	k_mass(atof(arg[9])),tc(atof((arg[10]))),test(atof(arg[11]))
{
  nevery = 1;
  time_integrate = 1;
  K.setK_mass(k_mass);
  std::cout << "k/m = "<<K.getK_mass()<<std::endl;
  std::cout << "Virtual-Real Transition initial with dt = "<<update->dt << std::endl;
  dynamic_group_allow = 1;
}

//---------------------------------------------------------------------------//

int
FixVRTransition::setmask ()
{
  using namespace FixConst;
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
//	mask |= PRE_FORCE;
//	mask |= POST_FORCE;
	mask |= FINAL_INTEGRATE;
  return mask;
}
}
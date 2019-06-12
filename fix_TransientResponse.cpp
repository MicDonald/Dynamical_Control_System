#include "fix_TransientResponse.h"


using namespace Eigen;
using namespace std;
namespace LAMMPS_NS {

//       0        1          2                    3                   4                 5             6        7      8     9
//fix fix_ID ALLgroup_id TransientResponse     VgroupID       spring_const/mass    angle_constant   angle     tc    dtau   mode
//    10
//  Nthreads

FixTransientResponse::FixTransientResponse (
  class LAMMPS *lmp,
  int narg,
  char **arg
) : FixNVE(lmp, narg, arg),
  k_mass(atof(arg[4])), angle_constant(atof(arg[5])), angle(atof(arg[6])), tc(atof(arg[7])), dtau(atof(arg[8])), mode(*arg[9])
{
  if (mode == 'u') cout << "Mode: Displacement" << endl;
  else cout << "Mode: Velocity" << endl;
  if (tc==0) tc = std::numeric_limits<double>::infinity();
  nevery = 1;
  time_integrate = 1;
  K.setK_mass(k_mass);
  std::cout << "k/m = " << K.getK_mass() << std::endl;
  
  int Vgroup = group->find(arg[3]);
  std::cout << "dt = " << update->dt << std::endl;
  if (Vgroup == -1) error->all(FLERR, "Could not find Virtual group ID.");
  Vgroupbit = group->bitmask[Vgroup];
  std::cout << "tc = " << tc << std::endl;
}

//---------------------------------------------------------------------------//

int
FixTransientResponse::setmask ()
{
  using namespace FixConst;
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
//  mask |= PRE_FORCE;
//  mask |= POST_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}


void
FixTransientResponse::recount_topology() {
  bigint nbonds = 0;
  bigint nangles = 0;
  bigint ndihedrals = 0;
  bigint nimpropers = 0;

  if (atom->molecular == 1) {
    int *num_bond = atom->num_bond;
    int *num_angle = atom->num_angle;
    int *num_dihedral = atom->num_dihedral;
    int *num_improper = atom->num_improper;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (num_bond) nbonds += num_bond[i];
      if (num_angle) nangles += num_angle[i];
      if (num_dihedral) ndihedrals += num_dihedral[i];
      if (num_improper) nimpropers += num_improper[i];
    }

  } else if (atom->molecular == 2) {
    Molecule **onemols = atom->avec->onemols;
    int *molindex = atom->molindex;
    int *molatom = atom->molatom;
    int nlocal = atom->nlocal;

    int imol, iatom;

    for (int i = 0; i < nlocal; i++) {
      imol = molindex[i];
      iatom = molatom[i];
      nbonds += onemols[imol]->num_bond[iatom];
      nangles += onemols[imol]->num_angle[iatom];
      ndihedrals += onemols[imol]->num_dihedral[iatom];
      nimpropers += onemols[imol]->num_improper[iatom];
    }
  }

  if (atom->avec->bonds_allow) {
    MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->nbonds /= 2;
  }
  if (atom->avec->angles_allow) {
    MPI_Allreduce(&nangles, &atom->nangles, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->nangles /= 3;
  }
  if (atom->avec->dihedrals_allow) {
    MPI_Allreduce(&ndihedrals, &atom->ndihedrals, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->ndihedrals /= 4;
  }
  if (atom->avec->impropers_allow) {
    MPI_Allreduce(&nimpropers, &atom->nimpropers, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->nimpropers /= 4;
  }
}

// void BondHarmonic::compute(int eflag, int vflag)
// void FixTransientResponseOMP::computeForce(){
//   force->bond->compute(1,1);
// }

}

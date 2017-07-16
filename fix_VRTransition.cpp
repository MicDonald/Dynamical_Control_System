#include "fix_VRTransition.h"


using namespace Eigen;
using namespace std;
namespace LAMMPS_NS{

//       0        1          2                 3                 4            5               6           
//fix fix_ID ALLgroup_id VRtransition  spring_constant/mass    tcutoff  TwoWay/Absorbing  group/length
//group        7   
//          VgroupID
//length       7    8    9   10   11   12
//           minX maxX minY maxY minZ maxZ

FixVRTransition::FixVRTransition (
        class LAMMPS *lmp,
        int narg,
        char **arg
) : FixNVE(lmp, narg, arg),
	k_mass(atof(arg[3])),tc(atof((arg[4]))),mode(atof(arg[5])),VRtype(arg[6])
{
  if (mode==0) cout<<"Simulation Mode: Two Way"<<endl;
  else if (mode==1) cout<<"Simulation Mode: Absorbing"<<endl;
  else error->all(FLERR, "Illegal fix vrtransitionOMP mode");
  nevery = 1;
  time_integrate = 1;
  K.setK_mass(k_mass);
  std::cout << "k/m = "<<K.getK_mass()<<std::endl;
  std::cout << "Virtual-Real Transition initial with dt = "<<update->dt << std::endl;
  if (VRtype=="group"){
    int Vgroup = group->find(arg[iarg++]);
    if (Vgroup == -1) error->all(FLERR,"Could not find Virtual group ID");
    Vgroupbit = group->bitmask[Vgroup];
  }
  else if (VRtype=="length"){
    minX=atof(arg[iarg++]);
    maxX=atof(arg[iarg++]);
    minY=atof(arg[iarg++]);
    maxY=atof(arg[iarg++]);
    minZ=atof(arg[iarg++]);
    maxZ=atof(arg[iarg++]);
  }
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


void
FixVRTransition::recount_topology(){
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

    int imol,iatom;

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
    MPI_Allreduce(&nbonds,&atom->nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (!force->newton_bond) atom->nbonds /= 2;
  }
  if (atom->avec->angles_allow) {
    MPI_Allreduce(&nangles,&atom->nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (!force->newton_bond) atom->nangles /= 3;
  }
  if (atom->avec->dihedrals_allow) {
    MPI_Allreduce(&ndihedrals,&atom->ndihedrals,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (!force->newton_bond) atom->ndihedrals /= 4;
  }
  if (atom->avec->impropers_allow) {
    MPI_Allreduce(&nimpropers,&atom->nimpropers,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (!force->newton_bond) atom->nimpropers /= 4;
  }
}

// void BondHarmonic::compute(int eflag, int vflag)
// void FixVRTransitionOMP::computeForce(){
//   force->bond->compute(1,1);
// }

}

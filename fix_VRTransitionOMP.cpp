#include "fix_VRTransitionOMP.h"
#include <cstdlib>
#include "lattice.h"
#include "neighbor.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "bond.h"
#include "pair.h"
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
using namespace std;
namespace LAMMPS_NS{

//       0        1          2         3  4  5  6  7  8           9               10         11          12           13 
//fix $fix_ID $group_id vrtransition   x  X  y  Y  z  Z  spring_constant/mass   t cutoff   testMode    Nthreads     TwoWay/Absorbing

FixVRTransitionOMP::FixVRTransitionOMP (
        class LAMMPS *lmp,
        int narg,
        char **arg
) :     FixVRTransition(lmp, narg, arg),
        Nthreads(atof(arg[12])),mode(atof(arg[13])),ifNVE(atof(arg[14]))
{       
  if(narg != 15) error->all(FLERR, "Illegal fix vrtransitionOMP command");
  if (mode==0) cout<<"Simulation Mode: Two Way"<<endl;
  else if (mode==1) cout<<"Simulation Mode: Absorbing"<<endl;
  else error->all(FLERR, "Illegal fix vrtransitionOMP mode");
  cout<<"using "<<Nthreads<<" threads"<<endl;
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(Nthreads);
  Eigen::initParallel();
  Eigen::setNbThreads(Nthreads);
  if(ifNVE==1) cout<<"co-NVE"<<endl;
}

//---------------------------------------------------------------------------//
void
FixVRTransitionOMP::init ()
{
#pragma omp single
{
  if (!m_initialized){
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  cout<<"START INIT"<<endl;
  m_initialized = true;
  int nlocal = atom->nlocal;
  int* mask = atom->mask;
  double** x = atom->x;
  double MIN_X=x[nlocal-1][0],MAX_X=x[0][0],MIN_Y=x[nlocal-1][1],MAX_Y=x[0][1],MIN_Z=x[nlocal-1][2],MAX_Z=x[0][2];
  cout<<"loop1 to find local min and MAX"<<endl;
  for ( int i=0; i<nlocal; ++i ){
    if ( mask[i] & groupbit ){
      double tempx=x[i][0],tempy=x[i][1],tempz=x[i][2];
      if (tempx>=MAX_X) MAX_X=tempx;
      if (tempx<=MIN_X) MIN_X=tempx;
      if (tempy>=MAX_Y) MAX_Y=tempy;
      if (tempy<=MIN_Y) MIN_Y=tempy;
      if (tempz>=MAX_Z) MAX_Z=tempz;
      if (tempz<=MIN_Z) MIN_Z=tempz;
    }
  }
  cout<<MIN_X<<" "<<MAX_X<<" "<<MIN_Y<<" "<<MAX_Y<<" "<<MIN_Z<<" "<<MAX_Z<<endl;
  cout<<"loop2 to set model information"<<endl;
  int j=0;
  for ( int i=0; i<nlocal; ++i ){
    if ( mask[i] & groupbit ){
      K.model.atomGID.push_back(atom->tag[i]);
      K.model.atomCoord.push_back({x[i][0],x[i][1],x[i][2]});
      if (x[i][0]>= MIN_X+minX && x[i][0] <= MAX_X-maxX &&
          x[i][1]>= MIN_Y+minY && x[i][1] <= MAX_Y-maxY &&
          x[i][2]>= MIN_Z+minZ && x[i][2] <= MAX_Z-maxZ){
          if (mode==0) K.model.atomVirtual.push_back(j);
          else if(mode==1) K.model.atomReal.push_back(j);
      } else {
        if (mode==0) K.model.atomReal.push_back(j);
		    else if(mode==1) K.model.atomVirtual.push_back(j);
		  }
		++j;
    }
  }
  K.bondNeighborIdendifier();
  cout<<"bonds have been identified"<<endl;
  pr.setZero(K.model.atomR2v.size()*3,1);
  pv.setZero(K.model.atomV2r.size()*3,1);
  cout<<"pos real has been resized"<<endl;
  cout<<"Atoms Global ID in the group: "<<K.model.atomGID.size()<<endl;
  for (auto &e: K.model.atomGID) cout<<e<<" ";
  cout<<endl;
  cout<<"Atoms init-local ID in Real: "<<K.model.atomReal.size()<<endl;
  for (auto &e: K.model.atomReal)    cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
  cout<<endl;

  cout<<"Atoms init-local ID in Virtual: "<<K.model.atomVirtual.size()<<endl;
  for (auto &e: K.model.atomVirtual) cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
  cout<<endl;

  cout<<"Atoms init-local ID in V of VR: "<<K.model.atomV2r.size()<<endl;
  for (auto &e: K.model.atomV2r)     cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
  cout<<endl;

  cout<<"Atoms init-local ID in R of VR: "<<K.model.atomR2v.size()<<endl;
  for (auto &e: K.model.atomR2v)     cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
  cout<<endl;

  for (const auto& e :K.model.atomR2v){
    int gid = K.model.atomGID[e];
    int i = atom->map(gid); // NEED ATOM_STYLE BOND
    vector<int>::iterator iter = find( K.model.atomR2v.begin(), K.model.atomR2v.end(), e );
    int jr = distance( K.model.atomR2v.begin(), iter );
    // cout<<"gid: "<<gid<<" i: "<<i<<" iter: "<<*iter<<" jr: "<<jr<<endl;
    pr(K.bondOrAtom2MatrixDof(jr)[0],0)=x[i][0];
    pr(K.bondOrAtom2MatrixDof(jr)[1],0)=x[i][1];
    pr(K.bondOrAtom2MatrixDof(jr)[2],0)=x[i][2];
  }
  cout<<"Pos[0] of R2v: "<<pr.size()<<"x1\n"<<pr.transpose()<<endl;

  for (const auto& e :K.model.atomV2r){
    int gid = K.model.atomGID[e];
    int i = atom->map(gid);
    vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
    int jv = distance( K.model.atomV2r.begin(), iter );
    pv(K.bondOrAtom2MatrixDof(jv)[0],0)=x[i][0];
    pv(K.bondOrAtom2MatrixDof(jv)[1],0)=x[i][1];
    pv(K.bondOrAtom2MatrixDof(jv)[2],0)=x[i][2];
  }
  cout<<"Pos[0] of V2r: "<<pv.size()<<"x1\n"<<pv.transpose()<<endl;

  int* dlist;
  int n = 0;
  memory->create(dlist,nlocal,"delete_atoms:dlist");
  for (int i = 0; i < nlocal; ++i) dlist[i] = 0;
  for (const auto& e: K.model.atomVirtual)
    if (count(K.model.atomV2r.begin(),
      K.model.atomV2r.end(),e)==0) {
      int gid = K.model.atomGID[e];
      dlist[atom->map(gid)] = 1;
      ++n;
    }

  AtomVec *avec = atom->avec;
  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i,1);
      dlist[i] = dlist[nlocal-1];
      --nlocal;
      } else ++i;
  }
  atom->nlocal = nlocal;
  memory->destroy(dlist);
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
  recount_topology();
  delete dlist;
  //delete avec;
	cout<<"delete "<<n<<" atoms"<<endl;
	K.calculateEigen();
	cout<<"Initialization done"<<endl;
  }
}
}

//---------------------------------------------------------------------------//

void
FixVRTransitionOMP::initial_integrate (int vflag)
{
  double** x = atom->x;
  double** f = atom->f;
  double** v = atom->v;
  int nlocal = atom->nlocal;
  int* mask = atom->mask;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
	int *type = atom->type;
  double dtfm;
  t+=update->dt;
  MatrixXd uv,tempPr,tempVr;
  tempPr.setZero(K.model.atomR2v.size()*3,1);
  tempVr.setZero(K.model.atomR2v.size()*3,1);
  uv.setZero(K.model.atomV2r.size()*3,1);


//Normal NVE for real atom 
#pragma omp single
{
  if(ifNVE==1){
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
    if (rmass) {
      //for (const auto& e :K.model.atomGID){
      for (const auto& e :K.model.atomReal){
        int gid = K.model.atomGID[e];
        //int gid = e;
        int i = atom->map(gid);
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    } else {
      //for (const auto& e :K.model.atomGID){
        for (const auto& e :K.model.atomReal){
          int gid = K.model.atomGID[e];
          //int gid = e;
          int i = atom->map(gid);
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
    }
  }
}

//get ur(t+dt)
#pragma omp single
{
  if (test && t==update->dt) cout<<"i_integrate"<<endl;
  for (const auto& e :K.model.atomR2v){
    int gid = K.model.atomGID[e];
    int i = atom->map(gid);
    vector<int>::iterator iter = find( K.model.atomR2v.begin(), K.model.atomR2v.end(), e );
    int jr = distance( K.model.atomR2v.begin(), iter );
    // cout<<"gid: "<<gid<<" i: "<<i<<" iter: "<<*iter<<" jr: "<<jr<<endl;
    tempPr(K.bondOrAtom2MatrixDof(jr)[0],0)=x[i][0];
    tempPr(K.bondOrAtom2MatrixDof(jr)[1],0)=x[i][1];
    tempPr(K.bondOrAtom2MatrixDof(jr)[2],0)=x[i][2];
  }
  MatrixXd tempUr=tempPr-pr;
  ur.push_back(tempUr);
  KM.push_back(K.calculateKernelMatrix(t));
}
#pragma omp parallel
{
#pragma omp for
  for (int i=t/update->dt;i>=1;--i){
    int ii = t/update->dt-i;
    MatrixXd tempUv=update->dt*KM[i-1]*ur[ii];
    // double ti = i*update->dt;
    // MatrixXd tempUv=update->dt*K.calculateKernelMatrix(ti)*ur[ii];
    #pragma omp critical
    {
      uv+=tempUv;
    }
  }
  
}


//update uv
#pragma omp single
{
  for (const auto& e :K.model.atomV2r){
    int gid = K.model.atomGID[e];
    int i = atom->map(gid);
    vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
    int jv = distance( K.model.atomV2r.begin(), iter );
    x[i][0]=pv(K.bondOrAtom2MatrixDof(jv)[0],0)+uv(K.bondOrAtom2MatrixDof(jv)[0],0);
    x[i][1]=pv(K.bondOrAtom2MatrixDof(jv)[1],0)+uv(K.bondOrAtom2MatrixDof(jv)[1],0);
    x[i][2]=pv(K.bondOrAtom2MatrixDof(jv)[2],0)+uv(K.bondOrAtom2MatrixDof(jv)[2],0);
  }
}

}
//---------------------------------------------------------------------------//

void
FixVRTransitionOMP::final_integrate()
{

  double** f = atom->f;
  double** v = atom->v;
  int nlocal = atom->nlocal;
  int* mask = atom->mask;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;

#pragma omp single
{
  if(ifNVE==1){
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
    double dtfm;

    if (rmass) {
      for (const auto& e :K.model.atomReal){
      int gid = K.model.atomGID[e];
      int i = atom->map(gid);
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }

    } else {
        for (const auto& e :K.model.atomReal){
          int gid = K.model.atomGID[e];
          int i = atom->map(gid);
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
    }
  }
}
}
//---------------------------------------------------------------------------//

void
FixVRTransitionOMP::recount_topology(){
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

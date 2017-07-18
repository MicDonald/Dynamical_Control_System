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

//    8 or 13
//   Nthreads  

FixVRTransitionOMP::FixVRTransitionOMP (
        class LAMMPS *lmp,
        int narg,
        char **arg
) :     FixVRTransition(lmp, narg, arg),
        Nthreads(atof(arg[iarg]))
{      
  if(iarg != 8 && iarg != 13) error->all(FLERR, "Illegal fix vrtransitionOMP command");
  cout<<"using "<<Nthreads<<" threads"<<endl;
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(Nthreads);
  Eigen::initParallel();
  Eigen::setNbThreads(Nthreads);

  if (!equi_initialized){
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
    cout<<"START INIT"<<endl;
    int nlocal = atom->nlocal;
    int* mask = atom->mask;
    double** x = atom->x;
    int j=0;
    if(VRtype=="length"){
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
    }

    else if(VRtype=="group"){
      for ( int i=0; i<nlocal; ++i ){
        if ( mask[i] & groupbit ){
          K.model.atomGID.push_back(atom->tag[i]);
          K.model.atomCoord.push_back({x[i][0],x[i][1],x[i][2]});
          if(mask[i] & Vgroupbit){
            if (mode==0) K.model.atomVirtual.push_back(j);
            else if(mode==1) K.model.atomReal.push_back(j);
          }
          else{
            if (mode==0) K.model.atomReal.push_back(j);
            else if(mode==1) K.model.atomVirtual.push_back(j);
          }
          ++j;
        }
      }
    }

    K.bondNeighborIdendifier();
    cout<<"bonds have been identified"<<endl;
    pr.setZero(K.model.atomR2v.size()*3,1);
    pv.setZero(K.model.atomV2r.size()*3,1);
    pv_all.setZero(K.model.atomVirtual.size()*3,1);
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
    cout<<"Pos[0] of Virtual: "<<pv_all.size()<<"x1\n"<<endl;
    cout<<"Pos[0] of V2r: "<<pv.size()<<"x1\n"<<pv.transpose()<<endl;

    for (const auto& e :K.model.atomVirtual){
      int gid = K.model.atomGID[e];
      int i = atom->map(gid);
      vector<int>::iterator iter = find( K.model.atomVirtual.begin(), K.model.atomVirtual.end(), e );
      int jv = distance( K.model.atomVirtual.begin(), iter );
      pv_all(K.bondOrAtom2MatrixDof(jv)[0],0)=x[i][0];
      pv_all(K.bondOrAtom2MatrixDof(jv)[1],0)=x[i][1];
      pv_all(K.bondOrAtom2MatrixDof(jv)[2],0)=x[i][2];
    }
    K.calculateEigen();
    cout<<"Initialization done"<<endl;
    equi_initialized = true;
  }
}

//---------------------------------------------------------------------------//
void
FixVRTransitionOMP::init ()
{
// zero state
  if (!zero_initialized){
    double** x = atom->x;
    double** v = atom->v;
    int nlocal = atom->nlocal;

    MatrixXd tempPv,tempVv;
    tempPv.setZero(K.model.atomVirtual.size()*3,1);
    tempVv.setZero(K.model.atomVirtual.size()*3,1);
    for (const auto& e :K.model.atomVirtual){
      int gid = K.model.atomGID[e];
      int i = atom->map(gid);
      vector<int>::iterator iter = find( K.model.atomVirtual.begin(), K.model.atomVirtual.end(), e );
      int jv = distance( K.model.atomVirtual.begin(), iter );
      tempPv(K.bondOrAtom2MatrixDof(jv)[0],0)=x[i][0];
      tempPv(K.bondOrAtom2MatrixDof(jv)[1],0)=x[i][1];
      tempPv(K.bondOrAtom2MatrixDof(jv)[2],0)=x[i][2];

      tempVv(K.bondOrAtom2MatrixDof(jv)[0],0)=v[i][0];
      tempVv(K.bondOrAtom2MatrixDof(jv)[1],0)=v[i][1];
      tempVv(K.bondOrAtom2MatrixDof(jv)[2],0)=v[i][2];
    }
    MatrixXd uv0=pv_all-tempPv;
    K.set0state(uv0,tempVv);

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
    zero_initialized=true;
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
  
  MatrixXd uv,tempPr,tempVr;
  tempPr.setZero(K.model.atomR2v.size()*3,1);
  tempVr.setZero(K.model.atomR2v.size()*3,1);
  uv.setZero(K.model.atomV2r.size()*3,1);
  t+=update->dt;
  
  clock_t t_all_temp=clock();
//Normal NVE for real atom
#pragma omp single
{
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
  } 
  else {
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
  clock_t Ttemp=clock();
  t_nve+=Ttemp-t_all_temp;
  
    
//get ur(t+dt)
#pragma omp single
{
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

  clock_t t_temp_KF=clock();
  auto KF=K.calculateKernelMatrix(t);
  t_temp_KF=clock()-t_temp_KF;

  clock_t t_temp_add=clock();
  //KM[t/update->dt]=KF;
  KM.push_back(KF);
  t_temp_add=clock()-t_temp_add;

  t_KF+=t_temp_KF;
  t_add+=t_temp_add;
}
  clock_t t_conv_temp=clock();
#pragma omp parallel
{
#pragma omp for
  for (int i=t/update->dt;i>=1;--i){
    int ii = t/update->dt-i;
    // cout <<"i: " <<i<<endl;
    // cout <<"t: " <<t<<endl;
    // cout <<"I: " <<t/update->dt << endl;
    MatrixXd tempUv=update->dt*KM[i-1]*ur[ii];
    // double ti = i*update->dt;
    // MatrixXd tempUv=update->dt*K.calculateKernelMatrix(ti)*ur[ii];
    #pragma omp critical
    {
      uv+=tempUv;
    }
  } 
}
  t_conv+=clock()-t_conv_temp;

#pragma omp single
{
//update uv
  for (const auto& e :K.model.atomV2r){
    int gid = K.model.atomGID[e];
    int i = atom->map(gid);
    vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
    int jv = distance( K.model.atomV2r.begin(), iter );
    x[i][0]=pv(K.bondOrAtom2MatrixDof(jv)[0],0)+uv(K.bondOrAtom2MatrixDof(jv)[0],0);
    x[i][1]=pv(K.bondOrAtom2MatrixDof(jv)[1],0)+uv(K.bondOrAtom2MatrixDof(jv)[1],0);
    x[i][2]=pv(K.bondOrAtom2MatrixDof(jv)[2],0)+uv(K.bondOrAtom2MatrixDof(jv)[2],0);
  }
  t_all+=clock()-t_all_temp;

  if(int(t/update->dt)%100==0){
    cout<<"\tNVE for real atom        : "<<1.*t_nve/CLOCKS_PER_SEC<<endl;
    cout<<"\tcalcuation of w          : "<<1.*K.t_w/CLOCKS_PER_SEC<<endl;
    cout<<"\tcalcuation of XVX'       : "<<1.*K.t_calK/CLOCKS_PER_SEC<<endl;
    cout<<"\tcalculateKernelMAtrix(t) : "<<1.*t_KF/CLOCKS_PER_SEC<<endl;
    cout<<"\tcollecting Kernel        : "<<1.*t_add/CLOCKS_PER_SEC<<endl;
    cout<<"\ttotal time               : " << 1.*t_all/CLOCKS_PER_SEC << endl;
    cout<<"\tconvolution              : " << 1.*t_conv/CLOCKS_PER_SEC << endl;
    cout<<"\tall but conv             : " << 1.*(t_all-t_conv)/CLOCKS_PER_SEC << endl;
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
  } 
  else {
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
//---------------------------------------------------------------------------//
}


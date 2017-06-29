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
) :     FixNVE(lmp, narg, arg),
	minX(atof(arg[3])),maxX(atof(arg[4])),minY(atof(arg[5])),maxY(atof(arg[6])),minZ(atof(arg[7])),maxZ(atof(arg[8])),
        m_initialized( false ),
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

//---------------------------------------------------------------------------//
void
FixVRTransition::init ()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

        if ( m_initialized ) return;
	std::cout<<"START INIT"<<std::endl;
        m_initialized = true;
	int nlocal = atom->nlocal;
	int* mask = atom->mask;
	double** x = atom->x;
	double MIN_X=x[nlocal-1][0],MAX_X=x[0][0],MIN_Y=x[nlocal-1][1],MAX_Y=x[0][1],MIN_Z=x[nlocal-1][2],MAX_Z=x[0][2];
        std::cout<<"loop1 to find local min and MAX"<<std::endl;
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
	std::cout<<MIN_X<<" "<<MAX_X<<" "<<MIN_Y<<" "<<MAX_Y<<" "<<MIN_Z<<" "<<MAX_Z<<std::endl;
	std::cout<<"loop2 to set model information"<<std::endl;
	int j=0;
	for ( int i=0; i<nlocal; ++i ){
		if ( mask[i] & groupbit ){
		K.model.atomGID.push_back(atom->tag[i]);
		K.model.atomCoord.push_back({x[i][0],x[i][1],x[i][2]});
		if (x[i][0]>= MIN_X+minX && x[i][0] <= MAX_X-maxX &&
		    x[i][1]>= MIN_Y+minY && x[i][1] <= MAX_Y-maxY &&
		    x[i][2]>= MIN_Z+minZ && x[i][2] <= MAX_Z-maxZ) 
			K.model.atomVirtual.push_back(j);
		else
			K.model.atomReal.push_back(j);
		++j;
		}
	}
	K.bondNeighborIdendifier();
    	std::cout<<"bonds have been identified"<<std::endl;
	K.model.atomPos.setZero(K.getgDof(),1);
	pr.setZero(K.getrDof(),1);
       	pv.setZero(K.getvDof(),1);
	std::cout<<"pos real has been resized"<<std::endl;
	K.calculateEigen();
	std::cout<<"eigen value calculation done"<<std::endl;
	std::cout<<"Atoms Global ID in the group: "<<K.model.atomGID.size()<<std::endl;
	for (auto &e: K.model.atomGID) std::cout<<e<<" ";
	std::cout<<std::endl;
	std::cout<<"Atoms init-local ID in Real: "<<K.model.atomReal.size()<<std::endl;
	for (auto &e: K.model.atomReal)    std::cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
        std::cout<<std::endl;

	std::cout<<"Atoms init-local ID in Virtual: "<<K.model.atomVirtual.size()<<std::endl;
	for (auto &e: K.model.atomVirtual) std::cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
        std::cout<<std::endl;

	std::cout<<"Atoms init-local ID in V of VR: "<<K.model.atomV2r.size()<<std::endl;
	for (auto &e: K.model.atomV2r) 	   std::cout<<"["<<e<<"]"<<K.model.atomGID[e]<<"  ";
	std::cout<<std::endl;
	
	for (const auto& e :K.model.atomReal){
                int gid = K.model.atomGID[e];
                int i = atom->map(gid);
                std::vector<int>::iterator iter = std::find( K.model.atomReal.begin(), K.model.atomReal.end(), e );
                int jr = std::distance( K.model.atomReal.begin(), iter );
                pr(K.bondOrAtom2MatrixDof(jr)[0],0)=x[i][0];
                pr(K.bondOrAtom2MatrixDof(jr)[1],0)=x[i][1];
                pr(K.bondOrAtom2MatrixDof(jr)[2],0)=x[i][2];
        }       
        std::cout<<"Pos[0] of real: "<<pr.size()<<"x1\n"<<pr.transpose()<<std::endl;
        
        for (const auto& e :K.model.atomVirtual){
                int gid = K.model.atomGID[e];
                int i = atom->map(gid);
                std::vector<int>::iterator iter = std::find( K.model.atomVirtual.begin(), K.model.atomVirtual.end(), e );
                int jv = std::distance( K.model.atomVirtual.begin(), iter );
                pv(K.bondOrAtom2MatrixDof(jv)[0],0)=x[i][0];
                pv(K.bondOrAtom2MatrixDof(jv)[1],0)=x[i][1];
                pv(K.bondOrAtom2MatrixDof(jv)[2],0)=x[i][2];
        }       
        std::cout<<"Pos[0] of virtual: "<<pv.size()<<"x1\n"<<pv.transpose()<<std::endl;
	j=0;
	for (const auto& gid :K.model.atomGID){
                int i = atom->map(gid);
                K.model.atomPos(K.bondOrAtom2MatrixDof(j)[0],0)=x[i][0];
                K.model.atomPos(K.bondOrAtom2MatrixDof(j)[1],0)=x[i][1];
                K.model.atomPos(K.bondOrAtom2MatrixDof(j)[2],0)=x[i][2];
        	++j;
	}

	int* dlist;
	int n = 0;
	memory->create(dlist,nlocal,"delete_atoms:dlist");
	std::cout<<"ID parsing test:"<<std::endl;
	for (int i = 0; i < nlocal; ++i) dlist[i] = 0;
	for (const auto& e: K.model.atomVirtual)
  		if (count(K.model.atomV2r.begin(),
			  K.model.atomV2r.end(),e)==0) {
			  int gid = K.model.atomGID[e];
			  dlist[atom->map(gid)] = 1;
			  ++n;
		}
	std::cout<<"delete "<<n<<" atoms"<<std::endl;
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
	
	std::cout<<"Initialization done"<<std::endl;

}

//---------------------------------------------------------------------------//

void
FixVRTransition::initial_integrate (int vflag)
{
	if (test)std::cout<<"i_integrate"<<std::endl;
        t+=update->dt;
        MatrixXd tempPr;
	tempPr.setZero(K.getrDof(),1);
	MatrixXd uv;
	uv.setZero(K.getvDof(),1);
        //std::cout<<"tempUr size:"<<tempUr.col(0).size()<<"x"<<tempUr.row(0).size();
        double** x = atom->x;
        double** f = atom->f;
        for (const auto& e :K.model.atomReal){
        	int gid = K.model.atomGID[e];
        	int i = atom->map(gid);
        	std::vector<int>::iterator iter = std::find( K.model.atomReal.begin(), K.model.atomReal.end(), e );
        	int jr = std::distance( K.model.atomReal.begin(), iter );   
        	tempPr(K.bondOrAtom2MatrixDof(jr)[0],0)=x[i][0];
        	tempPr(K.bondOrAtom2MatrixDof(jr)[1],0)=x[i][1];
        	tempPr(K.bondOrAtom2MatrixDof(jr)[2],0)=x[i][2];
        }
        //std::cout<<"pr:"<<pr.transpose()<<std::endl;
        //std::cout<<"tempPr:"<<tempPr.transpose()<<std::endl;
        MatrixXd tempUr=tempPr-pr;
	if (test)
		std::cout<<"tempUr: "<<tempUr.size()<<"x1\n"<<tempUr.transpose()<<std::endl;
	ur.push_back(tempUr);
	std::cout<<"ur size"<<ur.size()<<std::endl;
        int ii=0;
        for (double ti=t;ti>=update->dt;ti-=update->dt){
		if(test) std::cout<<"Kernel("<<ti<<")*ur("<<ii<<")*"<<test<<std::endl;
		uv-=K.calculateKernelMatrix(ti)*ur[ii++];
	}
	//uv*=update->dt;
	uv*=test;
	//for (const auto& e :K.model.atomVirtual){
	for (const auto& e :K.model.atomV2r){
		int gid = K.model.atomGID[e];
                int i = atom->map(gid);
		std::vector<int>::iterator iter = std::find( K.model.atomVirtual.begin(), K.model.atomVirtual.end(), e );
                int jv = std::distance( K.model.atomVirtual.begin(), iter );
		//if (iter!=K.model.atomV2r.end()){
                x[i][0]=pv(K.bondOrAtom2MatrixDof(jv)[0],0)+uv(K.bondOrAtom2MatrixDof(jv)[0],0);
                x[i][1]=pv(K.bondOrAtom2MatrixDof(jv)[1],0)+uv(K.bondOrAtom2MatrixDof(jv)[1],0);
                x[i][2]=pv(K.bondOrAtom2MatrixDof(jv)[2],0)+uv(K.bondOrAtom2MatrixDof(jv)[2],0);
		//}
		//else{
		//	x[i][0]=pv(K.bondOrAtom2MatrixDof(jv)[0],0);
                //      x[i][1]=pv(K.bondOrAtom2MatrixDof(jv)[1],0);
                //        x[i][2]=pv(K.bondOrAtom2MatrixDof(jv)[2],0);
		//}
	}
	if (test)std::cout<<"uv: "<<uv.size()<<"x1\n"<<uv.transpose()<<std::endl;
	computeForce();
        if (test){
	for (const auto& e :K.model.atomReal){
		int gid = K.model.atomGID[e];
                int i = atom->map(gid);
                std::cout<<"RGID: "<<gid<<" fx: "<<f[i][0]<<" fy: "<<f[i][1]<<" fz: "<<f[i][2]<<std::endl;
        }
	for (const auto& e :K.model.atomV2r){
                int gid = K.model.atomGID[e];
                int i = atom->map(gid);
                f[i][0]=f[i][1]=f[i][2]=0.;
		std::cout<<"VGID: "<<gid<<" fx: "<<f[i][0]<<" fy: "<<f[i][1]<<" fz: "<<f[i][2]<<std::endl;
        }

	}

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  double dtfm;

  // update v and x of atoms in group

  double **v = atom->v;
  int nlocal = atom->nlocal;
  int* mask = atom->mask;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
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
//---------------------------------------------------------------------------//

void
FixVRTransition::final_integrate ()
{
  computeForce();
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  
  if (test)std::cout<<"f_integrate"<<std::endl;
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (const auto& e :K.model.atomV2r){
        int gid = K.model.atomGID[e];
        int i = atom->map(gid);
        f[i][0]=f[i][1]=f[i][2]=0.;
  }

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }
}

//---------------------------------------------------------------------------//

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


void FixVRTransition::computeForce(){
	double** f = atom->f;

        for (const auto& e :K.model.atomReal){
                int gid = K.model.atomGID[e];
                int i = atom->map(gid);
                f[i][0]=f[i][1]=f[i][2]=0.;
        }
        //for (const auto& e :K.model.atomVirtual){
        for (const auto& e :K.model.atomV2r){
	        int gid = K.model.atomGID[e];
                int i = atom->map(gid);
                f[i][0]=f[i][1]=f[i][2]=0.;
        }

	force->bond->compute(1,1);

/*	MatrixXd tempU;
        double** x = atom->x;
	tempU.setZero(K.getgDof(),1);
        for (const auto& gid :K.model.atomGID){
                int i = atom->map(gid);
		int j = std::distance(K.model.atomGID.begin(),std::find(K.model.atomGID.begin(),K.model.atomGID.end(),gid));		
                tempU(K.bondOrAtom2MatrixDof(j)[0],0)=x[i][0]-K.model.atomPos(K.bondOrAtom2MatrixDof(j)[0],0);
                tempU(K.bondOrAtom2MatrixDof(j)[1],0)=x[i][1]-K.model.atomPos(K.bondOrAtom2MatrixDof(j)[1],0);
                tempU(K.bondOrAtom2MatrixDof(j)[2],0)=x[i][2]-K.model.atomPos(K.bondOrAtom2MatrixDof(j)[2],0);
	}
	std::cout<<"tempU: "<<tempU.transpose()<<std::endl;
	MatrixXd F=K.D*K.model.atomPos;
	for (const auto& gid :K.model.atomGID){
                int i = atom->map(gid);
                int j = std::distance(K.model.atomGID.begin(),std::find(K.model.atomGID.begin(),K.model.atomGID.end(),gid));
        	f[i][0]=F(K.bondOrAtom2MatrixDof(j)[0],0);
		f[i][1]=F(K.bondOrAtom2MatrixDof(j)[1],0);
		f[i][2]=F(K.bondOrAtom2MatrixDof(j)[2],0);
	}
	std::cout<<"forces computed"<<std::endl;
*/
}


}

#include "fix_TransientResponseOMP.h"

using namespace Eigen;
using namespace std;
namespace LAMMPS_NS {


FixTransientResponseOMP::FixTransientResponseOMP (
  class LAMMPS *lmp,
  int narg,
  char **arg
) : FixTransientResponse(lmp, narg, arg)
{
	if (narg != 9 && narg != 10) error -> all(FLERR, "Illegal fix TransientResponseOMP command");
	if (arg[9])
		Nthreads=atof(arg[9]);
	else
		Nthreads = 8;
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(Nthreads);
	Eigen::initParallel();
	Eigen::setNbThreads(Nthreads);
	
	if (!equi_initialized) {
		dtv = update->dt;
		dtf = 0.5 * update->dt * force->ftm2v;
		int nlocal = atom->nlocal;
		int* mask = atom->mask;
		double** x = atom->x;
		int j = 0;

		for ( int i = 0; i < nlocal; ++i ) {
			if ( mask[i] & groupbit ) {
				K.model.atomGID.push_back(atom->tag[i]);
				K.model.atomCoord.push_back({x[i][0], x[i][1], x[i][2]});
				if (mask[i] & Vgroupbit)
					K.model.atomVirtual.push_back(j);
				else
					K.model.atomReal.push_back(j);

				++j;
			}
		}

		K.bondIdendifier();
		K.angleIdendifier();

		K.printAngles();
		pr.setZero(K.getr2vDof(), 1);
		pv.setZero(K.getv2rDof(), 1);
		pv_all.setZero(K.getvDof(), 1);

		cout << "Atoms Global ID in the group: " << K.model.atomGID.size() << endl;
		for (auto &e : K.model.atomGID) cout << e << " ";
		cout << endl;
		cout << "Atoms init-local ID in Real: " << K.model.atomReal.size() << endl;
		for (auto &e : K.model.atomReal)    cout << "[" << e << "]" << K.model.atomGID[e] << "  ";
		cout << endl;

		cout << "Atoms init-local ID in Virtual: " << K.model.atomVirtual.size() << endl;
		for (auto &e : K.model.atomVirtual) cout << "[" << e << "]" << K.model.atomGID[e] << "  ";
		cout << endl;

		cout << "Atoms init-local ID in V of VR: " << K.model.atomV2r.size() << endl;
		for (auto &e : K.model.atomV2r)     cout << "[" << e << "]" << K.model.atomGID[e] << "  ";
		cout << endl;

		cout << "Atoms init-local ID in R of VR: " << K.model.atomR2v.size() << endl;
		for (auto &e : K.model.atomR2v)     cout << "[" << e << "]" << K.model.atomGID[e] << "  ";
		cout << endl;

		for (const auto& e : K.model.atomR2v) {
			int gid = K.model.atomGID[e];
			int i = atom->map(gid); // NEED ATOM_STYLE BOND
			vector<int>::iterator iter = find( K.model.atomR2v.begin(), K.model.atomR2v.end(), e );
			int jr = distance( K.model.atomR2v.begin(), iter );
			// cout<<"gid: "<<gid<<" i: "<<i<<" iter: "<<*iter<<" jr: "<<jr<<endl;
			pr(K.bondOrAtom2MatrixDof(jr)[0], 0) = x[i][0];
			pr(K.bondOrAtom2MatrixDof(jr)[1], 0) = x[i][1];
			pr(K.bondOrAtom2MatrixDof(jr)[2], 0) = x[i][2];
		}

		for (const auto& e : K.model.atomV2r) {
			int gid = K.model.atomGID[e];
			int i = atom->map(gid);
			vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
			int jv = distance( K.model.atomV2r.begin(), iter );
			pv(K.bondOrAtom2MatrixDof(jv)[0], 0) = x[i][0];
			pv(K.bondOrAtom2MatrixDof(jv)[1], 0) = x[i][1];
			pv(K.bondOrAtom2MatrixDof(jv)[2], 0) = x[i][2];
		}

		for (const auto& e : K.model.atomVirtual) {
			int gid = K.model.atomGID[e];
			int i = atom->map(gid);
			vector<int>::iterator iter = find( K.model.atomVirtual.begin(), K.model.atomVirtual.end(), e );
			int jv = distance( K.model.atomVirtual.begin(), iter );
			pv_all(K.bondOrAtom2MatrixDof(jv)[0], 0) = x[i][0];
			pv_all(K.bondOrAtom2MatrixDof(jv)[1], 0) = x[i][1];
			pv_all(K.bondOrAtom2MatrixDof(jv)[2], 0) = x[i][2];
		}
		K.calculateEigen();
		cout << "Initialization done" << endl;
		equi_initialized = true;
	}
}

//---------------------------------------------------------------------------//
void
FixTransientResponseOMP::init ()
{
// zero state
	if (!zero_initialized) {
		double** x = atom->x;
		double** v = atom->v;
		int nlocal = atom->nlocal;

		MatrixXd tempPv, tempVv;
		tempPv.setZero(K.getvDof(), 1);
		tempVv.setZero(K.getvDof(), 1);
		for (const auto& e : K.model.atomVirtual) {
			int gid = K.model.atomGID[e];
			int i = atom->map(gid);
			vector<int>::iterator iter = find( K.model.atomVirtual.begin(), K.model.atomVirtual.end(), e );
			int jv = distance( K.model.atomVirtual.begin(), iter );
			tempPv(K.bondOrAtom2MatrixDof(jv)[0], 0) = x[i][0];
			tempPv(K.bondOrAtom2MatrixDof(jv)[1], 0) = x[i][1];
			tempPv(K.bondOrAtom2MatrixDof(jv)[2], 0) = x[i][2];

			tempVv(K.bondOrAtom2MatrixDof(jv)[0], 0) = v[i][0];
			tempVv(K.bondOrAtom2MatrixDof(jv)[1], 0) = v[i][1];
			tempVv(K.bondOrAtom2MatrixDof(jv)[2], 0) = v[i][2];
		}
		MatrixXd uv0 = pv_all - tempPv;
		K.set0state(uv0, tempVv);

		int* dlist;
		int n = 0;
		memory->create(dlist, nlocal, "delete_atoms:dlist");
		for (int i = 0; i < nlocal; ++i) dlist[i] = 0;
		for (const auto& e : K.model.atomVirtual)
			if (count(K.model.atomV2r.begin(),
			          K.model.atomV2r.end(), e) == 0) {
				int gid = K.model.atomGID[e];
				dlist[atom->map(gid)] = 1;
				++n;
			}

		AtomVec *avec = atom->avec;
		int i = 0;
		while (i < nlocal) {
			if (dlist[i]) {
				avec->copy(nlocal - 1, i, 1);
				dlist[i] = dlist[nlocal - 1];
				--nlocal;
			} else ++i;
		}
		atom->nlocal = nlocal;
		memory->destroy(dlist);
		bigint nblocal = atom->nlocal;
		MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
		if (atom->map_style) {
			atom->nghost = 0;
			atom->map_init();
			atom->map_set();
		}
		recount_topology();
		delete dlist;
		//delete avec;
		cout << "delete " << n << " atoms" << endl;
		zero_initialized = true;
	}
}

//---------------------------------------------------------------------------//

void
FixTransientResponseOMP::initial_integrate (int vflag)
{
	double** x = atom->x;
	double** f = atom->f;
	double** v = atom->v;
	int nlocal = atom->nlocal ;
	int* mask = atom->mask;
	double *rmass = atom->rmass;
	double *mass = atom->mass;
	int *type = atom->type;
	double dtfm;

	clock_t t_all_temp = clock();

	//Normal NVE for real atom
	#pragma omp single
	{
		if (rmass) {
			// for (const auto& e :K.model.atomGID){
			for (const auto& e : K.model.atomReal) {
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
			if (mode == 'a')
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				//int gid = e;
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
				int jv = distance( K.model.atomV2r.begin(), iter );

				dtfm = dtf / rmass[i];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
				x[i][0] += dtv * v[i][0];
				x[i][1] += dtv * v[i][1];
				x[i][2] += dtv * v[i][2];
			}
			else if (mode == 'x')
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				//int gid = e;
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
				int jv = distance( K.model.atomV2r.begin(), iter );

				dtfm = dtf / rmass[i];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
			}
		}
		else {
			// for (const auto& e :K.model.atomGID){
			for (const auto& e : K.model.atomReal) {
				int gid = K.model.atomGID[e];
				// int gid = e;
				int i = atom->map(gid);
				dtfm = dtf / mass[type[i]];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
				x[i][0] += dtv * v[i][0];
				x[i][1] += dtv * v[i][1];
				x[i][2] += dtv * v[i][2];
			}
			if (mode == 'a')
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				// int gid = e;
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
				int jv = distance( K.model.atomV2r.begin(), iter );

				dtfm = dtf / mass[type[i]];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
				x[i][0] += dtv * v[i][0];
				x[i][1] += dtv * v[i][1];
				x[i][2] += dtv * v[i][2];
			}
			else if (mode == 'x')
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				// int gid = e;
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
				int jv = distance( K.model.atomV2r.begin(), iter );

				dtfm = dtf / mass[type[i]];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
			}
		}
	}


	clock_t Ttemp = clock();
	t_nve += Ttemp - t_all_temp;


	clock_t t_temp_KF = clock();
	if(t < tc)
	{
		auto KF = 0.5 * K.calculateKernelMatrix(t, false)	
		        + 0.5 * K.calculateKernelMatrix(t + 0.5 * dtv, false);
		KM.push_back(KF);
	}
	t_temp_KF = clock() - t_temp_KF;

	t_KF += t_temp_KF;

	if (mode == 'u'|| mode == 'x') {
		MatrixXd uv, tempPr;
		tempPr.setZero(K.getr2vDof(), 1);
		uv.setZero(K.getv2rDof(), 1);
		#pragma omp single
		{
			for (const auto& e : K.model.atomR2v) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomR2v.begin(), K.model.atomR2v.end(), e );
				int jr = distance( K.model.atomR2v.begin(), iter );
				tempPr(K.bondOrAtom2MatrixDof(jr)[0], 0) = x[i][0];
				tempPr(K.bondOrAtom2MatrixDof(jr)[1], 0) = x[i][1];
				tempPr(K.bondOrAtom2MatrixDof(jr)[2], 0) = x[i][2];
			}
			MatrixXd tempUr = tempPr - pr;
			ur.push_front(tempUr);
			
			if(t >= tc) ur.pop_back();
		}

		clock_t t_conv_temp = clock();

		// #pragma omp parallel
		// {
		// 	#pragma omp for
		// 	for (int vi = 0; vi < K.getv2rDof(); ++vi)
		// 		for (int i = 0; i < KM.size(); ++i)
		// 			uv.row(vi) += dtv * KM[i].row(vi) * ur[i].col(0);
		// }
		#pragma omp parallel
		{
			#pragma omp for
			for (int vi = 0; vi < K.getv2rDof(); ++vi)
				for (int i = 0; i < KM.size(); ++i)
					if (abs(i - int(KM.size()))%dtau == 0)
						uv.row(vi) += dtau * dtv * KM[i].row(vi) * ur[i].col(0);
		}
		t_conv += clock() - t_conv_temp;

		#pragma omp single
		{
			// update uv
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
				int jv = distance( K.model.atomV2r.begin(), iter );
				x[i][0] = pv(K.bondOrAtom2MatrixDof(jv)[0], 0) + uv(K.bondOrAtom2MatrixDof(jv)[0], 0);
				x[i][1] = pv(K.bondOrAtom2MatrixDof(jv)[1], 0) + uv(K.bondOrAtom2MatrixDof(jv)[1], 0);
				x[i][2] = pv(K.bondOrAtom2MatrixDof(jv)[2], 0) + uv(K.bondOrAtom2MatrixDof(jv)[2], 0);
			}
			t_all += clock() - t_all_temp;
		}
	}

}
//---------------------------------------------------------------------------//

void
FixTransientResponseOMP::final_integrate()
{

	double** f = atom->f;
	double** v = atom->v;
	int nlocal = atom->nlocal;
	int* mask = atom->mask;
	double *rmass = atom->rmass;
	double *mass = atom->mass;
	int *type = atom->type;

	clock_t t_all_temp = clock();

	#pragma omp single
	{
		double dtfm;
		if (rmass) {
			for (const auto& e : K.model.atomReal) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				dtfm = dtf / rmass[i];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
			}
		}
		else {
			for (const auto& e : K.model.atomReal) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				dtfm = dtf / mass[type[i]];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
			}
		}
	}
// 
	if (mode == 'a'|| mode == 'x') {
		MatrixXd av, tempAr;
		av.setZero(K.getv2rDof(), 1);
		tempAr.setZero(K.getr2vDof(), 1);
		#pragma omp single
		{
			for (const auto& e : K.model.atomR2v) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomR2v.begin(), K.model.atomR2v.end(), e );
				int jr = distance( K.model.atomR2v.begin(), iter );
				tempAr(K.bondOrAtom2MatrixDof(jr)[0], 0) = f[i][0] / mass[type[i]];
				tempAr(K.bondOrAtom2MatrixDof(jr)[1], 0) = f[i][1] / mass[type[i]];
				tempAr(K.bondOrAtom2MatrixDof(jr)[2], 0) = f[i][2] / mass[type[i]];
			}
			ar.push_front(tempAr);

			if(t >= tc) ar.pop_back();

		}
		clock_t t_conv_temp = clock();
		#pragma omp parallel
		{
			#pragma omp for
			for (int vi = 0; vi < K.getv2rDof(); ++vi)
				for (int i = 0; i < KM.size(); ++i)
					av.row(vi) += dtv * KM[i].row(vi) * ar[i].col(0);
		}
		t_conv += clock() - t_conv_temp;

		#pragma omp single
		{
			//update av
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				vector<int>::iterator iter = find( K.model.atomV2r.begin(), K.model.atomV2r.end(), e );
				int jv = distance( K.model.atomV2r.begin(), iter );
				f[i][0] = mass[type[i]] * av(K.bondOrAtom2MatrixDof(jv)[0], 0);
				f[i][1] = mass[type[i]] * av(K.bondOrAtom2MatrixDof(jv)[1], 0);
				f[i][2] = mass[type[i]] * av(K.bondOrAtom2MatrixDof(jv)[2], 0);
			}
		 t_all += clock() - t_all_temp;
		}
	}

	#pragma omp single
	{
		double dtfm;
		if (rmass) {
			if (mode == 'a'|| mode == 'x')
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				dtfm = dtf / rmass[i];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
			}
		}
		else {
			if (mode == 'a'|| mode == 'x')
			for (const auto& e : K.model.atomV2r) {
				int gid = K.model.atomGID[e];
				int i = atom->map(gid);
				dtfm = dtf / mass[type[i]];
				v[i][0] += dtfm * f[i][0];
				v[i][1] += dtfm * f[i][1];
				v[i][2] += dtfm * f[i][2];
			}
		}
	}


	if(t == 0 || t/dtv == 100){
		cout << "\tNVE for real atom        : " << 1.* t_nve/CLOCKS_PER_SEC << endl;
		cout << "\tcalcuation of w          : " << 1.* K.t_w/CLOCKS_PER_SEC << endl;
		cout << "\tcalcuation of XV, X'DVR  : " << 1.* K.t_calK/CLOCKS_PER_SEC << endl;
		cout << "\tcalculateKernelMatrix(t) : " << 1.* t_KF/CLOCKS_PER_SEC << endl;
		cout << "\ttotal time               : " << 1.* t_all/CLOCKS_PER_SEC << endl;
		cout << "\tconvolution              : " << 1.* t_conv/CLOCKS_PER_SEC << endl;
		cout << "\tall but conv             : " << 1.* (t_all-t_conv)/CLOCKS_PER_SEC << endl;
	}

	
	t += dtv;
}
//---------------------------------------------------------------------------//
}

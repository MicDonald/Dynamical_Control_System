/*
  F = -KU - K'U'
  K = X * Λ * X^T
  w = sqrt(Λ)
  Z = X^T * sin(w)t/w * X
  kernel function θ = -Z * K'
*/
#include "mkl_types.h"
#include <iostream>
#include <fstream>

#include "PotentialPack.h"

#include "ParallelMatrix.h"
#include "UnitCellStiff.h"
#include "PrepareStiffMatrix.h"
#include "KernelCalculator.h"
#include "print.h"

#include "Lattice_Triangular.h"
#include "Lattice_Triangular2.h"
#include "Lattice_Triangular3.h"
#include "Lattice_Hexagonal.h"
#include "Lattice_FCC.h"
#include "Lattice_Diamond.h"

#include "Topology_2D_Face.h"
#include "Topology_2D_Edge.h"
#include "Topology_3D_Face.h"
#include "Topology_3D_Edge.h"
#include "Topology_3D_Corner.h"


using namespace std;


constexpr MKL_INT width = 11;
constexpr MKL_INT nLayers = 50;
constexpr MKL_INT aWidth = 5;
constexpr MKL_INT cWidth = 0;
constexpr double dt =   0.01; // ps
constexpr double tmax = 100.0; // ps MAX:100
constexpr MKL_INT N = 1;
#define POTENTIAL potential(0.010340808, 3.4) //Ar
//#define POTENTIAL potential(0.2381, 3.4050)
//#define POTENTIAL potential


#ifndef PARALLEL
void DumpMatrix(const std::vector<double>& matrix, unsigned dof)
{
	ofstream fout("matrix.txt");
	for(unsigned i=0; i<matrix.size(); ++i)
	{
		fout << matrix[i] << " ";
		if ( i%dof == dof-1 )
			fout << "\n";
	}
	fout.flush();
}

template < typename POT, typename LATTICE, typename TOPOLOGY >
void SerialVersion (double latticeSpacing, double cutoff, double mass)
{
	using STIFF = UnitCellStiff<LATTICE, POT>;
	constexpr MKL_INT nUnits = TOPOLOGY::numOfUnits_Of( width ); // per layer
	constexpr MKL_INT dof = nLayers*nUnits*LATTICE::N*LATTICE::DIM;

	POT POTENTIAL;
	LATTICE lattice(latticeSpacing);

	STIFF stiff(lattice, potential, cutoff);
	TOPOLOGY topology(nUnits);
	auto latticeMatrix = PrepareStiffMatrix() (topology, nLayers, stiff, mass);
	DumpMatrix(latticeMatrix, dof);

	vector<double> eigenvalues(dof);
	eigensolver( dof, latticeMatrix, eigenvalues );

	for( double& eValue : eigenvalues )
	{
//cout << eValue<<" ";
		eValue = sqrt( eValue );
	}
//cout << endl;

	vector<double>& w = eigenvalues;
	vector<double>& eValues = latticeMatrix;
	PrekernelCalculator prekernel(eValues, w, dt, tmax);
	KernelCalculator<TOPOLOGY, STIFF> kernelCalculator(dof, prekernel, topology, stiff);
	kernelCalculator.calculate(aWidth, cWidth, N);
}
#else

template < typename POT, typename LATTICE, typename TOPOLOGY >
void ParallelVersion (double latticeSpacing, double cutoff, double mass)
{
	using STIFF = UnitCellStiff<LATTICE, POT>;
	MKL_INT myProc;
	MKL_INT nProcs;
	ParallelInit(myProc, nProcs);

	MKL_INT nprow;
	MKL_INT npcol;
	MKL_INT myrow;
	MKL_INT mycol;
	blacs_gridinfo_(&context, &nprow, &npcol, &myrow, &mycol);
	cout <<myProc<<"/"<<nProcs<<" "<< nprow<<" "<<npcol << " " <<myrow << " " <<mycol<<endl;

	constexpr MKL_INT nUnits = TOPOLOGY::numOfUnits_Of( width );
	POT POTENTIAL;
	LATTICE lattice(latticeSpacing);
	STIFF stiff(lattice, potential, cutoff);
	TOPOLOGY topology(nUnits);
	ParallelMatrix_t parallelMatrix;
cout<<"creating matrix"<<endl;
	auto latticeMatrix = PrepareStiffMatrix() (topology, nLayers, stiff, parallelMatrix, mass);

cout<<"solve eigen problem"<<endl;
	MKL_INT dof = nUnits * nLayers * lattice.N * lattice.DIM;
	vector<double> eigenvalues(dof);
	eigensolver( parallelMatrix, latticeMatrix, eigenvalues );
	for( double& eValue:eigenvalues )
		eValue = sqrt( eValue );

cout<<"calculate kernel functions"<<endl;
	auto& w = eigenvalues;
	std::vector<double>& eVector = latticeMatrix;
	PrekernelCalculator prekernel(eVector, w, &parallelMatrix, dt, tmax);
	KernelCalculator<TOPOLOGY, STIFF> kernelCalculator(dof, prekernel, topology, stiff);
	kernelCalculator.calculate(aWidth, cWidth);

	ParallelFinalize();
}
#endif

int main ()
{
//	using POT = Potential_StillingerWeber;
//	using POT = Potential_HarmonicBond;
//	using POT = Potential_HarmonicAngle;
	using POT = Potential_LennardJones;
//	using POT = Potential_EAM;
//	using POT = Potential_Exp6;
//	using POT = Potential_Tersoff;

//	using LATTICE = Lattice_Triangular3;
//	using LATTICE = Lattice_Hexagonal;
	using LATTICE = Lattice_FCC;
//	using LATTICE = Lattice_Diamond;

//	using TOPOLOGY = Topology_2D_Face;
//	using TOPOLOGY = Topology_2D_Edge;
	using TOPOLOGY = Topology_3D_Face;
//	using TOPOLOGY = Topology_3D_Edge;
//	using TOPOLOGY = Topology_3D_Corner;

//	double latticeSpacing = 1.;
//	double latticeSpacing = 4.08;
//	double latticeSpacing = 5.4309; // Si diamond
//	double latticeSpacing = 5.4309*0.25*sqrt(3.); // Si hexagonal
	double latticeSpacing = 5.256; // Ar FCC

//	double cutoff = 1.01;
//	double cutoff = sqrt(3.)*0.25*5.454; // Si
//	double cutoff = 8.5; // Ar LJ 2.5*sigma
	double cutoff = 6.8; // Ar LJ 2*sigma

//	double mass = 9648.50353525;
//	double mass = 196.97; // Au
//	double mass = 12.0107; // C
//	double mass = 28.0855; // Si
	double mass = 39.948; // Ar

#ifndef PARALLEL
	SerialVersion<POT, LATTICE, TOPOLOGY>(latticeSpacing, cutoff, mass);
#else
	ParallelVersion<POT, LATTICE, TOPOLOGY>(latticeSpacing, cutoff, mass);
#endif
}


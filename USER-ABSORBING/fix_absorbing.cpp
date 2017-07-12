#include "fix_absorbing.h"
#include "lattice.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "atom.h"
#include <iterator>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include "mpi.h"

namespace LAMMPS_NS {

LAMMPS_Difference::LAMMPS_Difference (
	const LAMMPS* lammpsPtr
) noexcept :
	lammps ( lammpsPtr )
{}


void
LAMMPS_Difference::SetLatticeSpacing (
	double latticeSpacing
) noexcept
{
	this->latticeSpacing = latticeSpacing;
	rLatticeSpacing = 1./latticeSpacing;
}


void
LAMMPS_Difference::Correction (
	double* rij
) const noexcept
{
	rij[0] *= latticeSpacing;
	rij[1] *= latticeSpacing;
	rij[2] *= latticeSpacing;
	lammps->domain->minimum_image( rij );
	rij[0] *= rLatticeSpacing;
	rij[1] *= rLatticeSpacing;
	rij[2] *= rLatticeSpacing;
}

//---------------------------------------------------------------------------//

//       0        1          2        3            4              5           6             7              8
//fix $fix_ID $group_id absorbing $latticeType  $lattice_spacing $time_factor $time_cutoff $space_cutoff $corner_cutoff
FixAbsorbing::FixAbsorbing (
	class LAMMPS *lmp,
	int narg,
	char **arg
) :	Fix(lmp, narg, arg),
	m_differenceFunc( lmp ),
	m_initialized( false )
{
	nevery = 1;
	time_integrate = 1;
	if(narg < 9) error->all(FLERR, "Illegal fix absorbing command");

	m_latticeSpacing = force->numeric(FLERR, arg[4]);
	m_differenceFunc.SetLatticeSpacing( m_latticeSpacing );
	double timeFactor = force->numeric(FLERR, arg[5]);

	double tc, rc, rCorner;
	if(strcmp(arg[6], "INF")==0)
		tc = -1.;
	else
		tc = force->numeric(FLERR, arg[6]);
	if(strcmp(arg[7], "INF")==0)
		rc = -1.;
	else
		rc = force->numeric(FLERR, arg[7]);
	if(strcmp(arg[8], "INF")==0)
		rCorner = -1.;
	else
		rCorner = force->numeric(FLERR, arg[8]);

	m_abc.Setup( arg[3], update->dt*timeFactor, tc*timeFactor, rc/m_latticeSpacing, rCorner/m_latticeSpacing );
std::cout << "absorbing boundary condition initial with dt = "<<update->dt << std::endl;
}

//---------------------------------------------------------------------------//

int
FixAbsorbing::setmask ()
{
	using namespace FixConst;
	int mask = 0;
	mask |= FINAL_INTEGRATE;
	mask |= POST_FORCE;
	return mask;
}

//---------------------------------------------------------------------------//

void
FixAbsorbing::init ()
{
	if ( m_initialized ) return;

	m_initialized = true;
	std::vector<int> myLocalOuterIDs = IdentifyOuterAtoms();
	IdentifyInnerAtoms( myLocalOuterIDs );

	std::vector<data3d> outerPos = AtomPositionsOf(m_bcGIDs);
	std::vector<data3d> innerPos = AtomPositionsOf(m_innerGIDs);
	std::cout << "outer pos:\n";
	for ( const auto& pos : outerPos )
		std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
	std::cout << "\ninner pos:\n";
	for ( const auto& pos : innerPos )
		std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
	std::cout << std::endl;
	m_abc.SetupInitialState(outerPos, innerPos, m_latticeSpacing, &m_differenceFunc);
}

//---------------------------------------------------------------------------//

void
FixAbsorbing::final_integrate ()
{
	std::vector<data3d> innerPos = AtomPositionsOf( m_innerGIDs );
	const std::vector<data3d>& currentOuterPos = CallABC( innerPos );

	int nlocal = atom->nlocal;
	int* mask = atom->mask;
	double** x = atom->x;
	for ( int i=0; i<nlocal; ++i )
		if ( mask[i] & groupbit )
		{
			int gid = atom->tag[i];
			std::vector<int>::iterator iter = std::find( m_bcGIDs.begin(), m_bcGIDs.end(), gid );
			int jndex = std::distance( m_bcGIDs.begin(), iter );
			x[i][0] = currentOuterPos[jndex][0];
			x[i][1] = currentOuterPos[jndex][1];
			x[i][2] = currentOuterPos[jndex][2];
		}
}

//---------------------------------------------------------------------------//

void
FixAbsorbing::post_force (int)
{
	double **f = atom->f;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	for ( int i=0; i<nlocal; ++i )
		if ( mask[i] & groupbit )
			f[i][0] = f[i][1] = f[i][2] = 0.;
}

//---------------------------------------------------------------------------//

std::vector<int>
FixAbsorbing::IdentifyOuterAtoms ()
{
	int nlocal = atom->nlocal;
	int *mask = atom->mask;
	std::vector<int> my_localIDs;
	std::vector<int> my_GIDs;
	for ( int i=0; i<nlocal; ++i )
	{
		if(! (mask[i] & groupbit) ) continue;
		my_localIDs.push_back( i );
		my_GIDs.push_back( atom->tag[i] );
	}
	int mySize = my_localIDs.size();

	// size and shift
	int nProcs;
	int myProc;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
	std::vector<int> processorSize( nProcs );
	MPI_Allgather(&mySize, 1, MPI_INT, &processorSize[0], 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<int> shift( nProcs+1, 0 );
	for ( int i=0; i<nProcs; ++i )
		shift[i+1] = shift[i] + processorSize[i];
	int nBCAtoms = shift[nProcs];

	// gid, mpi gather
	m_bcGIDs.resize( nBCAtoms );
	MPI_Allgatherv(&my_GIDs[0], mySize, MPI_INT, &m_bcGIDs[0], &processorSize[0], &shift[0], MPI_INT, MPI_COMM_WORLD);

	return my_localIDs;
}

//---------------------------------------------------------------------------//

void
FixAbsorbing::IdentifyInnerAtoms ( const std::vector<int>& my_localOuterIDs )
{
	int nlocal = atom->nlocal;// modify using neighborList
	double **x = atom->x;
	int *mask = atom->mask;
	std::vector<int> my_GIDs;
	double latticeSpacingSq = m_latticeSpacing * m_latticeSpacing;
	for ( int i=0; i<nlocal; ++i )
	{
		if( mask[i] & groupbit ) continue; // BC Atom
		for ( int j=0; j<my_localOuterIDs.size(); ++j )
		{
			double dx = x[ my_localOuterIDs[j] ][0] - x[i][0];
			double dy = x[ my_localOuterIDs[j] ][1] - x[i][1];
			double dz = x[ my_localOuterIDs[j] ][2] - x[i][2];
			domain->minimum_image(dx,dy,dz);
			double rSq = dx*dx + dy*dy + dz*dz;
			if(rSq/rSq> (latticeSpacingSq+0.01)/rSq) continue;
			my_GIDs.push_back( atom->tag[i] );
			break;
		}
	}
	int mySize = my_GIDs.size();

	// size and shift
	int nProcs;
	int myProc;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
	std::vector<int> processorSize( nProcs );
	MPI_Allgather(&mySize, 1, MPI_INT, &processorSize[0], 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<int> shift( nProcs+1, 0 );
	for ( int i=0; i<nProcs; ++i )
		shift[i+1] = shift[i] + processorSize[i];

	// gid, mpi gather
	int nNearbyAtoms = shift[nProcs];
	m_innerGIDs.resize( nNearbyAtoms );
	MPI_Allgatherv(&my_GIDs[0], mySize, MPI_INT, &m_innerGIDs[0], &processorSize[0], &shift[0], MPI_INT, MPI_COMM_WORLD);
}


std::vector<data3d>
FixAbsorbing::AtomPositionsOf ( const std::vector<int>& gids )
{
	int nlocal = atom->nlocal;
	double **x = atom->x;
	std::vector<int> my_indexInGIDs;
	std::vector<double> my_Pos;
	for ( int i=0; i<nlocal; ++i )
	{
		int gid = atom->tag[i];
		std::vector<int>::const_iterator iter = std::find( gids.begin(), gids.end(), gid );
		if ( iter == gids.end() ) continue;
		my_indexInGIDs.push_back( std::distance(gids.begin(), iter) );
		my_Pos.push_back( x[i][0] );
		my_Pos.push_back( x[i][1] );
		my_Pos.push_back( x[i][2] );
	}
	int mySize = my_indexInGIDs.size();

	// shift and size
	int nProcs;
	int myProc;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
	std::vector<int> processorSize( nProcs );
	MPI_Allgather(&mySize, 1, MPI_INT, &processorSize[0], 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<int> shift( nProcs, 0 );
	for ( int i=0; i<nProcs-1; ++i )
		shift[i+1] = shift[i] + processorSize[i];

	// index, mpi gather
	std::vector<int> indexInGIDs( gids.size() );
	MPI_Allgatherv(&my_indexInGIDs[0], mySize, MPI_INT, &indexInGIDs[0], &processorSize[0], &shift[0], MPI_INT, MPI_COMM_WORLD);

	// pos, mpi gather
	std::vector<double> tempPos( gids.size()*3 );
	for ( int i=0; i<nProcs; ++i )
	{
		processorSize[i] *= 3;
		shift[i] *= 3;
	}
	MPI_Allgatherv(&my_Pos[0], mySize*3, MPI_DOUBLE, &tempPos[0], &processorSize[0], &shift[0], MPI_DOUBLE, MPI_COMM_WORLD);

	std::vector<data3d> pos( gids.size() );
	for ( unsigned i=0; i<gids.size(); ++i )
	{
		int j = i*3;
		pos[ indexInGIDs[i] ] = {tempPos[j], tempPos[j+1], tempPos[j+2]};
	}
	return pos;
}


const std::vector<data3d>& FixAbsorbing::CallABC (
	const std::vector<data3d>& innerPos
) noexcept
{
	return m_abc.Update( innerPos, Serial_t() );
}

} // namespace LAMMPS


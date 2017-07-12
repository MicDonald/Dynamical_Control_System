#ifdef FIX_CLASS

FixStyle(absorbing/mpi, FixAbsorbingMPI)

#else

#ifndef LMP_FIX_ABSORBINGMPI_H
#define LMP_FIX_ABSORBINGMPI_H

#include "fix_absorbing.h"

namespace LAMMPS_NS {

class FixAbsorbingMPI: public FixAbsorbing {
public:
	FixAbsorbingMPI(class LAMMPS*, int, char** );

protected:
	const std::vector<data3d>& CallABC( const std::vector<data3d>& innerPos ) noexcept override;
};

class FixAbsorbingOMP: public FixAbsorbing {
public:
	FixAbsorbingOMP(class LAMMPS*, int, char** );

protected:
	const std::vector<data3d>& CallABC( const std::vector<data3d>& innerPos ) noexcept override;
};

}

#endif
#endif

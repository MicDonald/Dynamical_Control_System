#ifdef FIX_CLASS

FixStyle(absorbing/omp, FixAbsorbingOMP)

#else

#ifndef LMP_FIX_ABSORBINGOMP_H
#define LMP_FIX_ABSORBINGOMP_H

#include "fix_absorbing.h"

namespace LAMMPS_NS {

class FixAbsorbingOMP: public FixAbsorbing {
public:
	FixAbsorbingOMP(class LAMMPS*, int, char** );

protected:
	const std::vector<data3d>& CallABC( const std::vector<data3d>& innerPos ) noexcept override;
};

}

#endif
#endif

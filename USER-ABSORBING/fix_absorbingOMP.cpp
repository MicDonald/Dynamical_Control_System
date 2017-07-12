#include "fix_absorbingOMP.h"

namespace LAMMPS_NS {

FixAbsorbingOMP::FixAbsorbingOMP (
	class LAMMPS* lmp,
	int narg,
	char** arg
) : FixAbsorbing(lmp, narg, arg )
{}


const std::vector<data3d>& FixAbsorbingOMP::CallABC (
	const std::vector<data3d>& innerPos
) noexcept
{
	return m_abc.Update( innerPos, OpenMP_t() );
}


} // namespace LAMMPS


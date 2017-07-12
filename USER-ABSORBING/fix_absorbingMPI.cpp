#include "fix_absorbingMPI.h"

namespace LAMMPS_NS {

FixAbsorbingMPI::FixAbsorbingMPI (
	class LAMMPS *lmp,
	int narg,
	char **arg
) : FixAbsorbing(lmp, narg, arg)
{}


const std::vector<data3d>& FixAbsorbingMPI::CallABC (
	const std::vector<data3d>& innerPos
) noexcept
{
	return m_abc.Update( innerPos, MPI_t() );
}

} // namespace LAMMPS


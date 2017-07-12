#ifdef FIX_CLASS

FixStyle(absorbing, FixAbsorbing)

#else

#ifndef LMP_FIX_ABSORBING_H
#define LMP_FIX_ABSORBING_H

#include <vector>
#include "fix.h"
#include "math_vector.h"
#include "AbsorbingBoundaryCondition.h"
#include "PBCDifference.h"

namespace LAMMPS_NS {

class LAMMPS_Difference: public PBCDifference {
	const LAMMPS* lammps;
	double latticeSpacing;
	double rLatticeSpacing;

public:
	LAMMPS_Difference ( const LAMMPS* lammps ) noexcept;

	void SetLatticeSpacing ( double latticeSpacing ) noexcept;

	void Correction( double* rij ) const noexcept  override;
};


class FixAbsorbing: public Fix {
public:
	FixAbsorbing(class LAMMPS *, int, char **);
	int setmask();
	void init();

	void final_integrate();
	void post_force(int);

protected:
	std::vector<int> IdentifyOuterAtoms ();
	void IdentifyInnerAtoms ( const std::vector<int>& my_localIDs );
	std::vector<data3d> AtomPositionsOf ( const std::vector<int>& gids );
	virtual const std::vector<data3d>& CallABC( const std::vector<data3d>& innerPos ) noexcept;

	LAMMPS_Difference m_differenceFunc;
	AbsorbingBoundaryCondition m_abc;
	std::vector<int> m_bcGIDs;
	std::vector<int> m_innerGIDs;
	double m_latticeSpacing;
	bool m_initialized;
};

}

#endif
#endif

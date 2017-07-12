#ifndef ABSORBINGBOUNDARYCONDITION_H_INCLUDED
#define ABSORBINGBOUNDARYCONDITION_H_INCLUDED

#include "data3d.h"
#include "Lattice.h"
#include "KernelRelation.h"
#include "PBCDifference.h"
#include "PositionPredictor.h"
#include <vector>
#include <string>

class AbsorbingBoundaryCondition {
	unsigned _nStep;

	std::vector<data3d> _initialInnerPos;
	std::vector<data3d> _initialOuterPos;
	std::vector<data3d> _currentOuterPos;
	std::vector<std::vector<data3d>> _historyInnerDisplacement;
	std::unique_ptr<KernelRelationBase> _kernelRelationPtr;
	const PBCDifference* _diffFunc;

public:
	void
	Setup (	std::string latticeType, //  need dt, interpolate kernel function
		double timestep,
		double tCut,
		double rCut,
		double cornerCut
	) noexcept;


	void
	SetupInitialState (
		const std::vector<data3d>& equilibriumOuterPos,
		const std::vector<data3d>& equilibriumInnerPos,
		double latticeSpacing,
		const PBCDifference* diffFunc
	) noexcept;


	template < typename RunType >
	const std::vector<data3d>&
	Update (
		const std::vector<data3d>& innerPos,
		RunType runType
	) noexcept
	{
		_currentOuterPos = PositionPredictor() (
			_initialOuterPos,
			*_kernelRelationPtr,
			_historyInnerDisplacement,
			runType
		);
		BackupCurrentInnerPos( innerPos );
		return _currentOuterPos;
	}

private:
	void
	ReserveHistoryDisplacementMemory () noexcept;


	void
	BackupInitialPos (
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos
	) noexcept;


	void
	BackupCurrentInnerPos (
		const std::vector<data3d>& innerPos
	) noexcept;

	//-------------------------------------//

	void
	SetupKernelRelation (
		LatticeType latticeType,
		double tCut,
		double rCut,
		double cornerCut,
		double timestep
	) noexcept;


	void
	KernelRelationPairing (
		double latticeSpacing
	) noexcept;
};

#endif // ABSORBINGBOUNDARYCONDITION_H_INCLUDED

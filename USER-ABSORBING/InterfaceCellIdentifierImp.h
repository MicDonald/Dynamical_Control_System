#ifndef INTERFACECELLIDENTIFIERIMP_H_INCLUDED
#define INTERFACECELLIDENTIFIERIMP_H_INCLUDED

#include "data3d.h"
#include "Cell.h"
#include "PBCDifference.h"
#include <vector>
#include <map>
#include <queue>

class InterfaceCellIdentifierImp {
public:
	constexpr static bool IDENTIFIED = true;
	constexpr static bool UNIDENTIFIED = false;


	struct Neighbors {
		std::vector<unsigned> id;
		std::vector<data3d> diff;
	};


	struct CellCondition {
		std::vector<unsigned> nNeigh2Other;
		std::vector<std::vector<data3d>> neighDiff;
		std::vector<std::vector<data3d>> _outerCellDiff;


		bool
		IsSame (
			const CellCondition& other
		) const noexcept;


		bool
		IsSimilar (
			const CellCondition& other
		) const noexcept;
	};


	struct AtomInfo {
		unsigned cellID;
		unsigned local;
	};

	//-------------------------------------//

	std::pair<std::vector<Cell>, std::vector<Cell>>
	Identify (
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos
	) const noexcept;

protected:
	const PBCDifference* _DiffCorrection;

public:
	InterfaceCellIdentifierImp (
		const PBCDifference* _diffFunc
	) noexcept;


	virtual ~InterfaceCellIdentifierImp () noexcept = default;

private:
	std::vector<Neighbors>
	ParsingNeighboringStatus (
		const std::vector<data3d>& from,
		const std::vector<data3d>& to
	) const noexcept;


	std::vector<Cell>
	IdentifyInnerCells (
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2inner
	) const noexcept;


	std::vector<Cell>
	IdentifyOuterCells (
		const std::vector<Cell>& innerCells,
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos
	) const noexcept;

	//-------------------------------------//

	std::vector<Cell>
	IdentifyInnerCornerCells (
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		std::vector<bool>& atomIdentifiedStatus
	) const noexcept;


	virtual void
	IdentifyInnerEdgeCells (
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		std::vector<bool>& atomIdentifiedStatus,
		std::vector<Cell>& cells
	) const noexcept;


	void
	IdentifyInnerFaceCells (
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		std::vector<bool>& atomIdentifiedStatus,
		std::vector<Cell>& cells
	) const noexcept;

	//-------------------------------------//
protected:
	std::vector<Cell>
	IdentifyInnerCellBasedOn (
		unsigned seed,
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		CellType type,
		CellOrientation orientation = ANY
	) const noexcept;

private:
	std::vector<Cell>
	ProposeCandidates (
		unsigned orientation,
		const std::vector<Neighbors>& inner2Outer
	) const noexcept;


	bool
	ExtendCandidatesCell (
		std::vector<Cell>& candidate,
		const std::vector<Neighbors>& inner2Inner
	) const noexcept;


	void
	PruneCandidates (
		std::vector<Cell>& candidate,
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		CellOrientation orientation,
		CellType cellType = CellType::NTYPES
	) const noexcept;


	virtual bool
	IsValidInnerRelation (
		unsigned base,
		unsigned to,
		const data3d& innerDiff
	) const noexcept = 0;


	virtual bool
	IsValidInnerCell (
		Cell& cell,
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		CellOrientation orientation,
		CellType type
	) const noexcept = 0;


	void
	PrioritySort (
		std::vector<Cell>& innerCells,
		const std::vector<bool>& atomIdentifyStatus
	) const noexcept;


	void
	DifferenceHint (
		unsigned seed,
		const std::vector<Cell>& identifiedCells,
		const std::vector<Neighbors>& inner2Inner,
		std::vector<Cell>& candidateCells
	) const noexcept;

	//-------------------------------------//

	std::vector<Cell>
	OuterCellOf (
		const Cell& cell,
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos
	) const noexcept;


	void
	FindOuterCellAtoms (
		const Cell& cell,
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos,
		InterfaceCoordinate coord,
		std::vector<Cell>& outerCells
	) const noexcept;


	virtual const std::vector<std::vector<data3d>>&
	OuterCellDifference (
		InterfaceCoordinate coord,
		const Cell& cell
	) const noexcept = 0;

	//-------------------------------------//

	virtual void
	ChainingCells (
		std::vector<Cell>& cell,
		const std::vector<data3d>& pos,
		const Cell& refCell,
		const std::vector<data3d>& refPos
	) const noexcept = 0;


	virtual void
	CheckInnerCellChaining (
		std::vector<Cell>& cell,
		const std::vector<data3d>& pos,
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner
	) const noexcept = 0;

	//-------------------------------------//

	void
	AddSeeds (
		const std::vector<Cell>& cells,
		std::queue<std::pair<unsigned, CellOrientation>>& seeds,
		const std::vector<Neighbors>& inner2Inner
	) const noexcept;


	void
	AddSeeds (
		const Cell& cell,
		std::queue<std::pair<unsigned, CellOrientation>>& seeds,
		const std::vector<Neighbors>& inner2Inner
	) const noexcept;

	//-------------------------------------//

	bool
	AllAtomsIdentified (
		std::queue<std::pair<unsigned, CellOrientation>>& seeds,
		const std::vector<bool>& atomIdentifiedStatus
	) const noexcept;


	void
	MarkCellIdentified (
		const Cell& cell,
		std::vector<bool>& atomIdentifiedStatus
	) const noexcept;


	bool
	IsIdentified (
		unsigned seed,
		const std::vector<bool>& atomIdentifiedStatus
	) const noexcept;

	//-------------------------------------//

	virtual unsigned
	NBasis () const noexcept = 0;

protected:
	CellCondition
	ConvertToCellCondition (
		const Cell& cell,
		const std::vector<Neighbors>& inner2Outer
	) const noexcept;


	bool
	OutOfRange (
		const data3d& diff
	) const noexcept;


	data3d
	Difference (
		const data3d& to,
		const data3d& from
	) const noexcept;
};

#endif // INTERFACECELLIDENTIFIERIMP_H_INCLUDED

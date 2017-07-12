#ifndef INTERFACECELLIDENTIFIER2D_H_INCLUDED
#define INTERFACECELLIDENTIFIER2D_H_INCLUDED

#include "data3d.h"
#include "Cell.h"
#include "InterfaceCellIdentifierImp.h"
#include <memory>
#include <vector>
#include <set>
#include <map>

#include <iostream>

class InterfaceCellIdentifier2DImp : public InterfaceCellIdentifierImp {
public:
	enum InterfaceType {
		_X, X_, _Y, Y_, _X_Y, _XY_, X__Y, X_Y_,
		S_X, SX_, S_Y, SY_, S_X_Y, S_XY_, SX__Y, SX_Y_,
		NTYPES
	};

private:
	const std::array<data3d,2>& primitive;
	const std::array<data3d,2>& shift;
	std::map<InterfaceType, CellCondition> _innerCondition;
	std::vector<std::vector<data3d>> _innerCellDiff;
	std::vector<std::vector<data3d>> _SinnerCellDiff;

public:
	InterfaceCellIdentifier2DImp (
		const PBCDifference* diffFunc,
		const std::array<data3d,2>& primitive,
		const std::array<data3d,2>& shift,
		const data3d* basis,
		unsigned N
	) noexcept;

private:
	void
	CalculateNeighboringCondition (
		const std::vector<data3d>& basis
	) noexcept;


	unsigned
	CellNeighboring (
		int shiftX,
		int shiftY,
		const std::vector<data3d>& basis,
		unsigned basisI,
		unsigned basisJ,
		std::vector<data3d>& diff
	) const noexcept;

	//-------------------------------------//

	bool
	IsValidInnerRelation (
		unsigned base,
		unsigned to,
		const data3d& innerDiff
	) const noexcept override;


	bool
	IsValidInnerCell (
		Cell& cell,
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner,
		CellOrientation orientation,
		CellType type
	) const noexcept override;

	//-------------------------------------//

	const std::vector<std::vector<data3d>>&
	OuterCellDifference (
		InterfaceCoordinate coord,
		const Cell& cell
	) const noexcept override;

	//-------------------------------------//

	void
	ChainingCells (
		std::vector<Cell>& cell,
		const std::vector<data3d>& pos,
		const Cell& refCell,
		const std::vector<data3d>& refPos
	) const noexcept override;


	void
	CheckInnerCellChaining (
		std::vector<Cell>& cell,
		const std::vector<data3d>& pos,
		const std::vector<Neighbors>& inner2Outer,
		const std::vector<Neighbors>& inner2Inner
	) const noexcept override;

	//-------------------------------------//

	unsigned
	LatticeLocatingIndex (
		InterfaceType face
	) const noexcept;

	//-------------------------------------//

	void
	SetCellInterface (
		Cell& cell,
		const CellCondition& condition,
		CellOrientation orientation
	 ) const noexcept;


	InterfaceType
	DetailType (
		InterfaceDirection x,
		InterfaceDirection y,
		CellOrientation orientation
	) const noexcept;


	InterfaceType
	DetailType (
		const CellCondition& condition,
		CellOrientation orientation
	) const noexcept;


	CellType
	TypeConvert (
		InterfaceType interfaceType
	) const noexcept;


	CellOrientation
	Orientation (
		InterfaceType interfaceType
	) const noexcept;

	//-------------------------------------//

	unsigned
	NBasis () const noexcept override
	{
		return _innerCellDiff.size();
	}

	//-------------------------------------//

	data3d
	CellBasePosition (
		const Cell& cell,
		const std::vector<data3d>& pos
	) const noexcept;
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

template < typename LATTICE >
class InterfaceCellIdentifier2D {
	InterfaceCellIdentifier2DImp _Imp;

public:
	InterfaceCellIdentifier2D (
		const PBCDifference* _DiffFunc
	) noexcept :
		_Imp ( _DiffFunc, LATTICE::primitive, LATTICE::shift, LATTICE::basis.data(), LATTICE::N )
	{}


	decltype(auto)
	Identify (
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos
	) const noexcept
	{
		return _Imp.Identify( innerPos, outerPos );
	}

};


#endif // INTERFACECELLIDENTIFIER2D_H_INCLUDED

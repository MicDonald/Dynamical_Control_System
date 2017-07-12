#ifndef INTERFACECELLIDENTIFIER3D_H_INCLUDED
#define INTERFACECELLIDENTIFIER3D_H_INCLUDED

#include "InterfaceCellIdentifierImp.h"
#include "Cell.h"
#include "data3d.h"
#include <memory>
#include <vector>
#include <set>
#include <map>

#include <iostream>

class InterfaceCellIdentifier3DImp : public InterfaceCellIdentifierImp {
public:
	enum InterfaceType {
		_X, X_, _Y, Y_, _Z, Z_, S_X, SX_, S_Y, SY_, S_Z, SZ_,
		_X_Y, _XY_, X__Y, X_Y_, S_X_Y, S_XY_, SX__Y, SX_Y_,
		_X_Z, _XZ_, X__Z, X_Z_, S_X_Z, S_XZ_, SX__Z, SX_Z_,
		_Y_Z, _YZ_, Y__Z, Y_Z_, S_Y_Z, S_YZ_, SY__Z, SY_Z_,
		_X_Y_Z, _X_YZ_, _XY__Z, _XY_Z_, X__Y_Z, X__YZ_, X_Y__Z, X_Y_Z_,
		S_X_Y_Z, S_X_YZ_, S_XY__Z, S_XY_Z_, SX__Y_Z, SX__YZ_, SX_Y__Z, SX_Y_Z_,
		NTYPES
	};

private:
	const std::array<data3d,3>& primitive;
	const std::array<data3d,3>& shift;
	std::map<InterfaceType, CellCondition> _innerCondition;
	std::vector<std::vector<data3d>> _innerCellDiff;
	std::vector<std::vector<data3d>> _SinnerCellDiff;

public:
	InterfaceCellIdentifier3DImp (
		const PBCDifference* diffFunc,
		const std::array<data3d, 3>& primitive,
		const std::array<data3d, 3>& shift,
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
		int shiftZ,
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
		InterfaceDirection z,
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
class InterfaceCellIdentifier3D {
	InterfaceCellIdentifier3DImp _Imp;

public:
	InterfaceCellIdentifier3D (
		const PBCDifference* _DiffFunc
	) noexcept :
		_Imp( _DiffFunc, LATTICE::primitive, LATTICE::shift, LATTICE::basis.data(), LATTICE::N )
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


#endif // INTERFACECELLIDENTIFIER3D_H_INCLUDED

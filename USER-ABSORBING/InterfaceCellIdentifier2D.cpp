#include "InterfaceCellIdentifier2D.h"
#include "Lattice.h"
#include <iterator>
#include <list>

#include <iostream>

using namespace std;


InterfaceCellIdentifier2DImp::InterfaceCellIdentifier2DImp (
	const PBCDifference* diffFunc,
	const array<data3d,2>& primitive,
	const array<data3d,2>& shift,
	const data3d* basis,
	unsigned N
) noexcept :
	InterfaceCellIdentifierImp( diffFunc ),
	primitive( primitive ),
	shift( shift )
{
	vector<data3d> _basis( basis, basis+N );
	CalculateNeighboringCondition( _basis );
}

//---------------------------------------------//

void
InterfaceCellIdentifier2DImp::CalculateNeighboringCondition (
	const vector<data3d>& basis
) noexcept
{
	auto shiftBasis = ShiftLattice()( primitive, shift, basis );
	_innerCellDiff.resize( basis.size() );
	_SinnerCellDiff.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		_innerCellDiff[i].resize( basis.size() );
		_SinnerCellDiff[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			if ( i==j )
			{
				_innerCellDiff[i][j].fill( 0. );
				_SinnerCellDiff[i][j].fill( 0. );
				continue;
			}

			_innerCellDiff[i][j] = Substract( basis[j], basis[i] );
			_SinnerCellDiff[i][j] = Substract( shiftBasis[j], shiftBasis[i] );
		}
	}

	// LEFT
	auto& inner2OuterL = _innerCondition[_X].nNeigh2Other;
	auto& inner2OuterSL = _innerCondition[S_X].nNeigh2Other;
	auto& neighDiffL = _innerCondition[_X].neighDiff;
	auto& neighDiffSL = _innerCondition[S_X].neighDiff;
	inner2OuterL.assign( basis.size(), 0 );
	inner2OuterSL.assign( basis.size(), 0 );
	neighDiffL.resize( basis.size() );
	neighDiffSL.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for(	unsigned j=0; j<basis.size(); ++j)
		{
			inner2OuterL[i] += CellNeighboring(-1, -1, basis, i, j, neighDiffL[i]);
			inner2OuterL[i] += CellNeighboring(-1, 0, basis, i, j, neighDiffL[i]);
			inner2OuterL[i] += CellNeighboring(-1, 1, basis, i, j, neighDiffL[i]);
			inner2OuterSL[i] += CellNeighboring(-1, -1, shiftBasis, i, j, neighDiffSL[i]);
			inner2OuterSL[i] += CellNeighboring(-1, 0, shiftBasis, i, j, neighDiffSL[i]);
			inner2OuterSL[i] += CellNeighboring(-1, 1, shiftBasis, i, j, neighDiffSL[i]);
		}
	auto& outerCellDiffL = _innerCondition[_X]._outerCellDiff;
	auto& outerCellDiffSL = _innerCondition[S_X]._outerCellDiff;
	outerCellDiffL.resize( basis.size() );
	outerCellDiffSL.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffL[i].resize( basis.size() );
		outerCellDiffSL[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffL[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0], _innerCellDiff[i][j][1]-primitive[0][1] };
			outerCellDiffSL[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0], _SinnerCellDiff[i][j][1]-primitive[0][1] };
		}
	}

	// RIGHT
	auto& inner2OuterR = _innerCondition[X_].nNeigh2Other;
	auto& inner2OuterSR = _innerCondition[SX_].nNeigh2Other;
	auto& neighDiffR = _innerCondition[X_].neighDiff;
	auto& neighDiffSR = _innerCondition[SX_].neighDiff;
	inner2OuterR.assign( basis.size(), 0 );
	inner2OuterSR.assign( basis.size(), 0 );
	neighDiffR.resize( basis.size() );
	neighDiffSR.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterR[i] += CellNeighboring(1, -1, basis, i, j, neighDiffR[i]);
			inner2OuterR[i] += CellNeighboring(1, 0, basis, i, j, neighDiffR[i]);
			inner2OuterR[i] += CellNeighboring(1, 1, basis, i, j, neighDiffR[i]);
			inner2OuterSR[i] += CellNeighboring(1, -1, shiftBasis, i, j, neighDiffSR[i]);
			inner2OuterSR[i] += CellNeighboring(1, 0, shiftBasis, i, j, neighDiffSR[i]);
			inner2OuterSR[i] += CellNeighboring(1, 1, shiftBasis, i, j, neighDiffSR[i]);
		}
	auto& outerCellDiffR = _innerCondition[X_]._outerCellDiff;
	auto& outerCellDiffSR = _innerCondition[SX_]._outerCellDiff;
	outerCellDiffR.resize( basis.size() );
	outerCellDiffSR.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffR[i].resize( basis.size() );
		outerCellDiffSR[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffR[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0], _innerCellDiff[i][j][1]+primitive[0][1] };
			outerCellDiffSR[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0], _SinnerCellDiff[i][j][1]+primitive[0][1] };
		}
	}

	// UP
	auto& inner2OuterU = _innerCondition[Y_].nNeigh2Other;
	auto& inner2OuterSU = _innerCondition[SY_].nNeigh2Other;
	auto& neighDiffU = _innerCondition[Y_].neighDiff;
	auto& neighDiffSU = _innerCondition[SY_].neighDiff;
	inner2OuterU.assign( basis.size(), 0 );
	inner2OuterSU.assign( basis.size(), 0 );
	neighDiffU.resize( basis.size() );
	neighDiffSU.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterU[i] += CellNeighboring(1, 1, basis, i, j, neighDiffU[i]);
			inner2OuterU[i] += CellNeighboring(0, 1, basis, i, j, neighDiffU[i]);
			inner2OuterU[i] += CellNeighboring(-1, 1, basis, i, j, neighDiffU[i]);
			inner2OuterSU[i] += CellNeighboring(1, 1, shiftBasis, i, j, neighDiffSU[i]);
			inner2OuterSU[i] += CellNeighboring(0, 1, shiftBasis, i, j, neighDiffSU[i]);
			inner2OuterSU[i] += CellNeighboring(-1, 1, shiftBasis, i, j, neighDiffSU[i]);
		}
	auto& outerCellDiffU = _innerCondition[Y_]._outerCellDiff;
	auto& outerCellDiffSU = _innerCondition[SY_]._outerCellDiff;
	outerCellDiffU.resize( basis.size() );
	outerCellDiffSU.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffU[i].resize( basis.size() );
		outerCellDiffSU[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffU[i][j] = { _innerCellDiff[i][j][0]+primitive[1][0], _innerCellDiff[i][j][1]+primitive[1][1] };
			outerCellDiffSU[i][j] = { _SinnerCellDiff[i][j][0]+primitive[1][0], _SinnerCellDiff[i][j][1]+primitive[1][1] };
		}
	}

	// DOWN
	auto& inner2OuterD = _innerCondition[_Y].nNeigh2Other;
	auto& inner2OuterSD = _innerCondition[S_Y].nNeigh2Other;
	auto& neighDiffD = _innerCondition[_Y].neighDiff;
	auto& neighDiffSD = _innerCondition[S_Y].neighDiff;
	inner2OuterD.assign( basis.size(), 0 );
	inner2OuterSD.assign( basis.size(), 0 );
	neighDiffD.resize( basis.size() );
	neighDiffSD.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterD[i] += CellNeighboring(1, -1, basis, i, j, neighDiffD[i]);
			inner2OuterD[i] += CellNeighboring(0, -1, basis, i, j, neighDiffD[i]);
			inner2OuterD[i] += CellNeighboring(-1, -1, basis, i, j, neighDiffD[i]);
			inner2OuterSD[i] += CellNeighboring(1, -1, shiftBasis, i, j, neighDiffSD[i]);
			inner2OuterSD[i] += CellNeighboring(0, -1, shiftBasis, i, j, neighDiffSD[i]);
			inner2OuterSD[i] += CellNeighboring(-1, -1, shiftBasis, i, j, neighDiffSD[i]);
		}
	auto& outerCellDiffD = _innerCondition[_Y]._outerCellDiff;
	auto& outerCellDiffSD = _innerCondition[S_Y]._outerCellDiff;
	outerCellDiffD.resize( basis.size() );
	outerCellDiffSD.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffD[i].resize( basis.size() );
		outerCellDiffSD[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffD[i][j] = { _innerCellDiff[i][j][0]-primitive[1][0], _innerCellDiff[i][j][1]-primitive[1][1] };
			outerCellDiffSD[i][j] = { _SinnerCellDiff[i][j][0]-primitive[1][0], _SinnerCellDiff[i][j][1]-primitive[1][1] };
		}
	}

	//LU
	auto& inner2OuterLU = _innerCondition[_XY_].nNeigh2Other;
	auto& inner2OuterSLU = _innerCondition[S_XY_].nNeigh2Other;
	auto& neighDiffLU = _innerCondition[_XY_].neighDiff;
	auto& neighDiffSLU = _innerCondition[S_XY_].neighDiff;
	inner2OuterLU.insert( inner2OuterLU.end(), basis.size(), 0 );
	inner2OuterSLU.insert( inner2OuterSLU.end(), basis.size(), 0 );
	neighDiffLU.resize( basis.size() );
	neighDiffSLU.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterLU[i] += CellNeighboring(-1, -1, basis, i, j, neighDiffLU[i]);
			inner2OuterLU[i] += CellNeighboring(-1, 0, basis, i, j, neighDiffLU[i]);
			inner2OuterLU[i] += CellNeighboring(-1, 1, basis, i, j, neighDiffLU[i]);
			inner2OuterLU[i] += CellNeighboring(0, 1, basis, i, j, neighDiffLU[i]);
			inner2OuterLU[i] += CellNeighboring(1, 1, basis, i, j, neighDiffLU[i]);
			inner2OuterSLU[i] += CellNeighboring(-1, -1, shiftBasis, i, j, neighDiffSLU[i]);
			inner2OuterSLU[i] += CellNeighboring(-1, 0, shiftBasis, i, j, neighDiffSLU[i]);
			inner2OuterSLU[i] += CellNeighboring(-1, 1, shiftBasis, i, j, neighDiffSLU[i]);
			inner2OuterSLU[i] += CellNeighboring(0, 1, shiftBasis, i, j, neighDiffSLU[i]);
			inner2OuterSLU[i] += CellNeighboring(1, 1, shiftBasis, i, j, neighDiffSLU[i]);
		}
	auto& outerCellDiffLU = _innerCondition[_XY_]._outerCellDiff;
	auto& outerCellDiffSLU = _innerCondition[S_XY_]._outerCellDiff;
	outerCellDiffLU.resize( basis.size() );
	outerCellDiffSLU.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffLU[i].resize( basis.size() );
		outerCellDiffSLU[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffLU[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0], _innerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1] };
			outerCellDiffSLU[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0], _SinnerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1] };
		}
	}

	//RU
	auto& inner2OuterRU = _innerCondition[X_Y_].nNeigh2Other;
	auto& inner2OuterSRU = _innerCondition[SX_Y_].nNeigh2Other;
	auto& neighDiffRU = _innerCondition[X_Y_].neighDiff;
	auto& neighDiffSRU = _innerCondition[SX_Y_].neighDiff;
	inner2OuterRU.assign( basis.size(), 0 );
	inner2OuterSRU.assign( basis.size(), 0 );
	neighDiffRU.resize( basis.size() );
	neighDiffSRU.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterRU[i] += CellNeighboring(-1, 1, basis, i, j, neighDiffRU[i]);
			inner2OuterRU[i] += CellNeighboring(0, 1, basis, i, j, neighDiffRU[i]);
			inner2OuterRU[i] += CellNeighboring(1, 1, basis, i, j, neighDiffRU[i]);
			inner2OuterRU[i] += CellNeighboring(1, 0, basis, i, j, neighDiffRU[i]);
			inner2OuterRU[i] += CellNeighboring(1, -1, basis, i, j, neighDiffRU[i]);
			inner2OuterSRU[i] += CellNeighboring(-1, 1, shiftBasis, i, j, neighDiffSRU[i]);
			inner2OuterSRU[i] += CellNeighboring(0, 1, shiftBasis, i, j, neighDiffSRU[i]);
			inner2OuterSRU[i] += CellNeighboring(1, 1, shiftBasis, i, j, neighDiffSRU[i]);
			inner2OuterSRU[i] += CellNeighboring(1, 0, shiftBasis, i, j, neighDiffSRU[i]);
			inner2OuterSRU[i] += CellNeighboring(1, -1, shiftBasis, i, j, neighDiffSRU[i]);
		}
	auto& outerCellDiffRU = _innerCondition[X_Y_]._outerCellDiff;
	auto& outerCellDiffSRU = _innerCondition[SX_Y_]._outerCellDiff;
	outerCellDiffRU.resize( basis.size() );
	outerCellDiffSRU.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffRU[i].resize( basis.size() );
		outerCellDiffSRU[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffRU[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0], _innerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1] };
			outerCellDiffSRU[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0], _SinnerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1] };
		}
	}

	//LD
	auto& inner2OuterLD = _innerCondition[_X_Y].nNeigh2Other;
	auto& inner2OuterSLD = _innerCondition[S_X_Y].nNeigh2Other;
	auto& neighDiffLD = _innerCondition[_X_Y].neighDiff;
	auto& neighDiffSLD = _innerCondition[S_X_Y].neighDiff;
	inner2OuterLD.assign( basis.size(), 0 );
	inner2OuterSLD.assign( basis.size(), 0 );
	neighDiffLD.resize( basis.size() );
	neighDiffSLD.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterLD[i] += CellNeighboring(-1, 1, basis, i, j, neighDiffLD[i]);
			inner2OuterLD[i] += CellNeighboring(-1, 0, basis, i, j, neighDiffLD[i]);
			inner2OuterLD[i] += CellNeighboring(-1, -1, basis, i, j, neighDiffLD[i]);
			inner2OuterLD[i] += CellNeighboring(0, -1, basis, i, j, neighDiffLD[i]);
			inner2OuterLD[i] += CellNeighboring(1, -1, basis, i, j, neighDiffLD[i]);
			inner2OuterSLD[i] += CellNeighboring(-1, 1, shiftBasis, i, j, neighDiffSLD[i]);
			inner2OuterSLD[i] += CellNeighboring(-1, 0, shiftBasis, i, j, neighDiffSLD[i]);
			inner2OuterSLD[i] += CellNeighboring(-1, -1, shiftBasis, i, j, neighDiffSLD[i]);
			inner2OuterSLD[i] += CellNeighboring(0, -1, shiftBasis, i, j, neighDiffSLD[i]);
			inner2OuterSLD[i] += CellNeighboring(1, -1, shiftBasis, i, j, neighDiffSLD[i]);
		}
	auto& outerCellDiffLD = _innerCondition[_X_Y]._outerCellDiff;
	auto& outerCellDiffSLD = _innerCondition[S_X_Y]._outerCellDiff;
	outerCellDiffLD.resize( basis.size() );
	outerCellDiffSLD.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffLD[i].resize( basis.size() );
		outerCellDiffSLD[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffLD[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0], _innerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1] };
			outerCellDiffSLD[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0], _SinnerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1] };
		}
	}

	//RD
	auto& inner2OuterRD = _innerCondition[X__Y].nNeigh2Other;
	auto& inner2OuterSRD = _innerCondition[SX__Y].nNeigh2Other;
	auto& neighDiffRD = _innerCondition[X__Y].neighDiff;
	auto& neighDiffSRD = _innerCondition[SX__Y].neighDiff;
	inner2OuterRD.assign( basis.size(), 0 );
	inner2OuterSRD.assign( basis.size(), 0 );
	neighDiffRD.resize( basis.size() );
	neighDiffSRD.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterRD[i] += CellNeighboring(-1, -1, basis, i, j, neighDiffRD[i]);
			inner2OuterRD[i] += CellNeighboring(0, -1, basis, i, j, neighDiffRD[i]);
			inner2OuterRD[i] += CellNeighboring(1, -1, basis, i, j, neighDiffRD[i]);
			inner2OuterRD[i] += CellNeighboring(1, 0, basis, i, j, neighDiffRD[i]);
			inner2OuterRD[i] += CellNeighboring(1, 1, basis, i, j, neighDiffRD[i]);
			inner2OuterSRD[i] += CellNeighboring(-1, -1, shiftBasis, i, j, neighDiffSRD[i]);
			inner2OuterSRD[i] += CellNeighboring(0, -1, shiftBasis, i, j, neighDiffSRD[i]);
			inner2OuterSRD[i] += CellNeighboring(1, -1, shiftBasis, i, j, neighDiffSRD[i]);
			inner2OuterSRD[i] += CellNeighboring(1, 0, shiftBasis, i, j, neighDiffSRD[i]);
			inner2OuterSRD[i] += CellNeighboring(1, 1, shiftBasis, i, j, neighDiffSRD[i]);
		}
	auto& outerCellDiffRD = _innerCondition[X__Y]._outerCellDiff;
	auto& outerCellDiffSRD = _innerCondition[SX__Y]._outerCellDiff;
	outerCellDiffRD.resize( basis.size() );
	outerCellDiffSRD.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffRD[i].resize( basis.size() );
		outerCellDiffSRD[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffRD[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0], _innerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1] };
			outerCellDiffSRD[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0], _SinnerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1] };
		}
	}
}


unsigned
InterfaceCellIdentifier2DImp::CellNeighboring (
	int shiftX,
	int shiftY,
	const vector<data3d>& basis,
	unsigned basisI,
	unsigned basisJ,
	vector<data3d>& neighDiff
) const noexcept
{
	data3d shift{};
	shift[0] = primitive[0][0]*shiftX + primitive[1][0]*shiftY +
		basis[basisJ][0] - basis[basisI][0];
	shift[1] = primitive[0][1]*shiftX + primitive[1][1]*shiftY +
		basis[basisJ][1] - basis[basisI][1];
	if ( OutOfRange(shift) )
		return 0;
	neighDiff.emplace_back( move(shift) );
	return 1;
}

//---------------------------------------------------------------------------//

bool
InterfaceCellIdentifier2DImp::IsValidInnerRelation (
	unsigned base,
	unsigned to,
	const data3d& innerDiff
) const noexcept
{
	if ( Equivalent( _innerCellDiff[base][to], innerDiff ) )
		return true;
	if ( Equivalent( _SinnerCellDiff[base][to], innerDiff ) )
		return true;
	return false;
}


bool
InterfaceCellIdentifier2DImp::IsValidInnerCell (
	Cell& cell,
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner,
	CellOrientation orientation,
	CellType type
) const noexcept
{
	auto nBasis = NBasis();
	auto nInners = inner2Outer.size();
	bool regularCase = true;
	bool shiftCase = true;
	for ( unsigned i=0; i<nBasis; ++i )
	{
		if ( cell.id[i] >= nInners ) continue;
		for ( unsigned j=i+1; j<nBasis; ++j )
		{
			if ( cell.id[j] >= nInners ) continue;
			auto diffIter = inner2Inner[ cell.id[i] ].diff.begin();
			auto idIter = inner2Inner[ cell.id[i] ].id.begin();
			for (	; diffIter != inner2Inner[ cell.id[i] ].diff.end();
				++diffIter, ++idIter
			)
			{
				if ( *idIter != cell.id[j] ) continue;
				auto diff = Substract( *diffIter, _innerCellDiff[i][j] );
				auto diffS = Substract( *diffIter, _SinnerCellDiff[i][j] );
				if ( Square(diff) > 1.e-5 ) regularCase = false;
				if ( Square(diffS) > 1.e-5 ) shiftCase = false;
				break;
			}
		}
	}
	if ( !regularCase && !shiftCase ) return false;
	else if ( orientation == REGULAR && !regularCase )
		return false;
	else if ( orientation == SHIFT && !shiftCase )
		return false;

	auto condition = ConvertToCellCondition( cell, inner2Outer );
	if ( type != CellType::NTYPES )
	{
		for ( const auto& innerCondition : _innerCondition )
		{
			if ( TypeConvert(innerCondition.first) != type ) continue;
			if ( innerCondition.second.IsSame(condition) )
			{
				SetCellInterface( cell, condition, orientation );
				return true;
			}
		}
	}
	else
	{
		for ( const auto& innerCondition : _innerCondition )
		{
			if ( orientation != ANY )
				if ( Orientation(innerCondition.first) != orientation )
					continue;
			if ( innerCondition.second.IsSimilar(condition) )
				return true;
		}
	}
	return false;
}

//---------------------------------------------------------------------------//

const std::vector<std::vector<data3d>>&
InterfaceCellIdentifier2DImp::OuterCellDifference (
	InterfaceCoordinate coord,
	const Cell& cell
) const noexcept
{
	switch ( coord )
	{
	case X:
		return cell.x == MINUS ?
			( cell.orientation == REGULAR ?
				_innerCondition.at(_X)._outerCellDiff :
				_innerCondition.at(S_X)._outerCellDiff ) :
			( cell.orientation == REGULAR ?
				_innerCondition.at(X_)._outerCellDiff :
				_innerCondition.at(SX_)._outerCellDiff );
	case Y:
		return cell.y == MINUS ?
			( cell.orientation == REGULAR ?
				_innerCondition.at(_Y)._outerCellDiff :
				_innerCondition.at(S_Y)._outerCellDiff ) :
			( cell.orientation == REGULAR ?
				_innerCondition.at(Y_)._outerCellDiff :
				_innerCondition.at(SY_)._outerCellDiff );
	case XY:
		if ( cell.x == MINUS && cell.y == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_X_Y)._outerCellDiff :
				_innerCondition.at(S_X_Y)._outerCellDiff;
		else if ( cell.x == MINUS && cell.y == PLUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_XY_)._outerCellDiff :
				_innerCondition.at(S_XY_)._outerCellDiff;
		else if ( cell.x == PLUS && cell.y == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(X__Y)._outerCellDiff :
				_innerCondition.at(SX__Y)._outerCellDiff;
		else // cell.x == PLUS && cell.y == PLUS
			return cell.orientation == REGULAR ?
				_innerCondition.at(X_Y_)._outerCellDiff :
				_innerCondition.at(SX_Y_)._outerCellDiff;
	default :
		cerr << "check 1"<<endl;
		return _innerCellDiff;
	}
}

//---------------------------------------------------------------------------//

void
InterfaceCellIdentifier2DImp::ChainingCells (
	vector<Cell>& cells,
	const vector<data3d>& pos,
	const Cell& refCell,
	const vector<data3d>& refPos
) const noexcept
{
	if ( cells.empty() ) return;

	auto basePos = CellBasePosition( refCell, refPos );
	for ( auto& cell : cells )
	{
		auto cellPos = CellBasePosition( cell, pos );
		cell.position[0] = (cellPos[0]-basePos[0]) / primitive[0][0];
		cell.position[1] = (cellPos[1]-basePos[1]) / primitive[1][1];
		cell.position[2] = 0;
	}
}


void
InterfaceCellIdentifier2DImp::CheckInnerCellChaining (
	vector<Cell>& cells,
	const vector<data3d>& pos,
	const vector<Neighbors>& inner2Outer,
	const vector<Neighbors>& inner2Inner
) const noexcept
{
	for (	auto iter = cells.begin();
		iter != cells.end();
		++iter
	)
	{
		// find known face
		bool _x = iter->x == MINUS ? true : false;
		bool x_ = iter->x == PLUS ? true : false;
		bool _y = iter->y == MINUS ? true : false;
		bool y_ = iter->y == PLUS ? true : false;
		if ( iter->x == NONE )
			_y = y_ = true;
		if ( iter->y == NONE )
			_x = x_ = true;
		for (	auto iterJ = cells.begin();
			iterJ != cells.end();
			++iterJ
		)
		{
			if ( iter == iterJ ) continue;
			auto diff = Difference (
				data3d{iterJ->position[0]*primitive[0][0], iterJ->position[1]*primitive[1][1]},
				data3d{iter->position[0]*primitive[0][0], iter->position[1]*primitive[1][1]}
			);
			diff[0] /= primitive[0][0];
			diff[1] /= primitive[1][1];
			if ( OutOfRange(diff) ) continue;
			if ( diff[0] < 0. ) _x = true;
			if ( diff[0] > 0. ) x_ = true;
			if ( diff[1] < 0. ) _y = true;
			if ( diff[1] > 0. ) y_ = true;
		}
		if ( _x && x_ && _y && y_ ) continue; // connecting cells are known

		// find cell in the unknown cell face
		InterfaceDirection x = NONE;
		InterfaceDirection y = NONE;
		if ( !_x ) x = MINUS;
		else if ( !x_ ) x = PLUS;
		if ( !_y ) y = MINUS;
		else if ( !y_ ) y = PLUS;
		auto interfaceType = DetailType( x, y, iter->orientation );
		CellType cellType = TypeConvert( interfaceType );
		unsigned index = LatticeLocatingIndex( interfaceType );
		unsigned baseID = iter->id[ index ];
		CellOrientation orientation = iter->orientation == REGULAR ? SHIFT : REGULAR;
		if ( baseID >= inner2Inner.size() ) continue;
		auto extraCells = IdentifyInnerCellBasedOn (
				baseID,
				inner2Outer, inner2Inner,
				cellType,
				orientation
			);
		ChainingCells( extraCells, pos, cells[0], pos );
		cells.emplace_back( move(extraCells[0]) );
	}
}

//---------------------------------------------------------------------------//

unsigned
InterfaceCellIdentifier2DImp::LatticeLocatingIndex (
	InterfaceType type
) const noexcept
{
	auto orientation = Orientation( type );
	for (	unsigned index = 0;
		index < _innerCellDiff.size();
		++index
	)
	{
		unsigned _x = 0;
		unsigned x_ = 0;
		unsigned _y = 0;
		unsigned y_ = 0;
		switch ( orientation )
		{
		case REGULAR :
			for ( const auto& diff : _innerCellDiff[index] )
			{
				if ( diff[0] > 0. ) ++_x;
				else if ( diff[0] < 0. ) ++x_;
				if ( diff[1] > 0. ) ++_y;
				else if ( diff[1] < 0. ) ++y_;
			}
			break;
		case SHIFT :
			for ( const auto& diff : _SinnerCellDiff[index] )
			{
				if ( diff[0] > 0. ) ++_x;
				else if ( diff[0] < 0. ) ++x_;
				if ( diff[1] > 0. ) ++_y;
				else if ( diff[1] < 0. ) ++y_;
			}
			break;
		default : // ANY
			;
		}

		switch ( type )
		{
		case _X :
		case S_X :
			if ( _x > x_ ) return index;
			break;
		case X_ :
		case SX_ :
			if ( x_ > _x ) return index;
			break;
		case _Y :
		case S_Y :
			if ( _y > y_ ) return index;
			break;
		case Y_ :
		case SY_ :
			if ( y_ > _y ) return index;
			break;
		case _X_Y :
		case S_X_Y :
			if ( _x > x_ && _y > y_ ) return index;
			break;
		case _XY_ :
		case S_XY_ :
			if ( _x > x_ && y_ > _y ) return index;
			break;
		case X__Y :
		case SX__Y :
			if ( x_ > _x && _y > y_ ) return index;
			break;
		case X_Y_ :
		case SX_Y_ :
			if ( x_ > _x && y_ > _y ) return index;
			break;
		default :
			cout << "check,,,"<<endl;
			return _innerCellDiff.size();
		}
	}
	return _innerCellDiff.size();
}

//---------------------------------------------------------------------------//

void
InterfaceCellIdentifier2DImp::SetCellInterface (
	Cell& cell,
	const CellCondition& condition,
	CellOrientation orientation
) const noexcept
{
	auto detailType = DetailType( condition, orientation );
	cell.orientation = Orientation( detailType );
	switch ( detailType )
	{
	case S_X :
	case _X :
		cell.x = MINUS;
		break;
	case SX_ :
	case X_ :
		cell.x = PLUS;
		break;
	case S_Y :
	case _Y :
		cell.y = MINUS;
		break;
	case SY_ :
	case Y_ :
		cell.y = PLUS;
		break;
	case S_X_Y :
	case _X_Y :
		cell.x = MINUS;
		cell.y = MINUS;
		break;
	case S_XY_ :
	case _XY_ :
		cell.x = MINUS;
		cell.y = PLUS;
		break;
	case SX__Y :
	case X__Y :
		cell.x = PLUS;
		cell.y = MINUS;
		break;
	case SX_Y_ :
	case X_Y_ :
		cell.x = PLUS;
		cell.y = PLUS;
		break;
	default :
		cerr << "check 2"<<endl;
		break;
	}
}


InterfaceCellIdentifier2DImp::InterfaceType
InterfaceCellIdentifier2DImp::DetailType (
	InterfaceDirection x,
	InterfaceDirection y,
	CellOrientation orientation
) const noexcept
{
	if ( orientation == ANY )
	{
		cout << "check---" << endl;
		return NTYPES;
	}

	switch ( x )
	{
	case MINUS :
		switch ( y )
		{
		case NONE :
			return orientation==REGULAR ? _X : S_X;
		case MINUS :
			return orientation==REGULAR ? _X_Y : S_X_Y;
		case PLUS :
			return orientation==REGULAR ? _XY_ : S_XY_;
		}
	case PLUS :
		switch ( y )
		{
		case NONE :
			return orientation==REGULAR ? X_ : SX_;
		case MINUS :
			return orientation==REGULAR ? X__Y : SX__Y;
		case PLUS :
			return orientation==REGULAR ? X_Y_ : SX_Y_;
		}
	case NONE :
		switch ( y )
		{
		case NONE :
			cout << "check..." << endl;
			return NTYPES;
		case MINUS :
			return orientation==REGULAR ? _Y : S_Y;
		case PLUS :
			return orientation==REGULAR ? Y_ : SY_;
		}
	default :
		cout << "check..." << endl;
		return NTYPES;
	}
}


InterfaceCellIdentifier2DImp::InterfaceType
InterfaceCellIdentifier2DImp::DetailType (
	const CellCondition& condition,
	CellOrientation orientation
) const noexcept
{
	for (	auto iter = _innerCondition.begin();
		iter != _innerCondition.end();
		++iter
	)
	{
		if ( orientation != ANY )
			if ( Orientation(iter->first) != orientation )
				continue;
		if ( iter->second.IsSame(condition) )
			return iter->first;
	}
	return InterfaceType::NTYPES;
}


CellType
InterfaceCellIdentifier2DImp::TypeConvert (
	InterfaceType interfaceType
) const noexcept
{
	switch ( interfaceType )
	{
	case _X :
	case S_X :
	case X_ :
	case SX_ :
	case _Y :
	case S_Y :
	case Y_ :
	case SY_ :
		return FACE;
	case _X_Y :
	case S_X_Y :
	case _XY_ :
	case S_XY_ :
	case X__Y :
	case SX__Y :
	case X_Y_ :
	case SX_Y_ :
		return EDGE;
	default :
		cout<<"check 3"<<endl;
		return CellType::NTYPES;
	}
}


CellOrientation
InterfaceCellIdentifier2DImp::Orientation (
	InterfaceType interfaceType
) const noexcept
{
	switch ( interfaceType )
	{
	case _X :
	case X_ :
	case _Y :
	case Y_ :
	case _X_Y :
	case _XY_ :
	case X__Y :
	case X_Y_ :
		return REGULAR;
	case S_X :
	case SX_ :
	case S_Y :
	case SY_ :
	case S_X_Y :
	case S_XY_ :
	case SX__Y :
	case SX_Y_ :
		return SHIFT;
	default :
		cout<<"check 4"<<endl;
		return ANY;
	}
}


data3d
InterfaceCellIdentifier2DImp::CellBasePosition (
	const Cell& cell,
	const vector<data3d>& pos
) const noexcept
{
	auto nAtoms = pos.size();
	unsigned index = 0;
	while ( cell.id[index] >= nAtoms )
	{
		++index;
	}

	data3d sum{};
	switch ( cell.orientation )
	{
	case REGULAR :
		for ( const auto& diff : _innerCellDiff[index] )
			Increment( sum, diff );
		break;
	case SHIFT :
		for ( const auto& diff : _SinnerCellDiff[index] )
			Increment( sum, diff );
		break;
	default:
		return {};
	}

	auto nBasis = _innerCellDiff[index].size();
	auto inverse = 1./nBasis;
	sum[0] *= inverse;
	sum[1] *= inverse;
	Increment( sum, pos[ cell.id[index] ] );
	return sum;
}



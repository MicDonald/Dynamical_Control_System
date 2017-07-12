#include "InterfaceCellIdentifier3D.h"
#include "Lattice.h"
#include <iterator>
#include <list>

#include <iostream>

using namespace std;


InterfaceCellIdentifier3DImp::InterfaceCellIdentifier3DImp (
	const PBCDifference* diffFunc,
	const array<data3d, 3>& primitive,
	const array<data3d, 3>& shift,
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
InterfaceCellIdentifier3DImp::CalculateNeighboringCondition (
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

	// - 0 0
	auto& inner2Outer_X = _innerCondition[_X].nNeigh2Other;
	auto& inner2OuterS_X = _innerCondition[S_X].nNeigh2Other;
	auto& neighDiff_X = _innerCondition[_X].neighDiff;
	auto& neighDiffS_X = _innerCondition[S_X].neighDiff;
	inner2Outer_X.assign( basis.size(), 0 );
	inner2OuterS_X.assign( basis.size(), 0 );
	neighDiff_X.resize( basis.size() );
	neighDiffS_X.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_X[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_X[i] );
			inner2Outer_X[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_X[i] );

			inner2OuterS_X[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_X[i] );
			inner2OuterS_X[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_X[i] );
		}
	auto& outerCellDiff_X = _innerCondition[_X]._outerCellDiff;
	auto& outerCellDiffS_X = _innerCondition[S_X]._outerCellDiff;
	outerCellDiff_X.resize( basis.size() );
	outerCellDiffS_X.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_X[i].resize( basis.size() );
		outerCellDiffS_X[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_X[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0], _innerCellDiff[i][j][1]-primitive[0][1], _innerCellDiff[i][j][2]-primitive[0][2]  };
			outerCellDiffS_X[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0], _SinnerCellDiff[i][j][1]-primitive[0][1], _SinnerCellDiff[i][j][2]-primitive[0][2] };
		}
	}

	// + 0 0
	auto& inner2OuterX_ = _innerCondition[X_].nNeigh2Other;
	auto& inner2OuterSX_ = _innerCondition[SX_].nNeigh2Other;
	auto& neighDiffX_ = _innerCondition[X_].neighDiff;
	auto& neighDiffSX_ = _innerCondition[SX_].neighDiff;
	inner2OuterX_.assign( basis.size(), 0 );
	inner2OuterSX_.assign( basis.size(), 0 );
	neighDiffX_.resize( basis.size() );
	neighDiffSX_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX_[i] );
			inner2OuterX_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX_[i] );

			inner2OuterSX_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX_[i] );
			inner2OuterSX_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX_[i] );
		}
	auto& outerCellDiffX_ = _innerCondition[X_]._outerCellDiff;
	auto& outerCellDiffSX_ = _innerCondition[SX_]._outerCellDiff;
	outerCellDiffX_.resize( basis.size() );
	outerCellDiffSX_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX_[i].resize( basis.size() );
		outerCellDiffSX_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX_[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0], _innerCellDiff[i][j][1]+primitive[0][1], _innerCellDiff[i][j][2]+primitive[0][2]  };
			outerCellDiffSX_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0], _SinnerCellDiff[i][j][1]+primitive[0][1], _SinnerCellDiff[i][j][2]+primitive[0][2] };
		}
	}

	// 0 - 0
	auto& inner2Outer_Y = _innerCondition[_Y].nNeigh2Other;
	auto& inner2OuterS_Y = _innerCondition[S_Y].nNeigh2Other;
	auto& neighDiff_Y = _innerCondition[_Y].neighDiff;
	auto& neighDiffS_Y = _innerCondition[S_Y].neighDiff;
	inner2Outer_Y.assign( basis.size(), 0 );
	inner2OuterS_Y.assign( basis.size(), 0 );
	neighDiff_Y.resize( basis.size() );
	neighDiffS_Y.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_Y[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiff_Y[i] );
			inner2Outer_Y[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_Y[i] );

			inner2OuterS_Y[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffS_Y[i] );
			inner2OuterS_Y[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_Y[i] );
		}
	auto& outerCellDiff_Y = _innerCondition[_Y]._outerCellDiff;
	auto& outerCellDiffS_Y = _innerCondition[S_Y]._outerCellDiff;
	outerCellDiff_Y.resize( basis.size() );
	outerCellDiffS_Y.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_Y[i].resize( basis.size() );
		outerCellDiffS_Y[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_Y[i][j] = { _innerCellDiff[i][j][0]-primitive[1][0], _innerCellDiff[i][j][1]-primitive[1][1], _innerCellDiff[i][j][2]-primitive[1][2]  };
			outerCellDiffS_Y[i][j] = { _SinnerCellDiff[i][j][0]-primitive[1][0], _SinnerCellDiff[i][j][1]-primitive[1][1], _SinnerCellDiff[i][j][2]-primitive[1][2] };
		}
	}

	// 0 + 0
	auto& inner2OuterY_ = _innerCondition[Y_].nNeigh2Other;
	auto& inner2OuterSY_ = _innerCondition[SY_].nNeigh2Other;
	auto& neighDiffY_ = _innerCondition[Y_].neighDiff;
	auto& neighDiffSY_ = _innerCondition[SY_].neighDiff;
	inner2OuterY_.assign( basis.size(), 0 );
	inner2OuterSY_.assign( basis.size(), 0 );
	neighDiffY_.resize( basis.size() );
	neighDiffSY_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterY_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffY_[i] );
			inner2OuterY_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffY_[i] );

			inner2OuterSY_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSY_[i] );
			inner2OuterSY_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSY_[i] );
		}
	auto& outerCellDiffY_ = _innerCondition[Y_]._outerCellDiff;
	auto& outerCellDiffSY_ = _innerCondition[SY_]._outerCellDiff;
	outerCellDiffY_.resize( basis.size() );
	outerCellDiffSY_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffY_[i].resize( basis.size() );
		outerCellDiffSY_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffY_[i][j] = { _innerCellDiff[i][j][0]+primitive[1][0], _innerCellDiff[i][j][1]+primitive[1][1], _innerCellDiff[i][j][2]+primitive[1][2]  };
			outerCellDiffSY_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[1][0], _SinnerCellDiff[i][j][1]+primitive[1][1], _SinnerCellDiff[i][j][2]+primitive[1][2] };
		}
	}

	// 0 0 -
	auto& inner2Outer_Z = _innerCondition[_Z].nNeigh2Other;
	auto& inner2OuterS_Z = _innerCondition[S_Z].nNeigh2Other;
	auto& neighDiff_Z = _innerCondition[_Z].neighDiff;
	auto& neighDiffS_Z = _innerCondition[S_Z].neighDiff;
	inner2Outer_Z.assign( basis.size(), 0 );
	inner2OuterS_Z.assign( basis.size(), 0 );
	neighDiff_Z.resize( basis.size() );
	neighDiffS_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiff_Z[i] );
			inner2Outer_Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_Z[i] );

			inner2OuterS_Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffS_Z[i] );
			inner2OuterS_Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_Z[i] );
		}
	auto& outerCellDiff_Z = _innerCondition[_Z]._outerCellDiff;
	auto& outerCellDiffS_Z = _innerCondition[S_Z]._outerCellDiff;
	outerCellDiff_Z.resize( basis.size() );
	outerCellDiffS_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_Z[i].resize( basis.size() );
		outerCellDiffS_Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_Z[i][j] = { _innerCellDiff[i][j][0]-primitive[2][0], _innerCellDiff[i][j][1]-primitive[2][1], _innerCellDiff[i][j][2]-primitive[2][2]  };
			outerCellDiffS_Z[i][j] = { _SinnerCellDiff[i][j][0]-primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[2][2] };
		}
	}

	// 0 0 +
	auto& inner2OuterZ_ = _innerCondition[Z_].nNeigh2Other;
	auto& inner2OuterSZ_ = _innerCondition[SZ_].nNeigh2Other;
	auto& neighDiffZ_ = _innerCondition[Z_].neighDiff;
	auto& neighDiffSZ_ = _innerCondition[SZ_].neighDiff;
	inner2OuterZ_.assign( basis.size(), 0 );
	inner2OuterSZ_.assign( basis.size(), 0 );
	neighDiffZ_.resize( basis.size() );
	neighDiffSZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterZ_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffZ_[i] );
			inner2OuterZ_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffZ_[i] );

			inner2OuterSZ_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSZ_[i] );
			inner2OuterSZ_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSZ_[i] );
		}
	auto& outerCellDiffZ_ = _innerCondition[Z_]._outerCellDiff;
	auto& outerCellDiffSZ_ = _innerCondition[SZ_]._outerCellDiff;
	outerCellDiffZ_.resize( basis.size() );
	outerCellDiffSZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffZ_[i].resize( basis.size() );
		outerCellDiffSZ_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffZ_[i][j] = { _innerCellDiff[i][j][0]+primitive[2][0], _innerCellDiff[i][j][1]+primitive[2][1], _innerCellDiff[i][j][2]+primitive[2][2]  };
			outerCellDiffSZ_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[2][2] };
		}
	}

	// - - 0
	auto& inner2Outer_X_Y = _innerCondition[_X_Y].nNeigh2Other;
	auto& inner2OuterS_X_Y = _innerCondition[S_X_Y].nNeigh2Other;
	auto& neighDiff_X_Y = _innerCondition[_X_Y].neighDiff;
	auto& neighDiffS_X_Y = _innerCondition[S_X_Y].neighDiff;
	inner2Outer_X_Y.assign( basis.size(), 0 );
	inner2OuterS_X_Y.assign( basis.size(), 0 );
	neighDiff_X_Y.resize( basis.size() );
	neighDiffS_X_Y.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_X_Y[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiff_X_Y[i] );
			inner2Outer_X_Y[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_X_Y[i] );

			inner2OuterS_X_Y[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffS_X_Y[i] );
			inner2OuterS_X_Y[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_X_Y[i] );
		}
	auto& outerCellDiff_X_Y = _innerCondition[_X_Y]._outerCellDiff;
	auto& outerCellDiffS_X_Y = _innerCondition[S_X_Y]._outerCellDiff;
	outerCellDiff_X_Y.resize( basis.size() );
	outerCellDiffS_X_Y.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_X_Y[i].resize( basis.size() );
		outerCellDiffS_X_Y[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_X_Y[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0], _innerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1], _innerCellDiff[i][j][2]-primitive[0][2]-primitive[1][2]  };
			outerCellDiffS_X_Y[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0], _SinnerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1], _SinnerCellDiff[i][j][2]-primitive[0][2]-primitive[1][2] };
		}
	}

	// - + 0
	auto& inner2Outer_XY_ = _innerCondition[_XY_].nNeigh2Other;
	auto& inner2OuterS_XY_ = _innerCondition[S_XY_].nNeigh2Other;
	auto& neighDiff_XY_ = _innerCondition[_XY_].neighDiff;
	auto& neighDiffS_XY_ = _innerCondition[S_XY_].neighDiff;
	inner2Outer_XY_.assign( basis.size(), 0 );
	inner2OuterS_XY_.assign( basis.size(), 0 );
	neighDiff_XY_.resize( basis.size() );
	neighDiffS_XY_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_XY_[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiff_XY_[i] );
			inner2Outer_XY_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiff_XY_[i] );

			inner2OuterS_XY_[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffS_XY_[i] );
			inner2OuterS_XY_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffS_XY_[i] );
		}
	auto& outerCellDiff_XY_ = _innerCondition[_XY_]._outerCellDiff;
	auto& outerCellDiffS_XY_ = _innerCondition[S_XY_]._outerCellDiff;
	outerCellDiff_XY_.resize( basis.size() );
	outerCellDiffS_XY_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_XY_[i].resize( basis.size() );
		outerCellDiffS_XY_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_XY_[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0], _innerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1], _innerCellDiff[i][j][2]-primitive[0][2]+primitive[1][2]  };
			outerCellDiffS_XY_[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0], _SinnerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1], _SinnerCellDiff[i][j][2]-primitive[0][2]+primitive[1][2] };
		}
	}

	// + - 0
	auto& inner2OuterX__Y = _innerCondition[X__Y].nNeigh2Other;
	auto& inner2OuterSX__Y = _innerCondition[SX__Y].nNeigh2Other;
	auto& neighDiffX__Y = _innerCondition[X__Y].neighDiff;
	auto& neighDiffSX__Y = _innerCondition[SX__Y].neighDiff;
	inner2OuterX__Y.assign( basis.size(), 0 );
	inner2OuterSX__Y.assign( basis.size(), 0 );
	neighDiffX__Y.resize( basis.size() );
	neighDiffSX__Y.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX__Y[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX__Y[i] );
			inner2OuterX__Y[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX__Y[i] );

			inner2OuterSX__Y[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX__Y[i] );
			inner2OuterSX__Y[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX__Y[i] );
		}
	auto& outerCellDiffX__Y = _innerCondition[X__Y]._outerCellDiff;
	auto& outerCellDiffSX__Y = _innerCondition[SX__Y]._outerCellDiff;
	outerCellDiffX__Y.resize( basis.size() );
	outerCellDiffSX__Y.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX__Y[i].resize( basis.size() );
		outerCellDiffSX__Y[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX__Y[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0], _innerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1], _innerCellDiff[i][j][2]+primitive[0][2]-primitive[1][2]  };
			outerCellDiffSX__Y[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0], _SinnerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1], _SinnerCellDiff[i][j][2]+primitive[0][2]-primitive[1][2] };
		}
	}

	// + + 0
	auto& inner2OuterX_Y_ = _innerCondition[X_Y_].nNeigh2Other;
	auto& inner2OuterSX_Y_ = _innerCondition[SX_Y_].nNeigh2Other;
	auto& neighDiffX_Y_ = _innerCondition[X_Y_].neighDiff;
	auto& neighDiffSX_Y_ = _innerCondition[SX_Y_].neighDiff;
	inner2OuterX_Y_.assign( basis.size(), 0 );
	inner2OuterSX_Y_.assign( basis.size(), 0 );
	neighDiffX_Y_.resize( basis.size() );
	neighDiffSX_Y_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX_Y_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX_Y_[i] );
			inner2OuterX_Y_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX_Y_[i] );

			inner2OuterSX_Y_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX_Y_[i] );
			inner2OuterSX_Y_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX_Y_[i] );
		}
	auto& outerCellDiffX_Y_ = _innerCondition[X_Y_]._outerCellDiff;
	auto& outerCellDiffSX_Y_ = _innerCondition[SX_Y_]._outerCellDiff;
	outerCellDiffX_Y_.resize( basis.size() );
	outerCellDiffSX_Y_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX_Y_[i].resize( basis.size() );
		outerCellDiffSX_Y_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX_Y_[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0], _innerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1], _innerCellDiff[i][j][2]+primitive[0][2]+primitive[1][2]  };
			outerCellDiffSX_Y_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0], _SinnerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1], _SinnerCellDiff[i][j][2]+primitive[0][2]+primitive[1][2] };
		}
	}

	// - 0 -
	auto& inner2Outer_X_Z = _innerCondition[_X_Z].nNeigh2Other;
	auto& inner2OuterS_X_Z = _innerCondition[S_X_Z].nNeigh2Other;
	auto& neighDiff_X_Z = _innerCondition[_X_Z].neighDiff;
	auto& neighDiffS_X_Z = _innerCondition[S_X_Z].neighDiff;
	inner2Outer_X_Z.assign( basis.size(), 0 );
	inner2OuterS_X_Z.assign( basis.size(), 0 );
	neighDiff_X_Z.resize( basis.size() );
	neighDiffS_X_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_X_Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiff_X_Z[i] );
			inner2Outer_X_Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_X_Z[i] );

			inner2OuterS_X_Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
			inner2OuterS_X_Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_X_Z[i] );
		}
	auto& outerCellDiff_X_Z = _innerCondition[_X_Z]._outerCellDiff;
	auto& outerCellDiffS_X_Z = _innerCondition[S_X_Z]._outerCellDiff;
	outerCellDiff_X_Z.resize( basis.size() );
	outerCellDiffS_X_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_X_Z[i].resize( basis.size() );
		outerCellDiffS_X_Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_X_Z[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]-primitive[2][0], _innerCellDiff[i][j][1]-primitive[0][1]-primitive[2][1], _innerCellDiff[i][j][2]-primitive[0][2]-primitive[2][2]  };
			outerCellDiffS_X_Z[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]-primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[0][1]-primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[0][2]-primitive[2][2] };
		}
	}

	// - 0 +
	auto& inner2Outer_XZ_ = _innerCondition[_XZ_].nNeigh2Other;
	auto& inner2OuterS_XZ_ = _innerCondition[S_XZ_].nNeigh2Other;
	auto& neighDiff_XZ_ = _innerCondition[_XZ_].neighDiff;
	auto& neighDiffS_XZ_ = _innerCondition[S_XZ_].neighDiff;
	inner2Outer_XZ_.assign( basis.size(), 0 );
	inner2OuterS_XZ_.assign( basis.size(), 0 );
	neighDiff_XZ_.resize( basis.size() );
	neighDiffS_XZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_XZ_[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiff_XZ_[i] );
			inner2Outer_XZ_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiff_XZ_[i] );

			inner2OuterS_XZ_[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
			inner2OuterS_XZ_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffS_XZ_[i] );
		}
	auto& outerCellDiff_XZ_ = _innerCondition[_XZ_]._outerCellDiff;
	auto& outerCellDiffS_XZ_ = _innerCondition[S_XZ_]._outerCellDiff;
	outerCellDiff_XZ_.resize( basis.size() );
	outerCellDiffS_XZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_XZ_[i].resize( basis.size() );
		outerCellDiffS_XZ_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_XZ_[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]+primitive[2][0], _innerCellDiff[i][j][1]-primitive[0][1]+primitive[2][1], _innerCellDiff[i][j][2]-primitive[0][2]+primitive[2][2]  };
			outerCellDiffS_XZ_[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]+primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[0][1]+primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[0][2]+primitive[2][2] };
		}
	}

	// + 0 -
	auto& inner2OuterX__Z = _innerCondition[X__Z].nNeigh2Other;
	auto& inner2OuterSX__Z = _innerCondition[SX__Z].nNeigh2Other;
	auto& neighDiffX__Z = _innerCondition[X__Z].neighDiff;
	auto& neighDiffSX__Z = _innerCondition[SX__Z].neighDiff;
	inner2OuterX__Z.assign( basis.size(), 0 );
	inner2OuterSX__Z.assign( basis.size(), 0 );
	neighDiffX__Z.resize( basis.size() );
	neighDiffSX__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX__Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX__Z[i] );
			inner2OuterX__Z[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX__Z[i] );

			inner2OuterSX__Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX__Z[i] );
			inner2OuterSX__Z[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX__Z[i] );
		}
	auto& outerCellDiffX__Z = _innerCondition[X__Z]._outerCellDiff;
	auto& outerCellDiffSX__Z = _innerCondition[SX__Z]._outerCellDiff;
	outerCellDiffX__Z.resize( basis.size() );
	outerCellDiffSX__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX__Z[i].resize( basis.size() );
		outerCellDiffSX__Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX__Z[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]-primitive[2][0], _innerCellDiff[i][j][1]+primitive[0][1]-primitive[2][1], _innerCellDiff[i][j][2]+primitive[0][2]-primitive[2][2]  };
			outerCellDiffSX__Z[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]-primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[0][1]-primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[0][2]-primitive[2][2] };
		}
	}

	// + 0 +
	auto& inner2OuterX_Z_ = _innerCondition[X_Z_].nNeigh2Other;
	auto& inner2OuterSX_Z_ = _innerCondition[SX_Z_].nNeigh2Other;
	auto& neighDiffX_Z_ = _innerCondition[X_Z_].neighDiff;
	auto& neighDiffSX_Z_ = _innerCondition[SX_Z_].neighDiff;
	inner2OuterX_Z_.assign( basis.size(), 0 );
	inner2OuterSX_Z_.assign( basis.size(), 0 );
	neighDiffX_Z_.resize( basis.size() );
	neighDiffSX_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX_Z_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX_Z_[i] );
			inner2OuterX_Z_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX_Z_[i] );

			inner2OuterSX_Z_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX_Z_[i] );
			inner2OuterSX_Z_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX_Z_[i] );
		}
	auto& outerCellDiffX_Z_ = _innerCondition[X_Z_]._outerCellDiff;
	auto& outerCellDiffSX_Z_ = _innerCondition[SX_Z_]._outerCellDiff;
	outerCellDiffX_Z_.resize( basis.size() );
	outerCellDiffSX_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX_Z_[i].resize( basis.size() );
		outerCellDiffSX_Z_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX_Z_[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]+primitive[2][0], _innerCellDiff[i][j][1]+primitive[0][1]+primitive[2][1], _innerCellDiff[i][j][2]+primitive[0][2]+primitive[2][2]  };
			outerCellDiffSX_Z_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]+primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[0][1]+primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[0][2]+primitive[2][2] };
		}
	}

	// 0 - -
	auto& inner2Outer_Y_Z = _innerCondition[_Y_Z].nNeigh2Other;
	auto& inner2OuterS_Y_Z = _innerCondition[S_Y_Z].nNeigh2Other;
	auto& neighDiff_Y_Z = _innerCondition[_Y_Z].neighDiff;
	auto& neighDiffS_Y_Z = _innerCondition[S_Y_Z].neighDiff;
	inner2Outer_Y_Z.assign( basis.size(), 0 );
	inner2OuterS_Y_Z.assign( basis.size(), 0 );
	neighDiff_Y_Z.resize( basis.size() );
	neighDiffS_Y_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_Y_Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_Y_Z[i] );
			inner2Outer_Y_Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_Y_Z[i] );

			inner2OuterS_Y_Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
			inner2OuterS_Y_Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_Y_Z[i] );
		}
	auto& outerCellDiff_Y_Z = _innerCondition[_Y_Z]._outerCellDiff;
	auto& outerCellDiffS_Y_Z = _innerCondition[S_Y_Z]._outerCellDiff;
	outerCellDiff_Y_Z.resize( basis.size() );
	outerCellDiffS_Y_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_Y_Z[i].resize( basis.size() );
		outerCellDiffS_Y_Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_Y_Z[i][j] = { _innerCellDiff[i][j][0]-primitive[1][0]-primitive[2][0], _innerCellDiff[i][j][1]-primitive[1][1]-primitive[2][1], _innerCellDiff[i][j][2]-primitive[1][2]-primitive[2][2]  };
			outerCellDiffS_Y_Z[i][j] = { _SinnerCellDiff[i][j][0]-primitive[1][0]-primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[1][1]-primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[1][2]-primitive[2][2] };
		}
	}

	// 0 - +
	auto& inner2Outer_YZ_ = _innerCondition[_YZ_].nNeigh2Other;
	auto& inner2OuterS_YZ_ = _innerCondition[S_YZ_].nNeigh2Other;
	auto& neighDiff_YZ_ = _innerCondition[_YZ_].neighDiff;
	auto& neighDiffS_YZ_ = _innerCondition[S_YZ_].neighDiff;
	inner2Outer_YZ_.assign( basis.size(), 0 );
	inner2OuterS_YZ_.assign( basis.size(), 0 );
	neighDiff_YZ_.resize( basis.size() );
	neighDiffS_YZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_YZ_[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiff_YZ_[i] );
			inner2Outer_YZ_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiff_YZ_[i] );

			inner2OuterS_YZ_[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
			inner2OuterS_YZ_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffS_YZ_[i] );
		}
	auto& outerCellDiff_YZ_ = _innerCondition[_YZ_]._outerCellDiff;
	auto& outerCellDiffS_YZ_ = _innerCondition[S_YZ_]._outerCellDiff;
	outerCellDiff_YZ_.resize( basis.size() );
	outerCellDiffS_YZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_YZ_[i].resize( basis.size() );
		outerCellDiffS_YZ_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_YZ_[i][j] = { _innerCellDiff[i][j][0]-primitive[1][0]+primitive[2][0], _innerCellDiff[i][j][1]-primitive[1][1]+primitive[2][1], _innerCellDiff[i][j][2]-primitive[1][2]+primitive[2][2]  };
			outerCellDiffS_YZ_[i][j] = { _SinnerCellDiff[i][j][0]-primitive[1][0]+primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[1][1]+primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[1][2]+primitive[2][2] };
		}
	}

	// 0 + -
	auto& inner2OuterY__Z = _innerCondition[Y__Z].nNeigh2Other;
	auto& inner2OuterSY__Z = _innerCondition[SY__Z].nNeigh2Other;
	auto& neighDiffY__Z = _innerCondition[Y__Z].neighDiff;
	auto& neighDiffSY__Z = _innerCondition[SY__Z].neighDiff;
	inner2OuterY__Z.assign( basis.size(), 0 );
	inner2OuterSY__Z.assign( basis.size(), 0 );
	neighDiffY__Z.resize( basis.size() );
	neighDiffSY__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterY__Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffY__Z[i] );
			inner2OuterY__Z[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffY__Z[i] );

			inner2OuterSY__Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSY__Z[i] );
			inner2OuterSY__Z[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSY__Z[i] );
		}
	auto& outerCellDiffY__Z = _innerCondition[Y__Z]._outerCellDiff;
	auto& outerCellDiffSY__Z = _innerCondition[SY__Z]._outerCellDiff;
	outerCellDiffY__Z.resize( basis.size() );
	outerCellDiffSY__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffY__Z[i].resize( basis.size() );
		outerCellDiffSY__Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffY__Z[i][j] = { _innerCellDiff[i][j][0]+primitive[1][0]-primitive[2][0], _innerCellDiff[i][j][1]+primitive[1][1]-primitive[2][1], _innerCellDiff[i][j][2]+primitive[1][2]-primitive[2][2]  };
			outerCellDiffSY__Z[i][j] = { _SinnerCellDiff[i][j][0]+primitive[1][0]-primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[1][1]-primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[1][2]-primitive[2][2] };
		}
	}

	// 0 + +
	auto& inner2OuterY_Z_ = _innerCondition[Y_Z_].nNeigh2Other;
	auto& inner2OuterSY_Z_ = _innerCondition[SY_Z_].nNeigh2Other;
	auto& neighDiffY_Z_ = _innerCondition[Y_Z_].neighDiff;
	auto& neighDiffSY_Z_ = _innerCondition[SY_Z_].neighDiff;
	inner2OuterY_Z_.assign( basis.size(), 0 );
	inner2OuterSY_Z_.assign( basis.size(), 0 );
	neighDiffY_Z_.resize( basis.size() );
	neighDiffSY_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterY_Z_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffY_Z_[i] );
			inner2OuterY_Z_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffY_Z_[i] );

			inner2OuterSY_Z_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSY_Z_[i] );
			inner2OuterSY_Z_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSY_Z_[i] );
		}
	auto& outerCellDiffY_Z_ = _innerCondition[Y_Z_]._outerCellDiff;
	auto& outerCellDiffSY_Z_ = _innerCondition[SY_Z_]._outerCellDiff;
	outerCellDiffY_Z_.resize( basis.size() );
	outerCellDiffSY_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffY_Z_[i].resize( basis.size() );
		outerCellDiffSY_Z_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffY_Z_[i][j] = { _innerCellDiff[i][j][0]+primitive[1][0]+primitive[2][0], _innerCellDiff[i][j][1]+primitive[1][1]+primitive[2][1], _innerCellDiff[i][j][2]+primitive[1][2]+primitive[2][2]  };
			outerCellDiffSY_Z_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[1][0]+primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[1][1]+primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[1][2]+primitive[2][2] };
		}
	}

	// - - -
	auto& inner2Outer_X_Y_Z = _innerCondition[_X_Y_Z].nNeigh2Other;
	auto& inner2OuterS_X_Y_Z = _innerCondition[S_X_Y_Z].nNeigh2Other;
	auto& neighDiff_X_Y_Z = _innerCondition[_X_Y_Z].neighDiff;
	auto& neighDiffS_X_Y_Z = _innerCondition[S_X_Y_Z].neighDiff;
	inner2Outer_X_Y_Z.assign( basis.size(), 0 );
	inner2OuterS_X_Y_Z.assign( basis.size(), 0 );
	neighDiff_X_Y_Z.resize( basis.size() );
	neighDiffS_X_Y_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_X_Y_Z[i] );
			inner2Outer_X_Y_Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_X_Y_Z[i] );

			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
			inner2OuterS_X_Y_Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_X_Y_Z[i] );
		}
	auto& outerCellDiff_X_Y_Z = _innerCondition[_X_Y_Z]._outerCellDiff;
	auto& outerCellDiffS_X_Y_Z = _innerCondition[S_X_Y_Z]._outerCellDiff;
	outerCellDiff_X_Y_Z.resize( basis.size() );
	outerCellDiffS_X_Y_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_X_Y_Z[i].resize( basis.size() );
		outerCellDiffS_X_Y_Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_X_Y_Z[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0]-primitive[2][0], _innerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1]-primitive[2][1], _innerCellDiff[i][j][2]-primitive[0][2]-primitive[1][2]-primitive[2][2]  };
			outerCellDiffS_X_Y_Z[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0]-primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1]-primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[0][2]-primitive[1][2]-primitive[2][2] };
		}
	}

	// - - +
	auto& inner2Outer_X_YZ_ = _innerCondition[_X_YZ_].nNeigh2Other;
	auto& inner2OuterS_X_YZ_ = _innerCondition[S_X_YZ_].nNeigh2Other;
	auto& neighDiff_X_YZ_ = _innerCondition[_X_YZ_].neighDiff;
	auto& neighDiffS_X_YZ_ = _innerCondition[S_X_YZ_].neighDiff;
	inner2Outer_X_YZ_.assign( basis.size(), 0 );
	inner2OuterS_X_YZ_.assign( basis.size(), 0 );
	neighDiff_X_YZ_.resize( basis.size() );
	neighDiffS_X_YZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiff_X_YZ_[i] );
			inner2Outer_X_YZ_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiff_X_YZ_[i] );

			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
			inner2OuterS_X_YZ_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffS_X_YZ_[i] );
		}
	auto& outerCellDiff_X_YZ_ = _innerCondition[_X_YZ_]._outerCellDiff;
	auto& outerCellDiffS_X_YZ_ = _innerCondition[S_X_YZ_]._outerCellDiff;
	outerCellDiff_X_YZ_.resize( basis.size() );
	outerCellDiffS_X_YZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_X_YZ_[i].resize( basis.size() );
		outerCellDiffS_X_YZ_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_X_YZ_[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0]+primitive[2][0], _innerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1]+primitive[2][1], _innerCellDiff[i][j][2]-primitive[0][2]-primitive[1][2]+primitive[2][2]  };
			outerCellDiffS_X_YZ_[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]-primitive[1][0]+primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[0][1]-primitive[1][1]+primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[0][2]-primitive[1][2]+primitive[2][2] };
		}
	}

	// - + -
	auto& inner2Outer_XY__Z = _innerCondition[_XY__Z].nNeigh2Other;
	auto& inner2OuterS_XY__Z = _innerCondition[S_XY__Z].nNeigh2Other;
	auto& neighDiff_XY__Z = _innerCondition[_XY__Z].neighDiff;
	auto& neighDiffS_XY__Z = _innerCondition[S_XY__Z].neighDiff;
	inner2Outer_XY__Z.assign( basis.size(), 0 );
	inner2OuterS_XY__Z.assign( basis.size(), 0 );
	neighDiff_XY__Z.resize( basis.size() );
	neighDiffS_XY__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_XY__Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiff_XY__Z[i] );
			inner2Outer_XY__Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiff_XY__Z[i] );

			inner2OuterS_XY__Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
			inner2OuterS_XY__Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffS_XY__Z[i] );
		}
	auto& outerCellDiff_XY__Z = _innerCondition[_XY__Z]._outerCellDiff;
	auto& outerCellDiffS_XY__Z = _innerCondition[S_XY__Z]._outerCellDiff;
	outerCellDiff_XY__Z.resize( basis.size() );
	outerCellDiffS_XY__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_XY__Z[i].resize( basis.size() );
		outerCellDiffS_XY__Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_XY__Z[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0]-primitive[2][0], _innerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1]-primitive[2][1], _innerCellDiff[i][j][2]-primitive[0][2]+primitive[1][2]-primitive[2][2]  };
			outerCellDiffS_XY__Z[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0]-primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1]-primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[0][2]+primitive[1][2]-primitive[2][2] };
		}
	}

	// - + +
	auto& inner2Outer_XY_Z_ = _innerCondition[_XY_Z_].nNeigh2Other;
	auto& inner2OuterS_XY_Z_ = _innerCondition[S_XY_Z_].nNeigh2Other;
	auto& neighDiff_XY_Z_ = _innerCondition[_XY_Z_].neighDiff;
	auto& neighDiffS_XY_Z_ = _innerCondition[S_XY_Z_].neighDiff;
	inner2Outer_XY_Z_.assign( basis.size(), 0 );
	inner2OuterS_XY_Z_.assign( basis.size(), 0 );
	neighDiff_XY_Z_.resize( basis.size() );
	neighDiffS_XY_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, 0, 0, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiff_XY_Z_[i] );
			inner2Outer_XY_Z_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiff_XY_Z_[i] );

			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, 0, 0, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
			inner2OuterS_XY_Z_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffS_XY_Z_[i] );
		}
	auto& outerCellDiff_XY_Z_ = _innerCondition[_XY_Z_]._outerCellDiff;
	auto& outerCellDiffS_XY_Z_ = _innerCondition[S_XY_Z_]._outerCellDiff;
	outerCellDiff_XY_Z_.resize( basis.size() );
	outerCellDiffS_XY_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiff_XY_Z_[i].resize( basis.size() );
		outerCellDiffS_XY_Z_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiff_XY_Z_[i][j] = { _innerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0]+primitive[2][0], _innerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1]+primitive[2][1], _innerCellDiff[i][j][2]-primitive[0][2]+primitive[1][2]+primitive[2][2]  };
			outerCellDiffS_XY_Z_[i][j] = { _SinnerCellDiff[i][j][0]-primitive[0][0]+primitive[1][0]+primitive[2][0], _SinnerCellDiff[i][j][1]-primitive[0][1]+primitive[1][1]+primitive[2][1], _SinnerCellDiff[i][j][2]-primitive[0][2]+primitive[1][2]+primitive[2][2] };
		}
	}

	// + - -
	auto& inner2OuterX__Y_Z = _innerCondition[X__Y_Z].nNeigh2Other;
	auto& inner2OuterSX__Y_Z = _innerCondition[SX__Y_Z].nNeigh2Other;
	auto& neighDiffX__Y_Z = _innerCondition[X__Y_Z].neighDiff;
	auto& neighDiffSX__Y_Z = _innerCondition[SX__Y_Z].neighDiff;
	inner2OuterX__Y_Z.assign( basis.size(), 0 );
	inner2OuterSX__Y_Z.assign( basis.size(), 0 );
	neighDiffX__Y_Z.resize( basis.size() );
	neighDiffSX__Y_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiffX__Y_Z[i] );
			inner2OuterX__Y_Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffX__Y_Z[i] );

			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
			inner2OuterSX__Y_Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSX__Y_Z[i] );
		}
	auto& outerCellDiffX__Y_Z = _innerCondition[X__Y_Z]._outerCellDiff;
	auto& outerCellDiffSX__Y_Z = _innerCondition[SX__Y_Z]._outerCellDiff;
	outerCellDiffX__Y_Z.resize( basis.size() );
	outerCellDiffSX__Y_Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX__Y_Z[i].resize( basis.size() );
		outerCellDiffSX__Y_Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX__Y_Z[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0]-primitive[2][0], _innerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1]-primitive[2][1], _innerCellDiff[i][j][2]+primitive[0][2]-primitive[1][2]-primitive[2][2]  };
			outerCellDiffSX__Y_Z[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0]-primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1]-primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[0][2]-primitive[1][2]-primitive[2][2] };
		}
	}

	// + - +
	auto& inner2OuterX__YZ_ = _innerCondition[X__YZ_].nNeigh2Other;
	auto& inner2OuterSX__YZ_ = _innerCondition[SX__YZ_].nNeigh2Other;
	auto& neighDiffX__YZ_ = _innerCondition[X__YZ_].neighDiff;
	auto& neighDiffSX__YZ_ = _innerCondition[SX__YZ_].neighDiff;
	inner2OuterX__YZ_.assign( basis.size(), 0 );
	inner2OuterSX__YZ_.assign( basis.size(), 0 );
	neighDiffX__YZ_.resize( basis.size() );
	neighDiffSX__YZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX__YZ_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 0, -1, 0, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( -1, -1, 0, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiffX__YZ_[i] );
			inner2OuterX__YZ_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffX__YZ_[i] );

			inner2OuterSX__YZ_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 0, -1, 0, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( -1, -1, 0, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
			inner2OuterSX__YZ_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSX__YZ_[i] );
		}
	auto& outerCellDiffX__YZ_ = _innerCondition[X__YZ_]._outerCellDiff;
	auto& outerCellDiffSX__YZ_ = _innerCondition[SX__YZ_]._outerCellDiff;
	outerCellDiffX__YZ_.resize( basis.size() );
	outerCellDiffSX__YZ_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX__YZ_[i].resize( basis.size() );
		outerCellDiffSX__YZ_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX__YZ_[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0]+primitive[2][0], _innerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1]+primitive[2][1], _innerCellDiff[i][j][2]+primitive[0][2]-primitive[1][2]+primitive[2][2]  };
			outerCellDiffSX__YZ_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]-primitive[1][0]+primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[0][1]-primitive[1][1]+primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[0][2]-primitive[1][2]+primitive[2][2] };
		}
	}

	// + + -
	auto& inner2OuterX_Y__Z = _innerCondition[X_Y__Z].nNeigh2Other;
	auto& inner2OuterSX_Y__Z = _innerCondition[SX_Y__Z].nNeigh2Other;
	auto& neighDiffX_Y__Z = _innerCondition[X_Y__Z].neighDiff;
	auto& neighDiffSX_Y__Z = _innerCondition[SX_Y__Z].neighDiff;
	inner2OuterX_Y__Z.assign( basis.size(), 0 );
	inner2OuterSX_Y__Z.assign( basis.size(), 0 );
	neighDiffX_Y__Z.resize( basis.size() );
	neighDiffSX_Y__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 0, 0, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( 0, -1, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( -1, 0, -1, basis, i, j, neighDiffX_Y__Z[i] );
			inner2OuterX_Y__Z[i] += CellNeighboring( -1, -1, -1, basis, i, j, neighDiffX_Y__Z[i] );

			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 0, 0, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( 0, -1, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( -1, 0, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
			inner2OuterSX_Y__Z[i] += CellNeighboring( -1, -1, -1, shiftBasis, i, j, neighDiffSX_Y__Z[i] );
		}
	auto& outerCellDiffX_Y__Z = _innerCondition[X_Y__Z]._outerCellDiff;
	auto& outerCellDiffSX_Y__Z = _innerCondition[SX_Y__Z]._outerCellDiff;
	outerCellDiffX_Y__Z.resize( basis.size() );
	outerCellDiffSX_Y__Z.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX_Y__Z[i].resize( basis.size() );
		outerCellDiffSX_Y__Z[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX_Y__Z[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0]-primitive[2][0], _innerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1]-primitive[2][1], _innerCellDiff[i][j][2]+primitive[0][2]+primitive[1][2]-primitive[2][2]  };
			outerCellDiffSX_Y__Z[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0]-primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1]-primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[0][2]+primitive[1][2]-primitive[2][2] };
		}
	}

	// + + +
	auto& inner2OuterX_Y_Z_ = _innerCondition[X_Y_Z_].nNeigh2Other;
	auto& inner2OuterSX_Y_Z_ = _innerCondition[SX_Y_Z_].nNeigh2Other;
	auto& neighDiffX_Y_Z_ = _innerCondition[X_Y_Z_].neighDiff;
	auto& neighDiffSX_Y_Z_ = _innerCondition[SX_Y_Z_].neighDiff;
	inner2OuterX_Y_Z_.assign( basis.size(), 0 );
	inner2OuterSX_Y_Z_.assign( basis.size(), 0 );
	neighDiffX_Y_Z_.resize( basis.size() );
	neighDiffSX_Y_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, -1, -1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, -1, 0, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, -1, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, 0, -1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, 0, 0, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, 0, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, 1, -1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, 1, 0, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 1, 1, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 0, 1, -1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 0, 1, 0, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 0, 1, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( -1, 1, -1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( -1, 1, 0, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( -1, 1, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 0, 0, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( 0, -1, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( -1, 0, 1, basis, i, j, neighDiffX_Y_Z_[i] );
			inner2OuterX_Y_Z_[i] += CellNeighboring( -1, -1, 1, basis, i, j, neighDiffX_Y_Z_[i] );

			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, -1, -1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, -1, 0, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, -1, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, 0, -1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, 0, 0, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, 0, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, 1, -1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, 1, 0, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 1, 1, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 0, 1, -1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 0, 1, 0, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 0, 1, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( -1, 1, -1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( -1, 1, 0, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( -1, 1, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 0, 0, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( 0, -1, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( -1, 0, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
			inner2OuterSX_Y_Z_[i] += CellNeighboring( -1, -1, 1, shiftBasis, i, j, neighDiffSX_Y_Z_[i] );
		}
	auto& outerCellDiffX_Y_Z_ = _innerCondition[X_Y_Z_]._outerCellDiff;
	auto& outerCellDiffSX_Y_Z_ = _innerCondition[SX_Y_Z_]._outerCellDiff;
	outerCellDiffX_Y_Z_.resize( basis.size() );
	outerCellDiffSX_Y_Z_.resize( basis.size() );
	for (	unsigned i=0; i<basis.size(); ++i )
	{
		outerCellDiffX_Y_Z_[i].resize( basis.size() );
		outerCellDiffSX_Y_Z_[i].resize( basis.size() );
		for (	unsigned j=0; j<basis.size(); ++j )
		{
			outerCellDiffX_Y_Z_[i][j] = { _innerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0]+primitive[2][0], _innerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1]+primitive[2][1], _innerCellDiff[i][j][2]+primitive[0][2]+primitive[1][2]+primitive[2][2]  };
			outerCellDiffSX_Y_Z_[i][j] = { _SinnerCellDiff[i][j][0]+primitive[0][0]+primitive[1][0]+primitive[2][0], _SinnerCellDiff[i][j][1]+primitive[0][1]+primitive[1][1]+primitive[2][1], _SinnerCellDiff[i][j][2]+primitive[0][2]+primitive[1][2]+primitive[2][2] };
		}
	}
}


unsigned
InterfaceCellIdentifier3DImp::CellNeighboring (
	int shiftX,
	int shiftY,
	int shiftZ,
	const vector<data3d>& basis,
	unsigned basisI,
	unsigned basisJ,
	vector<data3d>& neighDiff
) const noexcept
{
	data3d shift{};
	shift[0] = primitive[0][0]*shiftX + primitive[1][0]*shiftY + primitive[2][0]*shiftZ +
		basis[basisJ][0] - basis[basisI][0];
	shift[1] = primitive[0][1]*shiftX + primitive[1][1]*shiftY + primitive[2][1]*shiftZ +
		basis[basisJ][1] - basis[basisI][1];
	shift[2] = primitive[0][2]*shiftX + primitive[1][2]*shiftY + primitive[2][2]*shiftZ +
		basis[basisJ][2] - basis[basisI][2];
	if ( OutOfRange(shift) )
		return 0;
	neighDiff.emplace_back( move(shift) );
	return 1;
}


bool
InterfaceCellIdentifier3DImp::IsValidInnerRelation (
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
InterfaceCellIdentifier3DImp::IsValidInnerCell (
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
				if ( !regularCase && !shiftCase ) return false;
				break;
			}
		}
	}
	if( orientation == REGULAR && !regularCase )
		return false;
	else if ( orientation == SHIFT && ! shiftCase )
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

const vector<vector<data3d>>&
InterfaceCellIdentifier3DImp::OuterCellDifference (
	InterfaceCoordinate coord,
	const Cell& cell
) const noexcept
{
	switch ( coord )
	{
	case X :
		return cell.x == MINUS ?
			( cell.orientation == REGULAR ?
				_innerCondition.at(_X)._outerCellDiff :
				_innerCondition.at(S_X)._outerCellDiff ) :
			( cell.orientation == REGULAR ?
				_innerCondition.at(X_)._outerCellDiff :
				_innerCondition.at(SX_)._outerCellDiff );
	case Y :
		return cell.y == MINUS ?
			( cell.orientation == REGULAR ?
				_innerCondition.at(_Y)._outerCellDiff :
				_innerCondition.at(S_Y)._outerCellDiff ) :
			( cell.orientation == REGULAR ?
				_innerCondition.at(Y_)._outerCellDiff :
				_innerCondition.at(SY_)._outerCellDiff );
	case Z :
		return cell.z == MINUS ?
			( cell.orientation == REGULAR ?
				_innerCondition.at(_Z)._outerCellDiff :
				_innerCondition.at(S_Z)._outerCellDiff ) :
			( cell.orientation == REGULAR ?
				_innerCondition.at(Z_)._outerCellDiff :
				_innerCondition.at(SZ_)._outerCellDiff );
	case XY :
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
	case XZ :
		if ( cell.x == MINUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_X_Z)._outerCellDiff :
				_innerCondition.at(S_X_Z)._outerCellDiff;
		else if ( cell.x == MINUS && cell.z == PLUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_XZ_)._outerCellDiff :
				_innerCondition.at(S_XZ_)._outerCellDiff;
		else if ( cell.x == PLUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(X__Z)._outerCellDiff :
				_innerCondition.at(SX__Z)._outerCellDiff;
		else // cell.x == PLUS && cell.z == PLUS
			return cell.orientation == REGULAR ?
				_innerCondition.at(X_Z_)._outerCellDiff :
				_innerCondition.at(SX_Z_)._outerCellDiff;
	case YZ :
		if ( cell.y == MINUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_Y_Z)._outerCellDiff :
				_innerCondition.at(S_Y_Z)._outerCellDiff;
		else if ( cell.y == MINUS && cell.z == PLUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_YZ_)._outerCellDiff :
				_innerCondition.at(S_YZ_)._outerCellDiff;
		else if ( cell.y == PLUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(Y__Z)._outerCellDiff :
				_innerCondition.at(SY__Z)._outerCellDiff;
		else // cell.y == PLUS && cell.z == PLUS
			return cell.orientation == REGULAR ?
				_innerCondition.at(Y_Z_)._outerCellDiff :
				_innerCondition.at(SY_Z_)._outerCellDiff;
	case XYZ :
		if ( cell.x == MINUS && cell.y == MINUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_X_Y_Z)._outerCellDiff :
				_innerCondition.at(S_X_Y_Z)._outerCellDiff;
		else if ( cell.x == MINUS && cell.y == MINUS && cell.z == PLUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_X_YZ_)._outerCellDiff :
				_innerCondition.at(S_X_YZ_)._outerCellDiff;
		else if ( cell.x == MINUS && cell.y == PLUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_XY__Z)._outerCellDiff :
				_innerCondition.at(S_XY__Z)._outerCellDiff;
		else if ( cell.x == MINUS && cell.y == PLUS && cell.z == PLUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(_XY_Z_)._outerCellDiff :
				_innerCondition.at(S_XY_Z_)._outerCellDiff;
		else if ( cell.x == PLUS && cell.y == MINUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(X__Y_Z)._outerCellDiff :
				_innerCondition.at(SX__Y_Z)._outerCellDiff;
		else if ( cell.x == PLUS && cell.y == MINUS && cell.z == PLUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(X__YZ_)._outerCellDiff :
				_innerCondition.at(SX__YZ_)._outerCellDiff;
		else if ( cell.x == PLUS && cell.y == PLUS && cell.z == MINUS )
			return cell.orientation == REGULAR ?
				_innerCondition.at(X_Y__Z)._outerCellDiff :
				_innerCondition.at(SX_Y__Z)._outerCellDiff;
		else // cell.x == PLUS && cell.y == PLUS && cell.z == PLUS
			return cell.orientation == REGULAR ?
				_innerCondition.at(X_Y_Z_)._outerCellDiff :
				_innerCondition.at(SX_Y_Z_)._outerCellDiff;
	default :
		cout << "check 1"<<endl;
		return _innerCellDiff;
	}
}


void
InterfaceCellIdentifier3DImp::ChainingCells (
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
		cell.position[2] = (cellPos[2]-basePos[2]) / primitive[2][2];
	}
}


void
InterfaceCellIdentifier3DImp::CheckInnerCellChaining (
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
		bool _x = iter->x == MINUS ? true : false;
		bool x_ = iter->x == PLUS ? true : false;
		bool _y = iter->y == MINUS ? true : false;
		bool y_ = iter->y == PLUS ? true : false;
		bool _z = iter->z == MINUS ? true : false;
		bool z_ = iter->z == PLUS ? true : false;
		if ( iter->x == NONE && iter->y == NONE )
			_z = z_ = true;
		if ( iter->x == NONE && iter->z == NONE )
			_y = y_ = true;
		if ( iter->y == NONE && iter->z == NONE )
			_x = x_ = true;
		for (	auto iterJ = cells.begin();
			iterJ != cells.end();
			++iterJ
		)
		{
			if ( iter == iterJ ) continue;
			auto diff = Difference (
				data3d{ iterJ->position[0]*primitive[0][0], iterJ->position[1]*primitive[1][1], iterJ->position[2]*primitive[2][2] },
				data3d{ iter->position[0]*primitive[0][0], iter->position[1]*primitive[1][1], iter->position[2]*primitive[2][2] }
			);
			diff[0] /= primitive[0][0];
			diff[1] /= primitive[1][1];
			diff[2] /= primitive[2][2];
			if ( OutOfRange(diff) ) continue;
			if ( diff[0] < 0. ) _x = true;
			if ( diff[0] > 0. ) x_ = true;
			if ( diff[1] < 0. ) _y = true;
			if ( diff[1] > 0. ) y_ = true;
			if ( diff[2] < 0. ) _z = true;
			if ( diff[2] > 0. ) z_ = true;
		}
		if ( _x && x_ && _y && y_ && _z && z_ ) continue;

		InterfaceDirection x = NONE;
		InterfaceDirection y = NONE;
		InterfaceDirection z = NONE;
		if ( !_x ) x = MINUS;
		else if ( !x_ ) x = PLUS;
		if ( !_y ) y = MINUS;
		else if ( !y_ ) x = PLUS;
		if ( !_z ) z = MINUS;
		else if ( !z_ ) z = PLUS;
		auto interfaceType = DetailType( x, y, z, iter->orientation );
		CellType cellType = TypeConvert( interfaceType );
		unsigned index = LatticeLocatingIndex( interfaceType );// need test
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
InterfaceCellIdentifier3DImp::LatticeLocatingIndex (
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
		unsigned _z = 0;
		unsigned z_ = 0;
		switch ( orientation )
		{
		case REGULAR :
			for ( const auto& diff : _innerCellDiff[index] )
			{
				if ( diff[0] > 0. ) ++_x;
				else if ( diff[0] < 0. ) ++x_;
				if ( diff[1] > 0. ) ++_y;
				else if ( diff[1] < 0.) ++y_;
				if ( diff[2] > 0. ) ++_z;
				else if ( diff[2] < 0. ) ++z_;
			}
			break;
		case SHIFT :
			for ( const auto& diff : _SinnerCellDiff[index] )
			{
				if ( diff[0] > 0. ) ++_x;
				else if ( diff[0] < 0. ) ++x_;
				if ( diff[1] > 0. ) ++_y;
				else if ( diff[1] < 0. ) ++y_;
				if ( diff[2] > 0. ) ++_z;
				else if ( diff[2] < 0. ) ++z_;
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
		case _Z :
		case S_Z :
			if ( _z > z_ ) return index;
			break;
		case Z_ :
		case SZ_ :
			if ( z_ > _z ) return index;
			break;
		case _X_Y :
		case S_X_Y :
			if ( _x > x_ && _y > y_ ) return index;
			break;
		case _XY_ :
		case S_XY_ :
			if ( _x > x_ && y_ > y_ ) return index;
			break;
		case X__Y :
		case SX__Y :
			if ( x_ > _x && _y > y_ ) return index;
			break;
		case X_Y_ :
		case SX_Y_ :
			if ( x_ > _x && y_ > _y ) return index;
			break;
		case _X_Z :
		case S_X_Z :
			if ( _x > x_ && _z > z_ ) return index;
			break;
		case _XZ_ :
		case S_XZ_ :
			if ( _x > x_ && z_ > _z ) return index;
			break;
		case X__Z :
		case SX__Z :
			if ( x_ > _x && _z > z_ ) return index;
			break;
		case X_Z_ :
		case SX_Z_ :
			if ( x_ > _x && z_ > _z ) return index;
			break;
		case _Y_Z :
		case S_Y_Z :
			if ( _y > y_ && _z > z_ ) return index;
			break;
		case _YZ_ :
		case S_YZ_ :
			if ( _y > y_ && z_ > _z ) return index;
			break;
		case Y__Z :
		case SY__Z :
			if ( y_ > _y && _z > z_ ) return index;
			break;
		case Y_Z_ :
		case SY_Z_ :
			if ( y_ > _y && z_ > _z ) return index;
			break;
		case _X_Y_Z :
		case S_X_Y_Z :
			if ( _x > x_ && _y > y_ && _z > z_ ) return index;
			break;
		case _X_YZ_ :
		case S_X_YZ_ :
			if ( _x > x_ && _y > y_ && z_ > _z ) return index;
			break;
		case _XY__Z :
		case S_XY__Z :
			if ( _x > x_ && y_ > _y && _z > z_ ) return index;
			break;
		case _XY_Z_ :
		case S_XY_Z_ :
			if ( _x > x_ && y_ > _y && z_ > _z ) return index;
			break;
		case X__Y_Z :
		case SX__Y_Z :
			if ( x_ > _x && _y > y_ && _z > z_ ) return index;
			break;
		case X__YZ_ :
		case SX__YZ_ :
			if ( x_ > _x && _y > y_ && z_ > _z ) return index;
			break;
		case X_Y__Z :
		case SX_Y__Z :
			if ( x_ > _x && y_ > _y && _z > z_ ) return index;
			break;
		case X_Y_Z_ :
		case SX_Y_Z_ :
			if ( x_ > _x && y_ > _y && z_ > _z ) return index;
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
InterfaceCellIdentifier3DImp::SetCellInterface (
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
	case S_Z :
	case _Z :
		cell.z = MINUS;
		break;
	case SZ_ :
	case Z_ :
		cell.z = PLUS;
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
	case S_X_Z :
	case _X_Z :
		cell.x = MINUS;
		cell.z = MINUS;
		break;
	case S_XZ_ :
	case _XZ_ :
		cell.x = MINUS;
		cell.z = PLUS;
		break;
	case SX__Z :
	case X__Z :
		cell.x = PLUS;
		cell.z = MINUS;
		break;
	case SX_Z_ :
	case X_Z_ :
		cell.x = PLUS;
		cell.z = PLUS;
		break;
	case S_Y_Z :
	case _Y_Z :
		cell.y = MINUS;
		cell.z = MINUS;
		break;
	case S_YZ_ :
	case _YZ_ :
		cell.y = MINUS;
		cell.z = PLUS;
		break;
	case SY__Z :
	case Y__Z :
		cell.y = PLUS;
		cell.z = MINUS;
		break;
	case SY_Z_ :
	case Y_Z_ :
		cell.y = PLUS;
		cell.z = PLUS;
		break;
	case S_X_Y_Z :
	case _X_Y_Z :
		cell.x = MINUS;
		cell.y = MINUS;
		cell.z = MINUS;
		break;
	case S_X_YZ_ :
	case _X_YZ_ :
		cell.x = MINUS;
		cell.y = MINUS;
		cell.z = PLUS;
		break;
	case S_XY__Z :
	case _XY__Z :
		cell.x = MINUS;
		cell.y = PLUS;
		cell.z = MINUS;
		break;
	case S_XY_Z_ :
	case _XY_Z_ :
		cell.x = MINUS;
		cell.y = PLUS;
		cell.z = PLUS;
		break;
	case SX__Y_Z :
	case X__Y_Z :
		cell.x = PLUS;
		cell.y = MINUS;
		cell.z = MINUS;
		break;
	case SX__YZ_ :
	case X__YZ_ :
		cell.x = PLUS;
		cell.y = MINUS;
		cell.z = PLUS;
		break;
	case SX_Y__Z :
	case X_Y__Z :
		cell.x = PLUS;
		cell.y = PLUS;
		cell.z = MINUS;
		break;
	case SX_Y_Z_ :
	case X_Y_Z_ :
		cell.x = PLUS;
		cell.y = PLUS;
		cell.z = PLUS;
		break;
	default :
		cerr << "check 2:"<<endl;
	}
}


InterfaceCellIdentifier3DImp::InterfaceType
InterfaceCellIdentifier3DImp::DetailType (
	InterfaceDirection x,
	InterfaceDirection y,
	InterfaceDirection z,
	CellOrientation orientation
) const noexcept
{
	if ( orientation == ANY )
	{
		cout << "check---"<<endl;
		return NTYPES;
	}

	switch ( x )
	{
	case MINUS :
		switch ( y )
		{
		case NONE :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? _X : S_X;
			case MINUS :
				return orientation == REGULAR ? _X_Z : S_X_Z;
			case PLUS :
				return orientation == REGULAR ? _XZ_ : S_XZ_;
			}
		case MINUS :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? _X_Y : S_X_Y;
			case MINUS :
				return orientation == REGULAR ? _X_Y_Z : S_X_Y_Z;
			case PLUS :
				return orientation == REGULAR ? _X_YZ_ : S_X_YZ_;
			}
		case PLUS :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? _XY_ : S_XY_;
			case MINUS :
				return orientation == REGULAR ? _XY__Z : S_XY__Z;
			case PLUS :
				return orientation == REGULAR ? _XY_Z_ : S_XY_Z_;
			}
		}
	case PLUS :
		switch ( y )
		{
		case NONE :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? X_ : SX_;
			case MINUS :
				return orientation == REGULAR ? X__Z : SX__Z;
			case PLUS :
				return orientation == REGULAR ? X_Z_ : SX_Z_;
			}
		case MINUS :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? X__Y : SX__Y;
			case MINUS :
				return orientation == REGULAR ? X__Y_Z : SX__Y_Z;
			case PLUS :
				return orientation == REGULAR ? X__YZ_ : SX__YZ_;
			}
		case PLUS :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? X_Y_ : SX_Y_;
			case MINUS :
				return orientation == REGULAR ? X_Y__Z : SX_Y__Z;
			case PLUS :
				return orientation == REGULAR ? X_Y_Z_ : SX_Y_Z_;
			}
		}
	case NONE :
		switch ( y )
		{
		case NONE :
			switch ( z )
			{
			case NONE :
				cout << "check ... " << endl;
				return NTYPES;
			case MINUS :
				return orientation == REGULAR ? _Z : S_Z;
			case PLUS :
				return orientation == REGULAR ? Z_ : SZ_;
			}
		case MINUS :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? _Y : S_Y;
			case MINUS :
				return orientation == REGULAR ? _Y_Z : S_Y_Z;
			case PLUS :
				return orientation == REGULAR ? _YZ_ : S_YZ_;
			}
		case PLUS :
			switch ( z )
			{
			case NONE :
				return orientation == REGULAR ? Y_ : SY_;
			case MINUS :
				return orientation == REGULAR ? Y__Z : SY__Z;
			case PLUS :
				return orientation == REGULAR ? Y_Z_ : SY_Z_;
			}
		}
	default :
		cout << "check..." << endl;
		return NTYPES;
	}
}


InterfaceCellIdentifier3DImp::InterfaceType
InterfaceCellIdentifier3DImp::DetailType (
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
InterfaceCellIdentifier3DImp::TypeConvert (
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
	case _Z :
	case S_Z :
	case Z_ :
	case SZ_ :
		return FACE;
	case _X_Y :
	case _XY_ :
	case X__Y :
	case X_Y_ :
	case S_X_Y :
	case S_XY_ :
	case SX__Y :
	case SX_Y_ :
	case _X_Z :
	case _XZ_ :
	case X__Z :
	case X_Z_ :
	case S_X_Z :
	case S_XZ_ :
	case SX__Z :
	case SX_Z_ :
	case _Y_Z :
	case _YZ_ :
	case Y__Z :
	case Y_Z_ :
	case S_Y_Z :
	case S_YZ_ :
	case SY__Z :
	case SY_Z_ :
		return EDGE;
	case _X_Y_Z :
	case _X_YZ_ :
	case _XY__Z :
	case _XY_Z_ :
	case X__Y_Z :
	case X__YZ_ :
	case X_Y__Z :
	case X_Y_Z_ :
	case S_X_Y_Z :
	case S_X_YZ_ :
	case S_XY__Z :
	case S_XY_Z_ :
	case SX__Y_Z :
	case SX__YZ_ :
	case SX_Y__Z :
	case SX_Y_Z_ :
		return CORNER;
	default :
		cout << "check 3"<<endl;
		return CellType::NTYPES;
	}
}


CellOrientation
InterfaceCellIdentifier3DImp::Orientation (
	InterfaceType interfaceType
) const noexcept
{
	switch ( interfaceType )
	{
	case _X :
	case X_ :
	case _Y :
	case Y_ :
	case _Z :
	case Z_ :
	case _X_Y :
	case _XY_ :
	case _X_Z :
	case _XZ_ :
	case _Y_Z :
	case _YZ_ :
	case X__Y :
	case X_Y_ :
	case X__Z :
	case X_Z_ :
	case Y__Z :
	case Y_Z_ :
	case _X_Y_Z:
	case _X_YZ_:
	case _XY__Z:
	case _XY_Z_:
	case X__Y_Z:
	case X__YZ_:
	case X_Y__Z:
	case X_Y_Z_:
		return REGULAR;
	case S_X :
	case SX_ :
	case S_Y :
	case SY_ :
	case S_Z :
	case SZ_ :
	case S_X_Y :
	case S_XY_ :
	case S_X_Z :
	case S_XZ_ :
	case S_Y_Z :
	case S_YZ_ :
	case SX__Y :
	case SX_Y_ :
	case SX__Z :
	case SX_Z_ :
	case SY__Z :
	case SY_Z_ :
	case S_X_Y_Z:
	case S_X_YZ_:
	case S_XY__Z:
	case S_XY_Z_:
	case SX__Y_Z:
	case SX__YZ_:
	case SX_Y__Z:
	case SX_Y_Z_:
		return SHIFT;
	default :
		cout << "check 4 " <<endl;
		return ANY;
	}
}


data3d
InterfaceCellIdentifier3DImp::CellBasePosition (
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
		for (const auto& diff : _innerCellDiff[index] )
			Increment( sum, diff );
		break;
	case SHIFT :
		for ( const auto& diff : _SinnerCellDiff[index] )
			Increment( sum, diff );
		break;
	default :
		return {};
	}

	auto nBasis = _innerCellDiff[index].size();
	auto inverse = 1./nBasis;
	sum[0] *= inverse;
	sum[1] *= inverse;
	sum[2] *= inverse;
	Increment( sum, pos[ cell.id[index] ] );
	return sum;
}


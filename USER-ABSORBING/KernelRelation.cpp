#include "KernelRelation.h"
#include <cstdlib>
#include <utility>

using namespace std;


KernelRelationImp::CellInfoPack::CellInfoPack (
	const std::vector<Cell>& outerCells
) noexcept :
	periodicStatus( {true, true, true} )
{
	lowerBound[0] = upperBound[0] = outerCells[0].position[0];
	lowerBound[1] = upperBound[1] = outerCells[0].position[1];
	lowerBound[2] = upperBound[2] = outerCells[0].position[2];
	for ( const auto& cell : outerCells )
	{
		if ( lowerBound[0] > cell.position[0] ) lowerBound[0] = cell.position[0];
		if ( upperBound[0] < cell.position[0] ) upperBound[0] = cell.position[0];
		if ( lowerBound[1] > cell.position[1] ) lowerBound[1] = cell.position[1];
		if ( upperBound[1] < cell.position[1] ) upperBound[1] = cell.position[1];
		if ( lowerBound[2] > cell.position[2] ) lowerBound[2] = cell.position[2];
		if ( upperBound[2] < cell.position[2] ) upperBound[2] = cell.position[2];

		if ( periodicStatus[0] && cell.x!=NONE ) periodicStatus[0] = false;
		if ( periodicStatus[1] && cell.y!=NONE ) periodicStatus[1] = false;
		if ( periodicStatus[2] && cell.z!=NONE ) periodicStatus[2] = false;
	}
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

KernelRelationImp::KernelRelationImp (
	unsigned validInnerRange,
	unsigned validOuterRange,
	const KernelFunctionsMap& kernelFunction,
	CellShiftFunc func,
	LatticeIndexRotationFunc rFunc
) noexcept :
	_validInnerIDRange ( validInnerRange ),
	_validOuterIDRange ( validOuterRange ),
	_kernelFunction ( kernelFunction ),
	CellShift ( func ),
	LatticeIndexRotation ( rFunc )
{}


void
KernelRelationImp::operator() (
	const vector<Cell>& innerCells,
	const vector<Cell>& outerCells,
	std::vector<RelationKernelFuncPack>& relationKernelFunc
) const noexcept
{
	CellInfoPack cellInfoPack( outerCells );
	for ( const auto& innerCell : innerCells )
		for ( const auto& outerCell : outerCells )
			InterRelation( innerCell, outerCell, cellInfoPack, relationKernelFunc );
}


void
KernelRelationImp::InterRelation (
	const Cell& innerCell,
	const Cell& outerCell,
	const CellInfoPack& cellInfoPack,
	std::vector<RelationKernelFuncPack>& relationKernelFunc
) const noexcept
{
	for ( unsigned j=0; j<outerCell.id.size(); ++j )
	{
		auto outerID = outerCell.id[j];
		if ( !IsValidOuterID(outerID) ) continue;
		for ( unsigned i=0; i<innerCell.id.size(); ++i )
		{
			auto innerID = innerCell.id[i];
			if ( !IsValidInnerID(innerID) ) continue;
			auto config = InterCellRelation( innerCell, i, outerCell, j, cellInfoPack );
			if ( !IsKernelFunctionExist(config) )
			{
				if ( !ConvertEdge2FaceRelation(config) )///need test
					continue;
				if ( !IsKernelFunctionExist(config) )
					continue;
			}
			AddPairConfig( relationKernelFunc, outerID, innerID, move(config) );
		}
	}
}


InterfacePairConfiguration
KernelRelationImp::InterCellRelation (
	const Cell& innerCell,
	unsigned innerLocal,
	const Cell& outerCell,
	unsigned outerLocal,
	const CellInfoPack& cellInfoPack
) const noexcept
{
	auto innerCellPos = innerCell.position;
	auto outerCellPos = outerCell.position;
	EdgeRelation( innerCellPos, outerCellPos, cellInfoPack );
	auto config = OrientationRelation (
		innerCell, innerLocal, move(innerCellPos),
		outerCell, outerLocal, move(outerCellPos)
	);
	SignFix( config, cellInfoPack );
	PeriodicRelation( config, cellInfoPack );
	InterfacePairConfigurationTransformer()( config, CellShift, LatticeIndexRotation );
	return config;
}


void
KernelRelationImp::EdgeRelation (
	data3d& innerCellPos,
	data3d& outerCellPos,
	const CellInfoPack& cellInfoPack
) const noexcept
{
	auto nearEdgePos = NearEdgeCellPos( outerCellPos, cellInfoPack );
	innerCellPos[0] -= nearEdgePos[0];
	innerCellPos[1] -= nearEdgePos[1];
	innerCellPos[2] -= nearEdgePos[2];
	outerCellPos[0] -= nearEdgePos[0];
	outerCellPos[1] -= nearEdgePos[1];
	outerCellPos[2] -= nearEdgePos[2];
}


InterfacePairConfiguration
KernelRelationImp::OrientationRelation (
	const Cell& innerCell,
	unsigned innerLocal,
	data3d&& innerCellPos,
	const Cell& outerCell,
	unsigned outerLocal,
	data3d&& outerCellPos
) const noexcept
{
	ShiftCellOrientation( innerCell, innerCellPos, innerLocal );
	auto shiftTimes = ShiftCellOrientation( outerCell, outerCellPos, outerLocal );

	InterfacePairConfiguration config;
	config.orientation = outerCell.orientation;
	if ( shiftTimes %2 != 0 )
		config.orientation = config.orientation==REGULAR ? SHIFT : REGULAR;
	config.outerLocal = outerLocal;
	config.innerLocal = innerLocal;
	config.outerCellID[0] = round( outerCellPos[0] );
	config.outerCellID[1] = round( outerCellPos[1] );
	config.outerCellID[2] = round( outerCellPos[2] );
	config.innerCellID[0] = round( innerCellPos[0] );
	config.innerCellID[1] = round( innerCellPos[1] );
	config.innerCellID[2] = round( innerCellPos[2] );
	return config;
}


unsigned
KernelRelationImp::ShiftCellOrientation (
	const Cell& cell,
	data3d& cellPos,
	unsigned& local
) const noexcept
{
	unsigned shiftTimes = 0;
	auto orientation = cell.orientation;
	for ( unsigned coord=0; coord<3; ++coord )
		if ( fabs( cellPos[coord] - round(cellPos[coord]) ) > 1.e-5 )
		{
			++shiftTimes;
			CellShift( orientation, cellPos, local, coord );
		}
	return shiftTimes;
}


void
KernelRelationImp::SignFix (
	InterfacePairConfiguration& config,
	const CellInfoPack& cellInfoPack
) const noexcept
{
	for ( unsigned coord=0; coord<3; ++coord )
		if ( !cellInfoPack.periodicStatus[coord] && config.innerCellID[coord] < 0 ) 
		{
			config.outerCellID[coord] = -config.outerCellID[coord];
			config.innerCellID[coord] = -config.innerCellID[coord];
			config.sign[coord] = true;
		}
}


void
KernelRelationImp::PeriodicRelation (
	InterfacePairConfiguration& config,
	const CellInfoPack& cellInfoPack
) const noexcept
{
	for ( unsigned coord=0; coord<3; ++coord )
	{
		if ( !cellInfoPack.periodicStatus[coord] ) continue;
		AliasConfiguration( config, coord );

		auto width = round( cellInfoPack.upperBound[coord] - cellInfoPack.lowerBound[coord] + 1. );
		auto& diff = config.innerCellID[coord];
		if ( abs(diff) > 0.5*width )
			diff += signbit(diff) ? width : -width;
		config.outerCellID[coord] = -1;
	}
}

//---------------------------------------------------------------------------//

bool
KernelRelationImp::AliasConfiguration (
	InterfacePairConfiguration& config,
	unsigned coord
) const noexcept
{
	config.innerCellID[coord] -= config.outerCellID[coord];
	config.outerCellID[coord] = 0;

	return true;
}


data3d
KernelRelationImp::NearEdgeCellPos (
	const data3d& cell,
	const CellInfoPack& cellInfoPack
) const noexcept
{
	const auto& lowerBound = cellInfoPack.lowerBound;
	const auto& upperBound = cellInfoPack.upperBound;
	data3d nearEdge{};
	for ( unsigned coord=0; coord<3; ++coord )
		nearEdge[coord] = cellInfoPack.periodicStatus[coord] ? cell[coord] : (
			abs(cell[coord]-lowerBound[coord]) < abs(cell[coord]-upperBound[coord]) ? lowerBound[coord] : upperBound[coord]
		);
	return nearEdge;
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

bool
KernelRelationImp::IsKernelFunctionExist (
	const InterfacePairConfiguration& config
) const noexcept
{
//return true;
	return _kernelFunction.find( config ) != _kernelFunction.end();
}


void
KernelRelationImp::AddPairConfig (
	std::vector<RelationKernelFuncPack>& relationKernelFunc,
	unsigned outerID,
	unsigned innerID,
	InterfacePairConfiguration&& config
) const noexcept
{
	KernelFuncPack pack{forward<InterfacePairConfiguration>(config), nullptr};
	pack.funcPtr = _kernelFunction.at( pack.config );
	relationKernelFunc[outerID].emplace (
		make_pair( innerID, move(pack) )
	);
}


bool
KernelRelationImp::ConvertEdge2FaceRelation (
	InterfacePairConfiguration& config
) const noexcept
{
	bool converted = false;
	for ( unsigned coord=0; coord<3; ++coord )
		if ( config.outerCellID[coord] >= 0 && config.innerCellID[coord] > 1 )
		{
			config.innerCellID[coord] -= config.outerCellID[coord];
			config.outerCellID[coord] = -1;
			InterfacePairConfigurationTransformer()( config, CellShift, LatticeIndexRotation );
			ShiftConfigOrientationToRegular( config, coord );
			converted = true;
		}
	return converted;
}


void
KernelRelationImp::ShiftConfigOrientationToRegular (
	InterfacePairConfiguration& config,
	unsigned coord
) const noexcept
{
	if ( config.orientation == REGULAR ) return;

	auto orientation = SHIFT;
	data3d outerPos;
	CellShift(orientation, outerPos, config.outerLocal, coord);

	data3d innerPos;
	CellShift(config.orientation, innerPos, config.innerLocal, coord);

	if ( outerPos[coord] * innerPos[coord] < 0 )
	{
		if ( outerPos[coord] > 0 )
			--config.innerCellID[coord];
		else
			++config.innerCellID[coord];
	}
}


bool
KernelRelationImp::IsValidInnerID (
	unsigned id
) const noexcept
{
	return id < _validInnerIDRange;
}


bool
KernelRelationImp::IsValidOuterID (
	unsigned id
) const noexcept
{
	return id < _validOuterIDRange;
}


#include "KernelRelationComponent.h"
#include <string>
#include <iostream>

using namespace std;

bool
InterfacePairConfiguration::operator== (
	const InterfacePairConfiguration& other
) const noexcept
{
	return	orientation == other.orientation &&
		outerCellID == other.outerCellID &&
		outerLocal == other.outerLocal &&
		innerCellID == other.innerCellID &&
		innerLocal == other.innerLocal;
}


unsigned
InterfacePairConfiguration::ConfigureDistanceSq (
	const data3d& cellSize
) const noexcept
{
	array<int,3> difference;
	difference[0] = cellSize[0] * ( outerCellID[0]<0 ? innerCellID[0] : innerCellID[0]-outerCellID[0] );
	difference[1] = cellSize[1] * ( outerCellID[1]<0 ? innerCellID[1] : innerCellID[1]-outerCellID[1] );
	difference[2] = cellSize[2] * ( outerCellID[2]<0 ? innerCellID[2] : innerCellID[2]-outerCellID[2] );
	return difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2];
}


CellType
InterfacePairConfiguration::Type () const noexcept
{
	unsigned pbcCount = 0;
	pbcCount += outerCellID[0]==-1 ? 1 : 0;
	pbcCount += outerCellID[1]==-1 ? 1 : 0;
	pbcCount += outerCellID[2]==-1 ? 1 : 0;

	switch (pbcCount)
	{
	case 0:
		return CORNER;
	case 1:
		return EDGE;
	default: //case 2:
		return FACE;
	}
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

void
InterfacePairConfigurationTransformer::operator() (
	InterfacePairConfiguration& config,
	CellShiftFunc_t CellShiftFunc,
	LatticeIndexRotationFunc_t LatticeIndexRotationFunc
) const noexcept
{
	constexpr int cases[3] = {1, 2, 4};
	int direction = 0;
	int face = 0;
	for ( unsigned coord=0; coord<3; ++coord )
	{
		if ( config.outerCellID[coord] >= 0 )
			face += cases[coord];
		if ( config.innerCellID[coord] < 0 )
			direction += cases[coord];
	}

	if ( face==3 || face>=5 ) // 3 5 6 7
	{
		direction = 0;
		for ( int coord=0; coord<3; ++coord )
		{
			if ( config.sign[coord] || (config.outerCellID[coord]==0 && config.innerCellID[coord]<0) )
				direction += cases[coord];
			if ( config.sign[coord] )
			{
				config.innerCellID[coord] = -config.innerCellID[coord];
				config.sign[coord] = false;
			}
		}
	}
	int type = face*10 + direction;
//cout<<type<<endl;
	switch ( type )
	{
	case 11: // self - self
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		if ( config.innerCellID[1] < 0 )
			Mirror( config, CellShiftFunc, 1 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 1 );
		break;
	case 22:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		if ( config.innerCellID[0] < 0 )
			Mirror( config, CellShiftFunc, 0 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 0 );
		break;
	case 44:
		Rotate180( config, LatticeIndexRotationFunc, 1 );
		if ( config.innerCellID[0] < 0 )
			Mirror( config, CellShiftFunc, 0 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 0 );
	case 33:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		break;
	case 55:
		Rotate180( config, LatticeIndexRotationFunc, 1 );
		break;
	case 66:
		Rotate180( config, LatticeIndexRotationFunc, 0 );
		break;
	case 77:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		Mirror( config, CellShiftFunc, 2 );
		break;

	case 13:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 0 );
		break;
	case 23:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 1 );
		break;
	case 15:
		Rotate180( config, LatticeIndexRotationFunc, 1 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 0 );
		break;
	case 45:
		Rotate180( config, LatticeIndexRotationFunc, 1 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 2 );
		break;
	case 26:
		Rotate180( config, LatticeIndexRotationFunc, 0 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 1 );
		break;
	case 46:
		Rotate180( config, LatticeIndexRotationFunc, 0 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 2 );
		break;

	case 21:
	case 41:
		Mirror( config, CellShiftFunc, 0 );
		CellShift( config, CellShiftFunc, 0 );
		break;
	case 12:
	case 42:
		Mirror( config, CellShiftFunc, 1 );
		CellShift( config, CellShiftFunc, 1 );
		break;
	case 14:
	case 24:
		Mirror( config, CellShiftFunc, 2 );
		CellShift( config, CellShiftFunc, 2 );
		break;
	case 16:
		Mirror( config, CellShiftFunc, 1 );
		CellShift( config, CellShiftFunc, 1 );
		Mirror( config, CellShiftFunc, 2 );
		CellShift( config, CellShiftFunc, 2 );
		break;
	case 25:
		Mirror( config, CellShiftFunc, 0 );
		CellShift( config, CellShiftFunc, 0 );
		Mirror( config, CellShiftFunc, 2 );
		CellShift( config, CellShiftFunc, 2 );
		break;
	case 43 :
		Mirror( config, CellShiftFunc, 0 );
		CellShift( config, CellShiftFunc, 0 );
		Mirror( config, CellShiftFunc, 1 );
		CellShift( config, CellShiftFunc, 1 );
		break;

	case 31:
	case 51:
	case 61:
	case 71:
		Mirror( config, CellShiftFunc, 0 );
		break;
	case 32:
	case 52:
	case 62:
	case 72:
		Mirror( config, CellShiftFunc, 1 );
		break;
	case 53:
	case 63:
	case 73:
		Mirror( config, CellShiftFunc, 0 );
		Mirror( config, CellShiftFunc, 1 );
		break;
	case 34:
	case 54:
	case 64:
	case 74:
		Mirror( config, CellShiftFunc, 2 );
		break;
	case 35:
	case 65:
	case 75:
		Mirror( config, CellShiftFunc, 0 );
		Mirror( config, CellShiftFunc, 2 );
		break;
	case 36:
	case 56:
	case 76:
		Mirror( config, CellShiftFunc, 1 );
		Mirror( config, CellShiftFunc, 2 );
		break;
	case 37:
	case 57:
	case 67:
		Mirror( config, CellShiftFunc, 0 );
		Mirror( config, CellShiftFunc, 1 );
		Mirror( config, CellShiftFunc, 2 );
		break;

	case 17:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		Mirror( config, CellShiftFunc, 2 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 0 );
		break;
	case 27:
		Rotate180( config, LatticeIndexRotationFunc, 2 );
		Mirror( config, CellShiftFunc, 2 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 1 );
		break;
	case 47:
		Rotate180( config, LatticeIndexRotationFunc, 0 );
		Mirror( config, CellShiftFunc, 0 );
		if ( config.orientation == SHIFT )
			CellShift( config, CellShiftFunc, 2 );
		break;

	default:
		break;
	}
}


void
InterfacePairConfigurationTransformer::Rotate180 (
	InterfacePairConfiguration& config,
	LatticeIndexRotationFunc_t LatticeIndexRotationFunc,
	unsigned coord
) const noexcept
{
	for ( unsigned i=0; i<3; ++i )
	{
		if ( i == coord ) continue;
		config.innerCellID[i] = -config.innerCellID[i];
		config.rotateMatrix[i][i] *= -1.;
	}

	char aCoord = coord == 2 ? 'z' : ( coord == 1 ? 'y' : 'x' );
	LatticeIndexRotationFunc( config.outerLocal, aCoord, 180 );
	auto orientationShift = LatticeIndexRotationFunc( config.innerLocal, aCoord, 180 );

	if ( orientationShift )
		config.orientation = config.orientation == REGULAR ? SHIFT : REGULAR;
}


void
InterfacePairConfigurationTransformer::Mirror (
	InterfacePairConfiguration& config,
	CellShiftFunc_t CellShiftFunc,
	unsigned coord
) const noexcept
{
	data3d innerPos;
	CellOrientation innerOrientation = config.orientation;
	CellShiftFunc( innerOrientation, innerPos, config.innerLocal, coord );

	data3d outerPos;
	CellShiftFunc( config.orientation, outerPos, config.outerLocal, coord );

	config.innerCellID[coord] = -config.innerCellID[coord];
	config.rotateMatrix[coord][coord] *= -1.;
}


void
InterfacePairConfigurationTransformer::CellShift (
	InterfacePairConfiguration& config,
	CellShiftFunc_t CellShiftFunc,
	unsigned coord
) const noexcept
{
	data3d innerPos;
	CellOrientation innerOrientation = config.orientation;
	CellShiftFunc( innerOrientation, innerPos, config.innerLocal, coord );

	data3d outerPos;
	CellShiftFunc( config.orientation, outerPos, config.outerLocal, coord );

	if ( innerPos[coord] * outerPos[coord] < 0 )
	{
		if ( outerPos[coord] > 0 )
			--config.innerCellID[coord];
		else
			++config.innerCellID[coord];
	}
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

std::size_t
InterfacePairConfigurationHasher::operator() (
	const InterfacePairConfiguration& config
) const noexcept
{
	char str[35];
	sprintf ( str, "%d-%d_%d_%d_%d-%d_%d_%d_%d",
		int(config.orientation),
		config.outerCellID[0], config.outerCellID[1], config.outerCellID[2],
		config.outerLocal,
		config.innerCellID[0], config.innerCellID[1], config.innerCellID[2],
		config.innerLocal
	);
	return hash<string>()( str );
}


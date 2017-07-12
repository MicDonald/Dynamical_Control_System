#ifndef CELL_H_INCLUDED
#define CELL_H_INCLUDED

#include "data3d.h"
#include <vector>

enum InterfaceCoordinate {X, Y, Z, XY, XZ, YZ, XYZ};
enum CellType {FACE, EDGE, CORNER, NTYPES};
enum CellOrientation {REGULAR, SHIFT, ANY};
enum InterfaceDirection {PLUS, MINUS, NONE};


struct Cell {
	data3d position = {};
	CellOrientation orientation = ANY;
	std::vector<unsigned> id;
	InterfaceDirection x = NONE;
	InterfaceDirection y = NONE;
	InterfaceDirection z = NONE;
};


auto CellInterfaceAt = [] (
		const Cell& cell,
		InterfaceCoordinate coord
	) noexcept
	{
		switch ( coord )
		{
		case X:
			return cell.x != NONE;
		case Y:
			return cell.y != NONE;
		case Z:
			return cell.z != NONE;
		case XY:
			return (cell.x != NONE) && (cell.y != NONE);
		case XZ:
			return (cell.x != NONE) && (cell.z != NONE);
		case YZ:
			return (cell.y != NONE) && (cell.z != NONE);
		default: // XYZ
			return (cell.x != NONE) && (cell.y != NONE) && (cell.z != NONE);
		}
	};


auto ReciprocalCell = [] (
		const Cell& cell,
		InterfaceCoordinate coord,
		unsigned defaultID = 0
	) noexcept
	{
		Cell reciprocalCell;
		reciprocalCell.orientation = cell.orientation;
		reciprocalCell.id.assign( cell.id.size(), defaultID );
		switch ( coord )
		{
		case X:
			reciprocalCell.x = cell.x==PLUS ? MINUS : PLUS;
			break;
		case Y:
			reciprocalCell.y = cell.y==PLUS ? MINUS : PLUS;
			break;
		case Z:
			reciprocalCell.z = cell.z==PLUS ? MINUS : PLUS;
			break;
		case XY:
			reciprocalCell.x = cell.x==PLUS ? MINUS : PLUS;
			reciprocalCell.y = cell.y==PLUS ? MINUS : PLUS;
			break;
		case XZ:
			reciprocalCell.x = cell.x==PLUS ? MINUS : PLUS;
			reciprocalCell.z = cell.z==PLUS ? MINUS : PLUS;
			break;
		case YZ:
			reciprocalCell.y = cell.y==PLUS ? MINUS : PLUS;
			reciprocalCell.z = cell.z==PLUS ? MINUS : PLUS;
			break;
		default: //XYZ
			reciprocalCell.x = cell.x==PLUS ? MINUS : PLUS;
			reciprocalCell.y = cell.y==PLUS ? MINUS : PLUS;
			reciprocalCell.z = cell.z==PLUS ? MINUS : PLUS;
		}
		return reciprocalCell;
	};



#endif // CELL_H_INCLUDED

#ifndef KERNELRELATIONCOMPONENT_H_INCLUDED
#define KERNELRELATIONCOMPONENT_H_INCLUDED

#include "Cell.h"
#include "data3d.h"
#include <type_traits>
#include <array>
#include <unordered_map>
#include <memory>


struct InterfacePairConfiguration {
	CellOrientation orientation = ANY;
	std::array<int,3> outerCellID = {};
	unsigned outerLocal = 0;
	std::array<int,3> innerCellID = {};
	unsigned innerLocal = 0;
	std::array<bool,3> sign = {};// true = negative
	double rotateMatrix[3][3] = { {1.}, {0, 1.}, {0, 0, 1.} };


	bool
	operator== (
		const InterfacePairConfiguration& other
	) const noexcept;


	unsigned
	ConfigureDistanceSq (
		const data3d& cellSize
	) const noexcept;


	CellType Type () const noexcept;
};

//---------------------------------------------------------------------------//

class InterfacePairConfigurationTransformer {
public:
	using CellShiftFunc_t = void (*)( CellOrientation&, data3d&, unsigned&, unsigned );
	using LatticeIndexRotationFunc_t = bool (*)( unsigned&, char, int );


	void
	operator() (
		InterfacePairConfiguration& config,
		CellShiftFunc_t,
		LatticeIndexRotationFunc_t
	) const noexcept;

private:
	void
	Rotate180 (
		InterfacePairConfiguration& config,
		LatticeIndexRotationFunc_t,
		unsigned coord
	) const noexcept;


	void
	Mirror (
		InterfacePairConfiguration& config,
		CellShiftFunc_t,
		unsigned coord
	) const noexcept;


	void 
	CellShift (
		InterfacePairConfiguration& config,
		CellShiftFunc_t,
		unsigned coord
	) const noexcept;
};

//---------------------------------------------------------------------------//

struct InterfacePairConfigurationHasher {
	std::size_t
	operator() (
		const InterfacePairConfiguration& config
	) const noexcept;
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

struct KernelFunctions {
	std::unique_ptr<double[]> xx = nullptr;
	std::unique_ptr<double[]> xy = nullptr;
	std::unique_ptr<double[]> yx = nullptr;
	std::unique_ptr<double[]> yy = nullptr;
	std::unique_ptr<double[]> zz = nullptr;
	std::unique_ptr<double[]> xz = nullptr;
	std::unique_ptr<double[]> zx = nullptr;
	std::unique_ptr<double[]> yz = nullptr;
	std::unique_ptr<double[]> zy = nullptr;
};


using KernelFunctionsMap = std::unordered_map <
		InterfacePairConfiguration,
		std::shared_ptr<KernelFunctions>,
		InterfacePairConfigurationHasher
	>;

#endif // KERNELRELATIONCOMPONENT_H_INCLUDED

#ifndef KERNEL_RELATION_H_INCLUDED
#define KERNEL_RELATION_H_INCLUDED

#include "InterfaceCellIdentifier.h"
#include "KernelFileReader.h"
#include "KernelRelationComponent.h"
#include "data3d.h"
#include "Lattice.h"
#include "PBCDifference.h"
#include <vector>
#include <map>
#include <algorithm>


class KernelRelationImp {
public:
	struct KernelFuncPack {
		InterfacePairConfiguration config;
		std::shared_ptr<KernelFunctions> funcPtr;
	};
	using RelationKernelFuncPack = std::map<unsigned, KernelFuncPack>;
	using CellShiftFunc = void (*)(CellOrientation&, data3d&, unsigned&, unsigned);
	using LatticeIndexRotationFunc = bool (*)( unsigned&, char, int );

private:
	unsigned _validInnerIDRange;
	unsigned _validOuterIDRange;
	const KernelFunctionsMap& _kernelFunction;
	CellShiftFunc CellShift;
	LatticeIndexRotationFunc LatticeIndexRotation;

public:
	struct CellInfoPack {
		std::array<bool, 3> periodicStatus;
		data3d lowerBound;
		data3d upperBound;


		CellInfoPack (
			const std::vector<Cell>& outerCells
		) noexcept;
	};

	//-------------------------------------//

	KernelRelationImp (
		unsigned validInnerRange,
		unsigned validOuterRange,
		const KernelFunctionsMap& kernelFunction,
		CellShiftFunc func,
		LatticeIndexRotationFunc Func
	) noexcept;


	void
	operator() (
		const std::vector<Cell>& innerCells,
		const std::vector<Cell>& outerCells,
		std::vector<RelationKernelFuncPack>& relationKernelFunc
	) const noexcept;

private:
	void
	InterRelation (
		const Cell& innerCell,
		const Cell& outerCell,
		const CellInfoPack& cellInfoPack,
		std::vector<RelationKernelFuncPack>& relationKernelFunc
	) const noexcept;


	InterfacePairConfiguration
	InterCellRelation (
		const Cell& innerCell,
		unsigned innerLocal,
		const Cell& outerCell,
		unsigned outerLocal,
		const CellInfoPack& cellInfoPack
	) const noexcept;


	void
	EdgeRelation (
		data3d& innerCellPos,
		data3d& outerCellPos,
		const CellInfoPack& cellInfoPack
	) const noexcept;


	InterfacePairConfiguration
	OrientationRelation (
		const Cell& innerCell,
		unsigned innerLocal,
		data3d&& innerCellPos,
		const Cell& outerCell,
		unsigned outerLocal,
		data3d&& outerCellPos
	) const noexcept;


	unsigned
	ShiftCellOrientation (
		const Cell& cell,
		data3d& cellPos,
		unsigned& local
	) const noexcept;


	void
	SignFix (
		InterfacePairConfiguration& config,
		const CellInfoPack& cellInfoPack
	) const noexcept;


	void
	PeriodicRelation (
		InterfacePairConfiguration& config,
		const CellInfoPack& cellInfoPack
	) const noexcept;

	//-------------------------------------//

	data3d
	NearEdgeCellPos (
		const data3d& cellPos,
		const CellInfoPack& cellInfoPack
	) const noexcept;


	bool
	AliasConfiguration (
		InterfacePairConfiguration& config,
		unsigned coord
	) const noexcept;

	//-------------------------------------//

	bool
	IsKernelFunctionExist (
		const InterfacePairConfiguration& config
	) const noexcept;


	void
	AddPairConfig (
		std::vector<RelationKernelFuncPack>& relationKernelFunc,
		unsigned outerID,
		unsigned innerID,
		InterfacePairConfiguration&& config
	) const noexcept;


	bool
	ConvertEdge2FaceRelation (
		InterfacePairConfiguration& config
	) const noexcept;


	void
	ShiftConfigOrientationToRegular (
		InterfacePairConfiguration& config,
		unsigned coord
	) const noexcept;


	bool
	IsValidInnerID (
		unsigned id
	) const noexcept;


	bool
	IsValidOuterID (
		unsigned id
	) const noexcept;
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

class KernelRelationBase {
protected:
	using RelationKernelFuncPack = typename KernelRelationImp::RelationKernelFuncPack;
	KernelFunctionsMap _kernelFunction;
	std::vector<RelationKernelFuncPack> _relationKernelFunc;
	const PBCDifference* _diffFunc;

public:
	virtual ~KernelRelationBase () = default;


	virtual void
	SetInitialState (
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos,
		const PBCDifference* diffFunc
	) noexcept = 0;


	decltype(auto)
	KernelFunctionsRelationOf (
		unsigned outerID
	) const noexcept
	{
		return _relationKernelFunc[outerID];
	}


	const KernelFunctions&
	KernelFunctionsOf (
		const InterfacePairConfiguration& configuration
	) const noexcept
	{
		return *_kernelFunction.at( configuration );
	}
};

//---------------------------------------------------------------------------//

template < typename BASE, typename DERIVATIVE >
class KernelRelation_CRTP: public BASE {
	DERIVATIVE&
	MyType () noexcept
	{
		return static_cast<DERIVATIVE&>(*this);
	}

public:
	void
	SetInitialState (
		const std::vector<data3d>& innerPos,
		const std::vector<data3d>& outerPos,
		const PBCDifference* diffFunc
	) noexcept
	{
		BASE::_diffFunc = diffFunc;
		auto cellizeIdentifier = MyType().CellizeIdentifier( diffFunc );
		std::cout<<"CellizeIdentifier\n";
		auto cells = cellizeIdentifier->Identify( innerPos, outerPos );
		std::cout<<"Cellize\n";
		BASE::_relationKernelFunc.resize( outerPos.size() );
		std::cout<<"Kernel Resizing...\n";
		const auto& innerCells = cells.first;
		const auto& outerCells = cells.second;
		KernelRelationImp _Imp (
			innerPos.size(),
			outerPos.size(),
			BASE::_kernelFunction,
			&DERIVATIVE::ShiftFunc,
			&DERIVATIVE::LatticeIndexRotation
		);
		std::cout<<"Imp\n";
		_Imp( innerCells, outerCells, BASE::_relationKernelFunc );
		std::cout<<"after Imp\n";
		RemoveUnusedKernels();
	}

protected:
	void
	RemoveUnusedKernels () noexcept
	{
		for (	auto iter = BASE::_kernelFunction.begin();
			iter != BASE::_kernelFunction.end();
			++iter
		)
			if ( iter->second.use_count() == 1 )
				iter->second.reset();
	}
};

//---------------------------------------------------------------------------//

template < typename LATTICE >
class KernelRelation : public KernelRelation_CRTP<KernelRelationBase, KernelRelation<LATTICE>> {
public:
	KernelRelation (
		double tCut,
		double rCut,
		double cornerCut,
		double timeStep
	) noexcept
	{
		KernelFileReader<LATTICE> reader;
		reader.ReadKernelFile( tCut, timeStep, KernelRelationBase::_kernelFunction );
		RemoveOutOfRangeKernels( rCut*rCut );
		RemoveOutOfRangeCornerKernels( cornerCut*cornerCut );
	}


	std::unique_ptr<InterfaceCellIdentifier<LATTICE>>
	CellizeIdentifier (
		const PBCDifference* diffFunc
	) noexcept
	{
		std::unique_ptr<InterfaceCellIdentifier<LATTICE>> ptr (
			new InterfaceCellIdentifier<LATTICE>( diffFunc )
		);
		return ptr;
	}

	//-------------------------------------//

	static void
	ShiftFunc (
		CellOrientation& orientation,
		data3d& cellPos,
		unsigned& local,
		unsigned coord
	) noexcept
	{
		const auto& basis = LATTICE::basis;
		auto shiftedBasis = ShiftLattice()( LATTICE::primitive, LATTICE::shift, LATTICE::basis );
		auto halfCellDist = 0.5*LATTICE::primitive[coord][coord];
		auto pos = orientation == REGULAR ? basis[local] : shiftedBasis[local];
		if ( pos[coord] >= halfCellDist - 1.e-5 )
			cellPos[coord] += 0.5;
		else
			cellPos[coord] -= 0.5;
		pos[coord] += halfCellDist;
		if ( pos[coord] >= LATTICE::primitive[coord][coord] )
			pos[coord] -= LATTICE::primitive[coord][coord];

		switch ( orientation )
		{
		case REGULAR :
			orientation = SHIFT;
			for ( unsigned i=0; i<LATTICE::N; ++i )
				if ( Equivalent( shiftedBasis[i], pos ) )
				{
					local = i;
					break;
				}
			break;
		case SHIFT :
			orientation = REGULAR;
			for ( unsigned i=0; i<LATTICE::N; ++i )
				if ( Equivalent( basis[i], pos ) )
				{
					local = i;
					break;
				}
			break;
		default :
			break;
		}
	}


	static bool
	LatticeIndexRotation (
		unsigned& index,
		char axis,
		int degree
	) noexcept
	{
		return LATTICE::Rotate( index, axis, degree );
	}

private:
	void
	RemoveOutOfRangeKernels (
		double rCutSq
	) noexcept
	{
		auto& _kernelFunction = KernelRelationBase::_kernelFunction;
		data3d cellSize = {
			LATTICE::primitive[0][0],
			LATTICE::primitive[1][1],
			LATTICE::primitive[2][2]
		};
		for (	auto iter = _kernelFunction.begin();
			iter != _kernelFunction.end();
		)
			if ( iter->first.ConfigureDistanceSq(cellSize) > rCutSq )
				iter = _kernelFunction.erase( iter );
			else
				++iter;
	}


	void
	RemoveOutOfRangeCornerKernels (
		double rCutSq
	) noexcept
	{
		auto& _kernelFunction = KernelRelationBase::_kernelFunction;
		data3d cellSize = {
			LATTICE::primitive[0][0],
			LATTICE::primitive[1][1],
			LATTICE::primitive[2][2]
		};
		for (	auto iter = _kernelFunction.begin();
			iter != _kernelFunction.end();
		)
			if ( iter->first.Type() != FACE && iter->first.ConfigureDistanceSq(cellSize) > rCutSq )
				iter = _kernelFunction.erase( iter );
			else
				++iter;
	}
};

#endif // KERNEL_RELATION_H_INCLUDED

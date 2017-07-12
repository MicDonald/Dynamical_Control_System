#include "AbsorbingBoundaryCondition.h"
#include "PositionPredictor.h"
#include "Lattice.h"
#include <algorithm>
#include <vector>
#include <iostream>

using namespace std;


void
AbsorbingBoundaryCondition::Setup (
	string latticeType,
	double timestep,
	double tCut,
	double rCut,
	double cornerCut
) noexcept
{
	_nStep = round( tCut/timestep );
	auto type = String2LatticeType( latticeType );
	SetupKernelRelation( type, tCut, rCut, cornerCut, timestep );
}


void
AbsorbingBoundaryCondition::SetupInitialState (
	const vector<data3d>& outerPos,
	const vector<data3d>& innerPos,
	double latticeSpacing,
	const PBCDifference* diffFunc
) noexcept
{
	if ( outerPos.empty() || innerPos.empty() ) return;

	_diffFunc = diffFunc;
	//std::cout<<"_diffFunc = diffFunc\n";
	_currentOuterPos.resize( outerPos.size() );
	//std::cout<<"_currentOuterPos.resize( outerPos.size() )\n";
	BackupInitialPos( innerPos, outerPos );
	//std::cout<<"BackupInitialPos( innerPos, outerPos )\n";
	ReserveHistoryDisplacementMemory();
	//std::cout<<"ReserveHistoryDisplacementMemory()\n";
	KernelRelationPairing( latticeSpacing );
	std::cout<<"KernelRelationPairing( latticeSpacing )\n";
}



void
AbsorbingBoundaryCondition::ReserveHistoryDisplacementMemory () noexcept
{
	_historyInnerDisplacement.resize( _initialInnerPos.size() );
	for ( auto& historyInnerDisplacement : _historyInnerDisplacement )
		historyInnerDisplacement.assign( _nStep, {0.} );
}


void
AbsorbingBoundaryCondition::BackupInitialPos (
	const vector<data3d>& innerPos,
	const vector<data3d>& outerPos
) noexcept
{
	_initialInnerPos = innerPos;
	_initialOuterPos = outerPos;
}


void
AbsorbingBoundaryCondition::BackupCurrentInnerPos (
	const vector<data3d>& innerPos
) noexcept
{
	auto posIter = innerPos.begin();
	auto initPosIter = _initialInnerPos.begin();
	for ( auto& historyDisp : _historyInnerDisplacement )
	{
		rotate (
			historyDisp.begin(), historyDisp.end()-1,
			historyDisp.end()
		);
		historyDisp[0] = Substract(*posIter, *initPosIter);
		_diffFunc->Correction( historyDisp[0].data() );
		++posIter;
		++initPosIter;
	}
}


void
AbsorbingBoundaryCondition::SetupKernelRelation (
	LatticeType latticeType,
	double tCut,
	double rCut,
	double cornerCut,
	double timestep
) noexcept
{
	switch ( latticeType )
	{
	case TRIANGULAR :
		_kernelRelationPtr.reset( new KernelRelation<Lattice_Triangular>(tCut, rCut, cornerCut, timestep) );
		break;
//	case HEXAGONAL :
//		_kernelRelationPtr.reset( new KernelRelation<Lattice_Hexagonal>(tCut, rCut, cornerCut, timestep) );
//		break;
	case FCC :
		_kernelRelationPtr.reset( new KernelRelation<Lattice_FCC>(tCut, rCut, cornerCut, timestep) );
		break;
	case Si :
		_kernelRelationPtr.reset( new KernelRelation<Lattice_Si>(tCut, rCut, cornerCut, timestep) );
		break;
	case SQUARE :
	        _kernelRelationPtr.reset( new KernelRelation<Lattice_Square>(tCut, rCut, cornerCut, timestep) );
                cout<<"lattice: square"<<endl;
		break;
	default :
		_kernelRelationPtr.reset( nullptr );
		cout<<"no this type"<<endl;
	}
}


void
AbsorbingBoundaryCondition::KernelRelationPairing (
	double latticeSpacing
) noexcept
{
	vector<data3d> innerPos( _initialInnerPos );
	vector<data3d> outerPos( _initialOuterPos );
	std::cout<<"Initialize position...\n";
	double reverseSpacing = 1./latticeSpacing;
	for ( auto& pos : innerPos )
	{
		pos[0] *= reverseSpacing;
		pos[1] *= reverseSpacing;
		pos[2] *= reverseSpacing;
	}
	for ( auto& pos : outerPos )
	{
		pos[0] *= reverseSpacing;
		pos[1] *= reverseSpacing;
		pos[2] *= reverseSpacing;
	}
	std::cout<<"pos/spacing...\n";
	_kernelRelationPtr->SetInitialState( innerPos, outerPos, _diffFunc );
	std::cout<<"SetInitialState\n";
for(unsigned i=0; i<_initialOuterPos.size(); ++i)
{
	const auto& relation = _kernelRelationPtr->KernelFunctionsRelationOf(i);
	std::cout<<"Finding relation...\n";
	for(const auto& r : relation )
	{
		const auto& config = r.second.config;
	cout<<i<<" "<<r.first<<": ["<<config.orientation<<"]"<<config.outerCellID[0]<<"_"<<config.outerCellID[1]<<"_"<<config.outerCellID[2]<<"_"<<config.outerLocal<<" : "<<config.innerCellID[0]<<"_"<<config.innerCellID[1]<<"_"<<config.innerCellID[2]<<"_"<<config.innerLocal<<" ("<<config.sign[0]<<" "<<config.sign[1]<<" "<<config.sign[2]<<")"<<config.rotateMatrix[0][0]<<" "<<config.rotateMatrix[1][1]<<" "<<config.rotateMatrix[2][2]<<endl;
	}
}
}


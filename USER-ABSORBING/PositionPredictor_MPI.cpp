#include "PositionPredictor.h"
#include "mpi.h"
#include <iterator>
#include <numeric>
#include <functional>

using namespace std;

vector<data3d>
PositionPredictor::operator() (
	const vector<data3d>& initialOuterPosition,
	const KernelRelationBase& kernelRelation,
	const vector<vector<data3d>>& historyInnerDisplacement,
	MPI_t
) const noexcept
{
	int nProcs;
	int myProc;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
	vector<int> processorSize(nProcs);
	vector<int> shift(nProcs, 0);
	auto nAtoms = initialOuterPosition.size();
	for (	int i=1; i<nProcs; ++i )
	{
		shift[i] = nAtoms * i / nProcs;
		processorSize[i-1] = shift[i]-shift[i-1];
	}
	processorSize.back() = nAtoms - shift.back();

	auto& startIndex = shift[myProc];
	auto beginIter = initialOuterPosition.begin() + startIndex;
	vector<data3d> myPosition( beginIter, beginIter+processorSize[myProc] );
	for (	auto currentPosIter = myPosition.begin();
		currentPosIter != myPosition.end();
		++currentPosIter
	)
	{
		auto outerID = distance( myPosition.begin(), currentPosIter ) + startIndex;
		auto kernelFunctionsRelation = kernelRelation.KernelFunctionsRelationOf( outerID );
		for (	auto iter = kernelFunctionsRelation.begin();
			iter != kernelFunctionsRelation.end();
			++iter
		)
		{
			const auto& innerID = iter->first;
			const auto& config = iter->second.config;
			const auto& kernelFunctions = *iter->second.funcPtr;
			auto xy = config.rotateMatrix[0][0] * config.rotateMatrix[1][1];
			auto xz = config.rotateMatrix[0][0] * config.rotateMatrix[2][2];
			auto yz = config.rotateMatrix[1][1] * config.rotateMatrix[2][2];
			if ( kernelFunctions.xx != nullptr )
				(*currentPosIter)[0] += ComponentConvolution( kernelFunctions.xx.get(), historyInnerDisplacement[innerID], 0 );
			if ( kernelFunctions.yx != nullptr )
				(*currentPosIter)[1] += xy*ComponentConvolution( kernelFunctions.yx.get(), historyInnerDisplacement[innerID], 0 );
			if ( kernelFunctions.zx != nullptr )
				(*currentPosIter)[2] += xz*ComponentConvolution( kernelFunctions.zx.get(), historyInnerDisplacement[innerID], 0 );
			if ( kernelFunctions.xy != nullptr )
				(*currentPosIter)[0] += xy*ComponentConvolution( kernelFunctions.xy.get(), historyInnerDisplacement[innerID], 1 );
			if ( kernelFunctions.yy != nullptr )
				(*currentPosIter)[1] += ComponentConvolution( kernelFunctions.yy.get(), historyInnerDisplacement[innerID], 1 );
			if ( kernelFunctions.zy != nullptr )
				(*currentPosIter)[2] += yz*ComponentConvolution( kernelFunctions.zy.get(), historyInnerDisplacement[innerID], 1 );
			if ( kernelFunctions.xz != nullptr )
				(*currentPosIter)[0] += xz*ComponentConvolution( kernelFunctions.xz.get(), historyInnerDisplacement[innerID], 2 );
			if ( kernelFunctions.yz != nullptr )
				(*currentPosIter)[1] += yz*ComponentConvolution( kernelFunctions.yz.get(), historyInnerDisplacement[innerID], 2 );
			if ( kernelFunctions.zz != nullptr )
				(*currentPosIter)[2] += ComponentConvolution( kernelFunctions.zz.get(), historyInnerDisplacement[innerID], 2 );
		}
	}

	for (	int i=0; i<nProcs; ++i )
	{
		processorSize[i] *= 3;
		shift[i] *= 3;
	}
	vector<data3d> currentPosition( nAtoms );
	MPI_Allgatherv(&myPosition[0], processorSize[myProc], MPI_DOUBLE, &currentPosition[0], &processorSize[0], &shift[0], MPI_DOUBLE, MPI_COMM_WORLD);
	return currentPosition;
}


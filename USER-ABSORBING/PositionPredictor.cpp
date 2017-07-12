#include "PositionPredictor.h"
#include <iterator>
#include <numeric>
#include <functional>

using namespace std;


vector<data3d>
PositionPredictor::operator() (
	const vector<data3d>& initialOuterPosition,
	const KernelRelationBase& kernelRelation,
	const vector<vector<data3d>>& historyInnerDisplacement,
	Serial_t
) const noexcept
{
	vector<data3d> currentPosition( initialOuterPosition );
	for (	auto currentPosIter = currentPosition.begin();
		currentPosIter != currentPosition.end();
		++currentPosIter
	)
	{
		auto outerID = distance( currentPosition.begin(), currentPosIter );
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
	return currentPosition;
}


double
PositionPredictor::ComponentConvolution (
	const double kernelFunction[],
	const vector<data3d>& historyInnerDisplacement,
	unsigned coord
) const noexcept
{
	return inner_product (
		historyInnerDisplacement.begin(), historyInnerDisplacement.end(),
		kernelFunction,
		0.,
		plus<double>(),
		[coord] (
			const data3d& a,
			const double& b
		) noexcept
		{
			return a[coord] * b;
		}
	);
}


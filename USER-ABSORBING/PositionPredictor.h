#ifndef POSITIONPREDICTOR_H_INCLUDED
#define POSITIONPREDICTOR_H_INCLUDED

#include "data3d.h"
#include "KernelRelation.h"
#include <vector>


struct Serial_t {};
struct MPI_t {};
struct OpenMP_t {};

class PositionPredictor {
public:
	std::vector<data3d>
	operator() (
		const std::vector<data3d>& initialOuterPosition,
		const KernelRelationBase& kernelRelation,
		const std::vector<std::vector<data3d>>& historyInnerDisplacement,
		Serial_t
	) const noexcept;


	std::vector<data3d>
	operator() (
		const std::vector<data3d>& initialOuterPosition,
		const KernelRelationBase& kernelRelation,
		const std::vector<std::vector<data3d>>& historyInnerDisplacement,
		MPI_t
	) const noexcept;


	std::vector<data3d>
	operator() (
		const std::vector<data3d>& initialOuterPosition,
		const KernelRelationBase& kernelRelation,
		const std::vector<std::vector<data3d>>& historyInnerDisplacement,
		OpenMP_t
	) const noexcept;

private:
	double
	ComponentConvolution (
		const double kernelFunction[],
		const std::vector<data3d>& historyInnerDisplacement,
		unsigned coord
	) const noexcept;


	double
	ComponentConvolution_OpenMP (
		const double kernelFunction[],
		const std::vector<data3d>& historyInnerDisplacement,
		unsigned coord
	) const noexcept;
};


#endif // POSITIONPREDICTOR_H_INCLUDED

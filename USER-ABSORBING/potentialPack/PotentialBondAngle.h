#ifndef POTENTIALBONDANGLE_H_INCLUDED
#define POTENTIALBONDANGLE_H_INCLUDED

#include "PotentialPair.h"
#include "PotentialAngle.h"

class PotentialBondAngle : public PotentialPair, public PotentialAngle {
public:
	double
	Energy (
		const array3d& rij
	) const noexcept override;


	double
	Energy (
		const array3d& rij,
		const array3d& rjk
	) const noexcept override;

	//-------------------------------------//

	array3d
	Force (
		const array3d& rij
	) const noexcept override;


	array3d
	Force (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept override;


	std::array<array3d, 2>
	Force (
		const array3d& rij,
		const array3d& rik
	) const noexcept override;


	std::array<array3d, 2>
	Force (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept override;

	//-------------------------------------//

	matrix3d
	Hessian (
		const array3d& rij
	) const noexcept override;


	matrix3d
	Hessian (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept override;


	std::array<matrix3d, 3>
	Hessian (
		const array3d& rij,
		const array3d& rik
	) const noexcept override;


	std::array<matrix3d, 3>
	Hessian (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept override;
};

#endif // POTENTIALBONDANGLE_H_INCLUDED

#include "PotentialBondAngle.h"

double
PotentialBondAngle::Energy (
	const array3d& rij
) const noexcept
{
	return PotentialPair::Energy(rij);
}


double
PotentialBondAngle::Energy (
	const array3d& rij,
	const array3d& rjk
) const noexcept
{
	return PotentialAngle::Energy(rij, rjk);
}


array3d
PotentialBondAngle::Force (
	const array3d& rij
) const noexcept
{
	return PotentialPair::Force(rij);
}


array3d
PotentialBondAngle::Force (
	const array3d& rij,
	FiniteDifference_t fd
) const noexcept
{
	return PotentialPair::Force(rij, fd);
}


std::array<array3d, 2>
PotentialBondAngle::Force (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	return PotentialAngle::Force(rij, rik);
}


std::array<array3d, 2>
PotentialBondAngle::Force (
	const array3d& rij,
	const array3d& rik,
	FiniteDifference_t fd
) const noexcept
{
	return PotentialAngle::Force(rij, rik, fd);
}


matrix3d
PotentialBondAngle::Hessian (
	const array3d& rij
) const noexcept
{
	return PotentialPair::Hessian(rij);
}


matrix3d
PotentialBondAngle::Hessian (
	const array3d& rij,
	FiniteDifference_t fd
) const noexcept
{
	return PotentialPair::Hessian(rij, fd);
}


std::array<matrix3d, 3>
PotentialBondAngle::Hessian (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	return PotentialAngle::Hessian(rij, rik);
}


std::array<matrix3d, 3>
PotentialBondAngle::Hessian (
	const array3d& rij,
	const array3d& rik,
	FiniteDifference_t fd
) const noexcept
{
	return PotentialAngle::Hessian(rij, rik, fd);
}


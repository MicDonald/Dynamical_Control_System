#include "Potential_LennardJones.h"

Potential_LennardJones::Potential_LennardJones (
	double epsilon,
	double sigma
) noexcept :
	m_epsilon( epsilon ),
	m_sigma( sigma )
{}


double
Potential_LennardJones::ObjectiveFunction (
	double rij
) const noexcept
{
	auto _rij = 1. / rij;
	auto factor = Factor6(_rij);
	return 4.*m_epsilon * factor * ( factor - 1. );
}


PotentialPair::_1stDerivative_t
Potential_LennardJones::_1stDerivative (
	double rij
) const noexcept
{
	auto _rij = 1. / rij;
	auto factor = Factor6(_rij);
	return -m_epsilon * _rij * factor * ( 48.*factor - 24. );
}


PotentialPair::_2ndDerivative_t
Potential_LennardJones::_2ndDerivative (
	double rij
) const noexcept
{
	auto _rij2 = 1./(rij*rij);
	auto factor = Factor6(_rij2, nullptr);
	return m_epsilon * _rij2 * factor * ( 624.*factor - 168. );
}


double
Potential_LennardJones::Factor6 (
	double inverseRij
) const noexcept
{
	auto sigma_rij = m_sigma * inverseRij;
	auto sigma_rij2 = sigma_rij * sigma_rij;
	return sigma_rij2 * sigma_rij2 * sigma_rij2;
}


double
Potential_LennardJones::Factor6 (
	double inverseRijSq,
	std::nullptr_t
) const noexcept
{
	auto sigma_rij2 = m_sigma * m_sigma * inverseRijSq;
	return sigma_rij2 * sigma_rij2 * sigma_rij2;
}


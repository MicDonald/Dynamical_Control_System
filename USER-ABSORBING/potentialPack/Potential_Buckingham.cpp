#include "Potential_Buckingham.h"
#include <cmath>

Potential_Buckingham::Potential_Buckingham (
	double A,
	double rho,
	double C
) noexcept :
	m_A( A ),
	m_B( 1./rho ),
	m_C( C )
{}


double
Potential_Buckingham::ObjectiveFunction (
	double rij
) const noexcept
{
	auto inverseR6 = 1. / R6( rij );
	return m_A*exp(-rij*m_B) - m_C*inverseR6;
}


PotentialPair::_1stDerivative_t
Potential_Buckingham::_1stDerivative (
	double rij
) const noexcept
{
	auto inverseR7 = 1. / (R6(rij)*rij);
	return -m_A*m_B*exp(-rij*m_B) + 6.*m_C*inverseR7;
}


PotentialPair::_2ndDerivative_t
Potential_Buckingham::_2ndDerivative (
	double rij
) const noexcept
{
	auto inverseR8 = 1. / R8(rij);
	return m_A*m_B*m_B*exp(-rij*m_B) - 42.*m_C*inverseR8;
}


double
Potential_Buckingham::R6 (
	double rij
) const noexcept
{
	auto rSq = rij * rij;
	return rSq * rSq * rSq;
}


double
Potential_Buckingham::R8 (
	double rij
) const noexcept
{
	auto rSq = rij * rij;
	auto r4 = rSq * rSq;
	return r4 * r4;
}


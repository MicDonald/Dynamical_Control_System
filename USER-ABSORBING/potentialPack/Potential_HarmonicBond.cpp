#include "Potential_HarmonicBond.h"

using namespace std;


Potential_HarmonicBond::Potential_HarmonicBond (
	double k,
	double r0
) noexcept :
	m_k( k ),
	m_r0( r0 )
{}


double
Potential_HarmonicBond::ObjectiveFunction (
	double rij
) const noexcept
{
	auto dr = rij - m_r0;
	return m_k * dr * dr;
}


PotentialPair::_1stDerivative_t
Potential_HarmonicBond::_1stDerivative (
	double rij
) const noexcept
{
	return 2.*m_k * (rij - m_r0);
}


PotentialPair::_2ndDerivative_t
Potential_HarmonicBond::_2ndDerivative (
	double /*rij*/
) const noexcept
{
	return 2.*m_k;
}


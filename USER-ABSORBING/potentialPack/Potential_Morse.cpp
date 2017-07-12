#include "Potential_Morse.h"
#include <cmath>

Potential_Morse::Potential_Morse (
	double D0,
	double alpha,
	double r0
) noexcept :
	m_D0( D0 ),
	m_alpha( alpha ),
	m_r0( r0 )
{}


double
Potential_Morse::ObjectiveFunction (
	double rij
) const noexcept
{
	auto expTerm = ExpTerm(rij);
	return m_D0 * expTerm * ( expTerm - 2. );
}


PotentialPair::_1stDerivative_t
Potential_Morse::_1stDerivative (
	double rij
) const noexcept
{
	auto expTerm = ExpTerm(rij);
	return 2.*m_alpha*m_D0 * expTerm * (1.-expTerm);
}


PotentialPair::_2ndDerivative_t
Potential_Morse::_2ndDerivative (
	double rij
) const noexcept
{
	auto expTerm = ExpTerm(rij);
	return m_alpha*m_alpha*m_D0 * expTerm * (4.*expTerm-2.);
}


double
Potential_Morse::ExpTerm (
	double rij
) const noexcept
{
	return exp(m_alpha*(m_r0-rij));
}


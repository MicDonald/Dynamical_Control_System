#include "Potential_HarmonicCosineAngle.h"
#include <cmath>

using namespace Angle_Namespace;

Potential_HarmonicCosineAngle::Potential_HarmonicCosineAngle (
	double k,
	double theta0
) noexcept :
	m_k( k ),
	m_cosTheta0( cos(theta0) )
{}


double
Potential_HarmonicCosineAngle::ObjectiveFunction (
	double cosTheta,
	double /*rij_*/,
	double /*rik_*/
) const noexcept
{
	auto dCosTheta = cosTheta - m_cosTheta0;
	return m_k*dCosTheta*dCosTheta;
}


_1stDerivative_t
Potential_HarmonicCosineAngle::_1stDerivative (
	double cosTheta,
	double /*rij_*/,
	double /*rik_*/
) const noexcept
{
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = 2.*m_k*(cosTheta-m_cosTheta0);
	return _1stDeri;
}


_2ndDerivative_t
Potential_HarmonicCosineAngle::_2ndDerivative (
	double /*cosTheta*/,
	double /*rij_*/,
	double /*rik_*/
) const noexcept
{
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] = 2.*m_k;
	return _2ndDeri;
}


#include "Potential_CosineAngle.h"
#include "Math_Triangular.h"
#include <cmath>

using namespace Angle_Namespace;


Potential_CosineAngle::Potential_CosineAngle (
	double k,
	unsigned m,
	double delta
) noexcept :
	m_k( k ),
	m_m( m ),
	m_delta( delta )
{}


double
Potential_CosineAngle::ObjectiveFunction (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	return m_k*(1.+CosNTheta_delta(m_m, cosTheta, m_delta));
}


_1stDerivative_t
Potential_CosineAngle::_1stDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = m_k*_1stDCosNTheta_delta(m_m, cosTheta, m_delta);
	return _1stDeri;
}


_2ndDerivative_t
Potential_CosineAngle::_2ndDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] = m_k*_2ndDCosNTheta_delta(m_m, cosTheta, m_delta);
	return _2ndDeri;
}


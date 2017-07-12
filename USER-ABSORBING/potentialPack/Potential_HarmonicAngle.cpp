#include "Potential_HarmonicAngle.h"
#include "Math_Triangular.h"
#include <cmath>

using namespace Angle_Namespace;
constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


Potential_HarmonicAngle::Potential_HarmonicAngle (
	double k,
	double theta0
) noexcept :
	m_k( k ),
	m_theta0( theta0 )
{}


double
Potential_HarmonicAngle::ObjectiveFunction (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;
	return m_k*dTheta*dTheta;
}


PotentialAngle::_1stDerivative_t
Potential_HarmonicAngle::_1stDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta-m_theta0;
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = 2.*m_k*dTheta*_1stDTheta(cosTheta);
	return _1stDeri;
}


PotentialAngle::_2ndDerivative_t
Potential_HarmonicAngle::_2ndDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta-m_theta0;
	auto _1stDTheta = ::_1stDTheta(cosTheta);
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] = 2.*m_k*(dTheta*_2ndDTheta(cosTheta)+_1stDTheta*_1stDTheta);
	return _2ndDeri;
}


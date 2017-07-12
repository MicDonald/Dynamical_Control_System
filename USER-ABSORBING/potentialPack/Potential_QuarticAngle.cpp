#include "Potential_QuarticAngle.h"
#include "Math_Triangular.h"

using namespace Angle_Namespace;


Potential_QuarticAngle::Potential_QuarticAngle (
	double theta0,
	double k2,
	double k3,
	double k4
) noexcept :
	m_theta0( theta0 ),
	m_k2( k2 ),
	m_k3( k3 ),
	m_k4( k4 )
{}


double
Potential_QuarticAngle::ObjectiveFunction (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;
	auto dThetaSq = dTheta*dTheta;
	auto dThetaCub = dThetaSq*dTheta;
	auto dThetaQad = dThetaSq*dThetaSq;
	return m_k2*dThetaSq + m_k3*dThetaCub + m_k4*dThetaQad;
}


_1stDerivative_t
Potential_QuarticAngle::_1stDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;
	auto dThetaSq = dTheta*dTheta;
	auto dThetaCub = dThetaSq*dTheta;
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = (2.*m_k2*dTheta + 3.*m_k3*dThetaSq + 4.*m_k4*dThetaCub)*_1stDTheta(cosTheta);
	return _1stDeri;
}


_2ndDerivative_t
Potential_QuarticAngle::_2ndDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;
	auto dThetaSq = dTheta*dTheta;
	auto dThetaCub = dThetaSq*dTheta;
	auto _1stDTheta = ::_1stDTheta(cosTheta);
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] =
		(2.*m_k2*dTheta+3.*m_k3*dThetaSq+4.*m_k4*dThetaCub)*_2ndDTheta(cosTheta)
		+ (2.*m_k2+6.*m_k3*dTheta+12.*m_k4*dThetaSq)*_1stDTheta*_1stDTheta;
	return _2ndDeri;
}


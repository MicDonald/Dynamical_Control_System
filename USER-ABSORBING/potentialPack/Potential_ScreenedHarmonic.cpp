#include "Potential_ScreenedHarmonic.h"
#include "Math_Triangular.h"
#include <cmath>

using namespace Angle_Namespace;


Potential_ScreenedHarmonic::Potential_ScreenedHarmonic (
	double k,
	double theta0,
	double rho1,
	double rho2
) noexcept :
	m_k( k ),
	m_theta0( theta0 ),
	m_inverseRho1( 1./rho1 ),
	m_inverseRho2( 1./rho2 )
{}


double
Potential_ScreenedHarmonic::ObjectiveFunction (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta-m_theta0;
	auto dThetaSq = dTheta*dTheta;
	auto expFactor = -rij_*m_inverseRho1 - rik_*m_inverseRho2;
	return m_k*dThetaSq*exp(expFactor);
}


_1stDerivative_t
Potential_ScreenedHarmonic::_1stDerivative (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta-m_theta0;
	auto expFactor = -rij_*m_inverseRho1 - rik_*m_inverseRho2;
	auto expTerm = exp(expFactor);
	auto common = m_k*dTheta*expTerm;
	auto commonR = -common*dTheta;
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = 2.*common*_1stDTheta(cosTheta);
	_1stDeri[DRij] = commonR*m_inverseRho1;
	_1stDeri[DRik] = commonR*m_inverseRho2;
	return _1stDeri;
}


_2ndDerivative_t
Potential_ScreenedHarmonic::_2ndDerivative (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	auto _1stDTheta = ::_1stDTheta(cosTheta);
	auto theta = acos(cosTheta);
	auto dTheta = theta-m_theta0;
	auto expFactor = -rij_*m_inverseRho1 - rik_*m_inverseRho2;
	auto expTerm = exp(expFactor);
	auto common = m_k*expTerm;
	auto common2 = 2.*common;
	auto commonT = -common2*dTheta*_1stDTheta;
	auto commonR = common*dTheta*dTheta;
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] = common2*(_1stDTheta*_1stDTheta+dTheta*_2ndDTheta(cosTheta));
	_2ndDeri[DCosTheta_DRij] = commonT*m_inverseRho1;
	_2ndDeri[DCosTheta_DRik] = commonT*m_inverseRho2;
	_2ndDeri[DRij_DRij] = commonR*m_inverseRho1*m_inverseRho1;
	_2ndDeri[DRij_DRik] = commonR*m_inverseRho1*m_inverseRho2;
	_2ndDeri[DRik_DRik] = commonR*m_inverseRho2*m_inverseRho2;
	return _2ndDeri;
}


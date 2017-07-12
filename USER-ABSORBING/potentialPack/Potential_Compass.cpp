#include "Potential_Compass.h"
#include "Math_Triangular.h"

using namespace Angle_Namespace;


Potential_Compass::Potential_Compass (
	double A,
	double B,
	double C,
	double theta0,
	double rij0,
	double rik0
) noexcept :
	m_A( A ),
	m_B( B ),
	m_C( C ),
	m_theta0( theta0 ),
	m_rij0( rij0 ),
	m_rik0( rik0 )
{}


double
Potential_Compass::ObjectiveFunction (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	auto dTheta = acos(cosTheta) - m_theta0;
	auto dRij = rij - m_rij0;
	auto dRik = rik - m_rik0;
	return m_A*dRij*dRik + dTheta*(m_B*dRij + m_C*dRik);
}


_1stDerivative_t
Potential_Compass::_1stDerivative (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	auto dTheta = acos(cosTheta) - m_theta0;
	auto dRij = rij - m_rij0;
	auto dRik = rik - m_rik0;

	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = _1stDTheta(cosTheta)*(m_B*dRij + m_C*dRik);
	_1stDeri[DRij] = m_A*dRik + m_B*dTheta;
	_1stDeri[DRik] = m_A*dRij + m_C*dTheta;
	return _1stDeri;
}


_2ndDerivative_t
Potential_Compass::_2ndDerivative (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	auto dTheta = acos(cosTheta) - m_theta0;
	auto dRij = rij - m_rij0;
	auto dRik = rik - m_rik0;
	auto _1stDTheta = ::_1stDTheta(cosTheta);

	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] = _2ndDTheta(cosTheta)*(m_B*dRij + m_C*dRik);
	_2ndDeri[DCosTheta_DRij] = m_B*_1stDTheta;
	_2ndDeri[DCosTheta_DRik] = m_C*_1stDTheta;
	_2ndDeri[DRij_DRik] = m_A;
	return _2ndDeri;
}


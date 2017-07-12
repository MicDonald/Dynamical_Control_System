#include "Potential_StillingerWeber.h"
#include <cmath>

using namespace std;

Potential_StillingerWeber::Potential_StillingerWeber (
	double epsilon,
	double sigma,
	double a,
	double lambda,
	double gamma,
	double cosTheta0,
	double A,
	double B,
	double p,
	double q
) noexcept :
	m_epsilon( epsilon ),
	m_sigma( sigma ),
	m_a( a ),
	m_lambda( lambda ),
	m_gamma( gamma ),
	m_cosTheta0( cosTheta0 ),
	m_A( A ),
	m_B( B ),
	m_p( p ),
	m_q( q )
{}


double
Potential_StillingerWeber::ObjectiveFunction (
	double rij
) const noexcept
{
	auto sigma_r = FactorPair(rij);
	auto sigma_rP = pow(sigma_r, m_p);
	auto sigma_rQ = pow(sigma_r, m_q);
	auto expTerm = exp( FactorPair(rij, m_a*m_sigma) );
	return m_A * m_epsilon * ( m_B*sigma_rP - sigma_rQ ) * expTerm;
}


double
Potential_StillingerWeber::ObjectiveFunction (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	auto aSigma = m_a * m_sigma;
	auto expIJTerm = exp( m_gamma*FactorPair(rij, aSigma) );
	auto expIKTerm = exp( m_gamma*FactorPair(rik, aSigma) );
	auto cosTerm = cosTheta - m_cosTheta0;
	cosTerm *= cosTerm;
	return m_lambda * m_epsilon * cosTerm * expIJTerm * expIKTerm;
}


PotentialPair::_1stDerivative_t
Potential_StillingerWeber::_1stDerivative (
	double rij
) const noexcept
{
	auto sigma_r = FactorPair(rij);
	auto sigma_rP = pow(sigma_r, m_p);
	auto sigma_rQ = pow(sigma_r, m_q);
	auto aSigma = m_a * m_sigma;
	auto sigma_r__aSigma = FactorPair(rij, aSigma);
	auto expTerm = exp(sigma_r__aSigma);
	return -m_A * m_epsilon * expTerm * (
		( m_B * m_p * sigma_rP - m_q * sigma_rQ ) / rij +
		( m_B * sigma_rP - sigma_rQ ) * sigma_r__aSigma / (rij - aSigma)
	);
}


PotentialAngle::_1stDerivative_t
Potential_StillingerWeber::_1stDerivative (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	using namespace Angle_Namespace;

	auto cosTerm = cosTheta - m_cosTheta0;
	auto aSigma = m_a * m_sigma;
	auto expIJFactor = m_gamma*FactorPair(rij, aSigma);
	auto expIKFactor = m_gamma*FactorPair(rik, aSigma);
	auto expIJTerm = exp( expIJFactor );
	auto expIKTerm = exp( expIKFactor );
	auto commonFactor = m_lambda * m_epsilon * expIJTerm * expIKTerm * cosTerm;
	auto commonFactorR = -commonFactor * cosTerm / (m_gamma*m_sigma);

	PotentialAngle::_1stDerivative_t _1stDeri;
	_1stDeri[DCosTheta] = commonFactor * 2.;
	_1stDeri[DRij] = commonFactorR * expIJFactor * expIJFactor;
	_1stDeri[DRik] = commonFactorR * expIKFactor * expIKFactor;
	return _1stDeri;
}


PotentialPair::_2ndDerivative_t
Potential_StillingerWeber::_2ndDerivative (
	double rij
) const noexcept
{
	auto _r = 1./rij;
	auto _rSq = _r*_r;
	auto sigma_r = FactorPair(rij);
	auto sigma_rP = pow(sigma_r, m_p);
	auto sigma_rQ = pow(sigma_r, m_q);
	auto aSigma = m_a * m_sigma;
	auto sigma_r__aSigma = FactorPair(rij, aSigma);
	auto expTerm = exp(sigma_r__aSigma);
	auto r__aSigma = rij - aSigma;
	auto _r__aSigma = 1./r__aSigma;
	auto sigma_r__aSigmaSq = sigma_r__aSigma * _r__aSigma;
	auto sigma_r__aSigmaTri = sigma_r__aSigmaSq * _r__aSigma;
	auto sigmaTerm = sigma_r__aSigmaTri*(sigma_r__aSigma+2.);
	return m_A * m_epsilon * expTerm * (
		( m_B*m_p*(m_p+1.)*sigma_rP - m_q*(m_q+1.)*sigma_rQ ) * _rSq +
		2. * ( m_B*m_p*sigma_rP - m_q*sigma_rQ ) * _r * sigma_r__aSigmaSq +
		( m_B*sigma_rP - sigma_rQ ) * sigmaTerm
	);
}


PotentialAngle::_2ndDerivative_t
Potential_StillingerWeber::_2ndDerivative (
	double cosTheta,
	double rij,
	double rik
) const noexcept
{
	using namespace Angle_Namespace;

	double cosTerm = cosTheta - m_cosTheta0;
	double aSigma = m_a * m_sigma;
	double expIJFactor = m_gamma * FactorPair(rij, aSigma);
	double expIKFactor = m_gamma * FactorPair(rik, aSigma);
	double expIJTerm = exp(expIJFactor);
	double expIKTerm = exp(expIKFactor);
	double _rij__aSigma = 1./(rij-aSigma);
	double _rik__aSigma = 1./(rik-aSigma);
	double ijFactor = expIJFactor * _rij__aSigma;
	double ikFactor = expIKFactor * _rik__aSigma;

	double commonFactor = m_lambda * m_epsilon * expIJTerm * expIKTerm;
	double commonCosTerm = commonFactor*cosTerm;
	double commonCosSqTerm = commonCosTerm*cosTerm;

	PotentialAngle::_2ndDerivative_t _2ndDeri;
	_2ndDeri[DCosTheta_DCosTheta] = 2.*commonFactor;
	_2ndDeri[DCosTheta_DRij] = -2.*commonCosTerm*ijFactor;
	_2ndDeri[DCosTheta_DRik] = -2.*commonCosTerm*ikFactor;
	_2ndDeri[DRij_DRij] = commonCosSqTerm*ijFactor*(ijFactor+2.*_rij__aSigma);
	_2ndDeri[DRij_DRik] = commonCosSqTerm*ijFactor*ikFactor;
	_2ndDeri[DRik_DRik] = commonCosSqTerm*ikFactor*(ikFactor+2.*_rik__aSigma);
	return _2ndDeri;
}


double
Potential_StillingerWeber::FactorPair (
	double rij,
	double minus
) const noexcept
{
	return m_sigma/(rij-minus);
}


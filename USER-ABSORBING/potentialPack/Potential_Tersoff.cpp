#include "Potential_Tersoff.h"

Potential_Tersoff::Potential_Tersoff (
	double R,
	double D,
	double A,
	double B,
	double c,
	double d,
	double lambda1,
	double lambda2,
	double lambda3,
	double beta,
	double gamma,
	double m,
	double n,
	double cosTheta0
) noexcept :
	m_R( R ),
	m_D( D ),
	m_A( A ),
	m_B( B ),
	m_cSq( c*c ),
	m_dSq( d*d ),
	m_lambda1( lambda1 ),
	m_lambda2( lambda2 ),
	m_lambda3( lambda3 ),
	m_beta( beta ),
	m_gamma( gamma ),
	m_m( m ),
	m_n( n ),
	m_cosTheta0( cosTheta0 ),
	m_lowerBoundR( R-D ),
	m_upperBoundR( R+D )
{}


double
Potential_Tersoff::ObjectiveFunction (
	double rij
) const noexcept
{
	return fc(rij) * fR(rij);
}


double
Potential_Tersoff::ObjectiveFunction (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	return fc(rij_) * b(rij_,rik_,cosTheta) * fA(rij_);
}


Potential::Pair_1stDerivative_t
Potential_Tersoff::_1stDerivative (
	double rij_
) const noexcept
{
	return _1stDfc(rij_)*fR(rij_) + fc(rij_)*_1stDfR(rij_);
}


ManyBody_Namespace::_1stDerivative_t
Potential_Tersoff::_1stDerivative (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	auto fc = this->fc(rij_);
	auto _1stDfc = this->_1stDfc(rij_);
	auto fA = this->fA(rij_);
	auto _1stDfA = this->_1stDfA(rij_);
	auto fcfA = fc * fA;
	auto b = this->b(rij_, rik_, cosTheta);
	auto _1stDb = this->_1stDb(rij_, rik_, cosTheta);
	auto _1stDbIter = _1stDb.begin();
	ManyBody_Namespace::_1stDerivative_t _1stDeri;
	_1stDeri.emplace_back( b*(_1stDfc*fA + fc*_1stDfA) + fcfA*_1stDbIter->first, 0. );
	++_1stDbIter;

	for ( auto rik_Iter=rik_.begin(), cosThetaIter=cosTheta.begin();
		rik_Iter!=rik_.end();
		++rik_Iter, ++cosThetaIter, ++_1stDbIter
	)
	{
		auto _1stDb = this->_1stDb(rij_, {*rik_Iter}, {*cosThetaIter});
		_1stDeri.emplace_back (
			fcfA*_1stDbIter->first,
			fcfA*_1stDbIter->second
		);
	}
	return _1stDeri;
}


Potential::Pair_2ndDerivative_t
Potential_Tersoff::_2ndDerivative (
	double rij_
) const noexcept
{
	return _2ndDfc(rij_)*fR(rij_) +
		2.*_1stDfc(rij_)*_1stDfR(rij_) +
		fc(rij_)*_2ndDfR(rij_);
}


ManyBody_Namespace::_2ndDerivative_t
Potential_Tersoff::_2ndDerivative (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	auto fc = this->fc(rij_);
	auto _1stDfc = this->_1stDfc(rij_);
	auto _2ndDfc = this->_2ndDfc(rij_);
	auto fA = this->fA(rij_);
	auto _1stDfA = this->_1stDfA(rij_);
	auto _2ndDfA = this->_2ndDfA(rij_);
/*	return _2ndDfc(rij_) * densityI.density * fA +
		2. * _1stDfc * (densityI._1stDdensity * fA +
				densityI.density * _1stDfA) +
		fc * ( densityI._2ndDdensity * fA +
			2. * densityI._1stDdensity * _1stDfA +
			densityI.density * _2ndDfA(rij_) );*/
}


double
Potential_Tersoff::fc (
	double rij_
) const noexcept
{
	constexpr double halfPI = 0.5 * PI;
	if ( rij_ < m_lowerBoundR )
		return 1.;
	else if ( rij_ > m_upperBoundR )
		return 0.;
	return 0.5*(1. - sin(halfPI * (rij_-m_R)/m_D));
}


double
Potential_Tersoff::_1stDfc (
	double rij_
) const noexcept
{
	if ( rij_ < m_lowerBoundR || rij_ > m_upperBoundR )
		return 0.;
	constexpr double halfPI = 0.5 * PI;
	constexpr double factor = -0.25*PI;
	auto inverseD = 1./m_D;
	return factor*inverseD * cos(halfPI*(rij_-m_R)*inverseD);
}


double
Potential_Tersoff::_2ndDfc (
	double rij_
) const noexcept
{
	if ( rij_ < m_lowerBoundR || rij_ > m_upperBoundR )
		return 0.;
	constexpr double halfPI = 0.5*PI;
	constexpr double factor = 0.125*PI*PI;
	double inverseD = 1./m_D;
	double inverseDSq = inverseD * inverseD;
	return factor*inverseDSq * sin(halfPI*(rij_-m_R)*inverseD);
}


double
Potential_Tersoff::fR (
	double rij_
) const noexcept
{
	return m_A * exp(-m_lambda1 * rij_);
}


double
Potential_Tersoff::_1stDfR (
	double rij_
) const noexcept
{
	return -m_lambda1 * m_A * exp(-m_lambda1 * rij_);
}


double
Potential_Tersoff::_2ndDfR (
	double rij_
) const noexcept
{
	return m_lambda1*m_lambda1 * m_A * exp(-m_lambda1 * rij_);
}


double
Potential_Tersoff::fA (
	double rij_
) const noexcept
{
	return -m_B * exp(-m_lambda2 * rij_);
}


double
Potential_Tersoff::_1stDfA (
	double rij_
) const noexcept
{
	return m_lambda2 * m_B * exp(-m_lambda2 * rij_);
}


double
Potential_Tersoff::_2ndDfA (
	double rij_
) const noexcept
{
	return -m_lambda2*m_lambda2 * m_B * exp(-m_lambda2 * rij_);
}


double
Potential_Tersoff::Zeta (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	return fc(rik_) * g(cosTheta) * exp(pow(m_lambda3*(rij_-rik_), m_m));
}


Angle_Namespace::_1stDerivative_t
Potential_Tersoff::_1stDZeta (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	using namespace Angle_Namespace;

	auto fc = this->fc(rik_);
	auto g = this->g(cosTheta);
	auto rijMinusRik = rij_ - rik_;
	auto expFactor = pow(m_lambda3 * rijMinusRik, m_m);
	auto factor = m_m * expFactor / rijMinusRik;
	auto expTerm = exp(expFactor);

	_1stDerivative_t _1stDeri;
	_1stDeri[DCosTheta] = fc * _1stDg(cosTheta) * expTerm;
	_1stDeri[DRij] = fc * g * expTerm * factor;
	_1stDeri[DRik] = g * expTerm * (_1stDfc(rik_) - fc*factor);
	return _1stDeri;
}


Angle_Namespace::_2ndDerivative_t
Potential_Tersoff::_2ndDZeta (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	using namespace Angle_Namespace;

	auto fc = this->fc(rik_);
	auto _1stDfc = this->_1stDfc(rik_);
	auto g = this->g(cosTheta);
	auto _1stDg = this->_1stDg(cosTheta);
	auto rijMinusRik = rij_ - rik_;
	auto _rijMinusRik = 1./rijMinusRik;
	auto expFactor = pow(m_lambda3 * rijMinusRik, m_m);
	auto factor = m_m * expFactor * _rijMinusRik;
	auto factor2 = (m_m-1.) * factor * _rijMinusRik;
	auto factor3 = factor*factor + factor2;
	auto expTerm = exp(expFactor);

	_2ndDerivative_t _2ndDeri;
	_2ndDeri[DCosTheta_DCosTheta] = fc * _2ndDg(cosTheta) * expTerm;
	_2ndDeri[DCosTheta_DRij] = fc * _1stDg * expTerm * factor;
	_2ndDeri[DCosTheta_DRik] = _1stDg * expTerm * (_1stDfc - fc*factor);
	_2ndDeri[DRij_DRij] = fc * g * expTerm * factor3;
	_2ndDeri[DRij_DRik] = g*expTerm * (_1stDfc*factor - fc*factor3);
	_2ndDeri[DRik_DRik] = g*expTerm * (_2ndDfc(rik_) - 2.*_1stDfc*factor + fc*factor3);
	return _2ndDeri;
}


double
Potential_Tersoff::g (
	double cosTheta
) const noexcept
{
	auto dcosTheta = cosTheta - m_cosTheta0;
	return m_gamma * ( 1. + m_cSq/m_dSq - m_cSq/(m_dSq+dcosTheta*dcosTheta) );
}


double
Potential_Tersoff::_1stDg (
	double cosTheta
) const noexcept
{
	auto dcosTheta = cosTheta - m_cosTheta0;
	auto denominator = m_dSq + dcosTheta*dcosTheta;
	auto inverseDenominatorSq = 1./(denominator*denominator);
	return m_gamma*2.*m_cSq*dcosTheta * inverseDenominatorSq;
}


double
Potential_Tersoff::_2ndDg (
	double cosTheta
) const noexcept
{
	auto dcosTheta = cosTheta - m_cosTheta0;
	auto dcosThetaSq = dcosTheta * dcosTheta;
	auto denominator = m_dSq + dcosThetaSq;
	auto inverseDenominator = 1./denominator;
	auto inverseDenominatorSq = inverseDenominator*inverseDenominator;
	return m_gamma*m_cSq*(2.-8.*inverseDenominator*dcosThetaSq)*inverseDenominatorSq;
}


double
Potential_Tersoff::b (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	auto rik_Iter = rik_.begin();
	auto cosThetaIter = cosTheta.begin();
	auto zetaI = 0.;
	for (; rik_Iter!=rik_.end(); ++rik_Iter, ++cosThetaIter)
		zetaI += Zeta(*cosThetaIter, rij_, *rik_Iter);
	return pow( 1.+pow(m_beta*zetaI, m_n), -0.5/m_n );
}


ManyBody_Namespace::_1stDerivative_t
Potential_Tersoff::_1stDb (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	ManyBody_Namespace::_1stDerivative_t _1stDeri;
	_1stDeri.emplace_back(0., 0.);
	auto zetaIJ = 0.;
	auto dZetaIJ = 0.;
	for ( auto rik_Iter=rik_.begin(), cosThetaIter=cosTheta.begin();
		rik_Iter!=rik_.end();
		++rik_Iter, ++cosThetaIter
	)
	{
		auto zeta = this->Zeta(*cosThetaIter, rij_, *rik_Iter);
		auto _1stDZeta = this->_1stDZeta(*cosThetaIter, rij_, *rik_Iter);
		auto betaZetaEN = pow(m_beta*zeta, m_n);
		auto factor = -0.5*pow(1.+betaZetaEN, -0.5/m_n-1.) * betaZetaEN/zeta;
		_1stDeri.emplace_back (
			factor*_1stDZeta[Angle_Namespace::DRik],
			factor*_1stDZeta[Angle_Namespace::DCosTheta]
		);
		zetaIJ += zeta;
		dZetaIJ += _1stDZeta[Angle_Namespace::DRij];
	}
	auto betaZetaEN = pow(m_beta*zetaIJ, m_n);
	auto factor = -0.5*pow(1.+betaZetaEN, -0.5/m_n-1.) * betaZetaEN/zetaIJ;
	_1stDeri.front().first = factor*dZetaIJ;

	return _1stDeri;
}


ManyBody_Namespace::_2ndDerivative_t
Potential_Tersoff::_2ndDb (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
}


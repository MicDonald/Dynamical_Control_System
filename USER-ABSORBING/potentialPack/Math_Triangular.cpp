#include "Math_Triangular.h"
#include <cmath>

//to do:
// need to check inverseSinTheta

double
Math_Triangular::_1stDTheta (
	double cosTheta
) const noexcept
{
	return -1./sqrt(1.-cosTheta*cosTheta);
}


double
Math_Triangular::_2ndDTheta (
	double cosTheta
) const noexcept
{
	auto sinThetaSq = 1.-cosTheta*cosTheta;
	auto sinTheta = sqrt(sinThetaSq);
	return -cosTheta/(sinThetaSq*sinTheta);
}


double
Math_Triangular::CosNTheta (
	unsigned n,
	double cosTheta
) const noexcept
{
	double cosSqTheta;
	switch( n )
	{
	case 0:
		return 1.;
	case 1:
		return cosTheta;
	case 2:
		return 2.*cosTheta*cosTheta-1.;
	case 3:
		return (4.*cosTheta*cosTheta-3.)*cosTheta;
	case 4:
		cosSqTheta = cosTheta * cosTheta;
		return 8.*cosSqTheta*(cosSqTheta-1) + 1.;
	case 5:
		cosSqTheta = cosTheta * cosTheta;
		return cosTheta*(cosSqTheta*(16.*cosSqTheta-20.)+5.);
	case 6:
		cosSqTheta = cosTheta * cosTheta;
		return cosSqTheta*(cosSqTheta*(32.*cosSqTheta-48.)+18.)-1.;
	}
	return cos( n*acos(cosTheta) );
}


double
Math_Triangular::_1stDCosNTheta (
	unsigned n,
	double cosTheta
) const noexcept
{
	double cosThetaSq;
	switch( n )
	{
	case 0:
		return 0.;
	case 1:
		return 1.;
	case 2:
		return 4.*cosTheta;
	case 3:
		return 12.*cosTheta*cosTheta - 3.;
	case 4:
		return cosTheta*(32.*cosTheta*cosTheta-16.);
	case 5:
		cosThetaSq = cosTheta * cosTheta;
		return 80.*cosThetaSq*cosThetaSq - 60.*cosThetaSq + 5.;
	case 6:
		cosThetaSq = cosTheta * cosTheta;
		return cosTheta*(cosThetaSq*(192.*cosThetaSq-192.)+36.);
	}
	auto phi = acos(cosTheta);
	return n*sin(n*phi)/sin(phi);
}


double
Math_Triangular::_2ndDCosNTheta (
	unsigned n,
	double cosTheta
) const noexcept
{
	double cosThetaSq;
	switch( n )
	{
	case 0:
	case 1:
		return 0.;
	case 2:
		return 4.;
	case 3:
		return 24.*cosTheta;
	case 4:
		return 96.*cosTheta*cosTheta-16.;
	case 5:
		cosThetaSq = cosTheta * cosTheta;
		return cosTheta*(320.*cosThetaSq - 120.);
	case 6:
		cosThetaSq = cosTheta * cosTheta;
		return cosThetaSq*(cosThetaSq*960.-576.)+36.;
	}
	auto phi = acos(cosTheta);
	auto sinTheta = sin(phi);
	auto sinNTheta = sin(n*phi);
	auto sinThetaSq = sinTheta*sinTheta;
	return -n*n*CosNTheta(n, cosTheta)/sinThetaSq + n*sinNTheta*cosTheta/(sinTheta*sinThetaSq);
}


double
Math_Triangular::CosNTheta_delta (
	unsigned n,
	double cosTheta,
	double delta
) const noexcept
{
	auto triPair = TriNTheta(n, cosTheta);
	return triPair.first*cos(delta) + triPair.second*sin(delta);
}


double
Math_Triangular::_1stDCosNTheta_delta (
	unsigned n,
	double cosTheta,
	double delta
) const noexcept
{
	auto _1stDTriPair = _1stDTriNTheta(n, cosTheta);
	return _1stDTriPair.first*cos(delta) + _1stDTriPair.second*sin(delta);
}


double
Math_Triangular::_2ndDCosNTheta_delta (
	unsigned n,
	double cosTheta,
	double delta
) const noexcept
{
	auto _2ndDTriPair = _2ndDTriNTheta(n, cosTheta);
	return _2ndDTriPair.first*cos(delta) + _2ndDTriPair.second*sin(delta);
}


std::pair<double, double>
Math_Triangular::TriNTheta (
	unsigned n,
	double cosTheta
) const noexcept
{
	auto cosNTheta = CosNTheta(n, cosTheta);

	double cosThetaSq = cosTheta*cosTheta;
	double sinTheta = sqrt(1.-cosThetaSq);
	switch( n )
	{
	case 0:
		return {cosNTheta, 0.};
	case 1:
		return {cosNTheta, sinTheta};
	case 2:
		return {cosNTheta, 2.*cosTheta*sinTheta};
	case 3:
		return {cosNTheta, sinTheta*(4.*cosThetaSq-1.)};
	case 4:
		return {cosNTheta, sinTheta*cosTheta*(8.*cosThetaSq-4.)};
	case 5:
		return {cosNTheta, sinTheta*((16.*cosThetaSq-12.)*cosTheta+1.)};
	case 6:
		return {cosNTheta, sinTheta*cosTheta*((32.*cosThetaSq-32.)*cosTheta-16.)};
	}
	auto theta = acos(cosTheta);
	return {cosNTheta, sin(n*theta)};
}


std::pair<double, double>
Math_Triangular::_1stDTriNTheta (
	unsigned n,
	double cosTheta
) const noexcept
{
	auto dCosNTheta = _1stDCosNTheta(n, cosTheta);

	auto cosNTheta = CosNTheta(n, cosTheta);
	auto cosThetaSq = cosTheta*cosTheta;
	auto sinTheta = sqrt(1.-cosThetaSq);
	auto inverseSinTheta = sinTheta<1.e-10 ? 1.e-10 : 1./sinTheta;
	return {dCosNTheta, -cosNTheta*inverseSinTheta*n};
}


std::pair<double, double>
Math_Triangular::_2ndDTriNTheta (
	unsigned n,
	double cosTheta
) const noexcept
{
	auto dSqCosTheta = _2ndDCosNTheta(n, cosTheta);

	auto dCosNTheta = _1stDCosNTheta(n, cosTheta);
	auto cosNTheta= CosNTheta(n, cosTheta);
	auto sinThetaSq = 1.-cosTheta*cosTheta;
	auto sinTheta = sqrt(sinThetaSq);
	auto inverseSinTheta = sinTheta<1.e-10 ? 1.e-10 : 1./sinTheta;
	return {dSqCosTheta, -(dCosNTheta + cosTheta*cosNTheta/sinThetaSq)*inverseSinTheta*n};
}


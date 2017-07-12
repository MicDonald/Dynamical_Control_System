#include "Potential_HarmonicDihedral.h"
#include "Math_Triangular.h"
#include <cmath>

using namespace Dihedral_Namespace;


Potential_HarmonicDihedral::Potential_HarmonicDihedral (
	double k,
	double phi0
) noexcept :
	m_k( k ),
	m_phi0( phi0 )
{}


double
Potential_HarmonicDihedral::ObjectiveFunction (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	auto dPhi = acos(cosPhi) - m_phi0;
	return m_k*dPhi*dPhi;
}


PotentialDihedral::_1stDerivative_t
Potential_HarmonicDihedral::_1stDerivative (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	auto dPhi = acos(cosPhi) - m_phi0;

	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosPhi] = 2.*m_k*dPhi*_1stDTheta(cosPhi);
	return _1stDeri;
}


PotentialDihedral::_2ndDerivative_t
Potential_HarmonicDihedral::_2ndDerivative (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	auto dPhi = acos(cosPhi) - m_phi0;
	auto _1stDPhi = _1stDTheta(cosPhi);

	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosPhi_DCosPhi] = 2.*m_k*(dPhi*_2ndDTheta(cosPhi) + _1stDPhi*_1stDPhi);
	return _2ndDeri;
}


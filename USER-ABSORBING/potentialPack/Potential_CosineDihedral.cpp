#include "Potential_CosineDihedral.h"
#include "Math_Triangular.h"
#include <cmath>

using namespace Dihedral_Namespace;


Potential_CosineDihedral::Potential_CosineDihedral (
	double k,
	int d,
	unsigned n
) noexcept :
	m_k( k ),
	m_d( d ),
	m_n( n )
{}


double
Potential_CosineDihedral::ObjectiveFunction (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	return m_k*( 1. + m_d*CosNTheta(m_n, cosPhi) );
}


PotentialDihedral::_1stDerivative_t
Potential_CosineDihedral::_1stDerivative (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosPhi] = m_k*m_d*_1stDCosNTheta(m_n, cosPhi);
	return _1stDeri;
}


PotentialDihedral::_2ndDerivative_t
Potential_CosineDihedral::_2ndDerivative (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosPhi_DCosPhi] = m_k*m_d*_2ndDCosNTheta(m_n, cosPhi);
	return _2ndDeri;
}


#include "Potential_HarmonicCosineDihedral.h"
#include "Math_Triangular.h"

using namespace Dihedral_Namespace;


Potential_HarmonicCosineDihedral::Potential_HarmonicCosineDihedral (
	double k,
	double cosPhi0
) noexcept :
	m_k( k ),
	m_cosPhi0( cosPhi0 )
{}


double
Potential_HarmonicCosineDihedral::ObjectiveFunction (
	double cosPhi,
	double /*rij*/,
	double /*rjk*/,
	double /*rkl*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	auto dCosPhi = cosPhi - m_cosPhi0;
	return m_k*dCosPhi*dCosPhi;
}


_1stDerivative_t
Potential_HarmonicCosineDihedral::_1stDerivative (
	double cosPhi,
	double /*rij*/,
	double /*rjk*/,
	double /*rkl*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	auto dCosPhi = cosPhi - m_cosPhi0;

	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosPhi] = 2.*m_k*dCosPhi;
	return _1stDeri;
}


_2ndDerivative_t
Potential_HarmonicCosineDihedral::_2ndDerivative (
	double /*cosPhi*/,
	double /*rij*/,
	double /*rjk*/,
	double /*rkl*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosPhi_DCosPhi] = 2.*m_k;
	return _2ndDeri;
}


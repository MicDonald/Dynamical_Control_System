#ifndef POTENTIAL_HARMONICDIHEDRAL_H_INCLUDED
#define POTENTIAL_HARMONICDIHEDRAL_H_INCLUDED

#include "PotentialDihedral.h"

class Potential_HarmonicDihedral : public PotentialDihedral {
public:
	Potential_HarmonicDihedral (
		double k = 0.5,
		double phi0 = 0.3
	) noexcept;

private:
	double
	ObjectiveFunction (
		double cosPhi,
		double rij_,
		double rjk_,
		double rkl_,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;


	_1stDerivative_t
	_1stDerivative (
		double cosPhi,
		double rij_,
		double rjk_,
		double rkl_,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;


	_2ndDerivative_t
	_2ndDerivative (
		double cosPhi,
		double rij_,
		double rjk_,
		double rkl_,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;

private:
	double m_k;
	double m_phi0;
};

#endif // POTENTIAL_HARMONICDIHEDRAL_H_INCLUDED

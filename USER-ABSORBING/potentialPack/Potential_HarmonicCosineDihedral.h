#ifndef POTENTIAL_HARMONICCOSINEDIHEDRAL_H_INCLUDED
#define POTENTIAL_HARMONICCOSINEDIHEDRAL_H_INCLUDED

#include "PotentialDihedral.h"


class Potential_HarmonicCosineDihedral : public PotentialDihedral {
public:
	Potential_HarmonicCosineDihedral (
		double k = 0.3,
		double cosPhi0 = 0.5
	) noexcept;

private:
	double
	ObjectiveFunction (
		double cosPhi,
		double rij,
		double rjk,
		double rkl,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;


	_1stDerivative_t
	_1stDerivative (
		double cosPhi,
		double rij,
		double rjk,
		double rkl,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;


	_2ndDerivative_t
	_2ndDerivative (
		double cosPhi,
		double rij,
		double rjk,
		double rkl,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;

private:
	double m_k;
	double m_cosPhi0;
};

#endif // POTENTIAL_HARMONICCOSINEDIHEDRAL_H_INCLUDED

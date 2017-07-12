#ifndef POTENTIAL_COSINEDIHEDRAL_H_INCLUDED
#define POTENTIAL_COSINEDIHEDRAL_H_INCLUDED

#include "PotentialDihedral.h"

class Potential_CosineDihedral : public PotentialDihedral {
public:
	Potential_CosineDihedral (
		double k = 0.5,
		int d = 1,
		unsigned n = 3
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
	int m_d;
	unsigned m_n;
};

#endif // POTENTIAL_COSINEDIHEDRAL_H_INCLUDED

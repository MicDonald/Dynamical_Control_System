#ifndef POTENTIAL_BUCKINGHAM_H_INCLUDED
#define POTENTIAL_BUCKINGHAM_H_INCLUDED

#include "PotentialPair.h"

class Potential_Buckingham : public PotentialPair {
public:
	Potential_Buckingham (
		double A = 0.13169812355066,
		double rho = 2.993,
		double C = 12.
	) noexcept;

private:
	double
	ObjectiveFunction (
		double rij
	) const noexcept override;


	_1stDerivative_t
	_1stDerivative (
		double rij
	) const noexcept override;


	_2ndDerivative_t
	_2ndDerivative (
		double rij
	) const noexcept override;

protected:
	double
	R6 (
		double rij
	) const noexcept;


	double
	R8 (
		double rij
	) const noexcept;

private:
	double m_A;
	double m_B;
	double m_C;
};

#endif // POTENTIAL_BUCKINGHAM_H_INCLUDED

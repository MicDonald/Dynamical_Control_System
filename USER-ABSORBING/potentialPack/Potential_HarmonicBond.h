#ifndef POTENTIAL_HARMONICBOND_H_INCLUDED
#define POTENTIAL_HARMONICBOND_H_INCLUDED

#include "PotentialPair.h"

class Potential_HarmonicBond : public PotentialPair {
public:
	Potential_HarmonicBond (
		double k = 0.5,
		double r0 = 1.
	) noexcept;

protected:
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

private:
	double m_k;
	double m_r0;
};

#endif // POTENTIAL_HARMONICBOND_H_INCLUDED

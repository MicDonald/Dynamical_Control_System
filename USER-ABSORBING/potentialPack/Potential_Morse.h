#ifndef POTENTIAL_MORSE_H_INCLUDED
#define POTENTIAL_MORSE_H_INCLUDED

#include "PotentialPair.h"

class Potential_Morse : public PotentialPair {
public:
	Potential_Morse (
		double D0 = 1.,
		double alpha = 1.,
		double r0 = 1.
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
	double ExpTerm (
		double rij
	) const noexcept;

private:
	double m_D0;
	double m_alpha;
	double m_r0;
};

#endif // POTENTIAL_MORSE_H_INCLUDED

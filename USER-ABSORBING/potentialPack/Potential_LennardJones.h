#ifndef POTENTIAL_LENNARDJONES_H_INCLUDED
#define POTENTIAL_LENNARDJONES_H_INCLUDED

#include "PotentialPair.h"

class Potential_LennardJones : public PotentialPair {
public:
	Potential_LennardJones (
		double epsilon = 1.,
		double sigma = 1.
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
	Factor6 (
		double inverseRij
	) const noexcept;


	double
	Factor6 (
		double inverseRijSq,
		std::nullptr_t
	) const noexcept;

private:
	double m_epsilon;
	double m_sigma;
};

#endif // POTENTIAL_LENNARDJONES_H_INCLUDED

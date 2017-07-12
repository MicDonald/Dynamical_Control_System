#ifndef POTENTIAL_COSINEANGLE_H_INCLUDED
#define POTENTIAL_COSINEANGLE_H_INCLUDED

#include "PotentialAngle.h"

class Potential_CosineAngle : public PotentialAngle {
public:
	Potential_CosineAngle (
		double k = 1.,
		unsigned m = 1,
		double delta = 0.33
	) noexcept;

private:
	double
	ObjectiveFunction (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;


	_1stDerivative_t
	_1stDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;


	_2ndDerivative_t
	_2ndDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;

private:
	double m_k;
	unsigned m_m;
	double m_delta;
};

#endif // POTENTIAL_COSINEANGLE_H_INCLUDED

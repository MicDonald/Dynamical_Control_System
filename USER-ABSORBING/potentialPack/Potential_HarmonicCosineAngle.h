#ifndef POTENTIAL_HARMONICCOSINEANGLE_H_INCLUDED
#define POTENTIAL_HARMONICCOSINEANGLE_H_INCLUDED

#include "PotentialAngle.h"

class Potential_HarmonicCosineAngle : public PotentialAngle {
public:
	Potential_HarmonicCosineAngle (
		double k = 1.,
		double theta0 = 1.04
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
	double m_cosTheta0;
};

#endif // POTENTIAL_HARMONICCOSINEANGLE_H_INCLUDED

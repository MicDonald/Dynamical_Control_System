#ifndef POTENTIAL_HARMONICANGLE_H_INCLUDED
#define POTENTIAL_HARMONICANGLE_H_INCLUDED

#include "PotentialAngle.h"

class Potential_HarmonicAngle : public PotentialAngle {
public:
	Potential_HarmonicAngle (
		double k = 0.5,
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
	double m_theta0;
};

#endif // POTENTIAL_HARMONICANGLE_H_INCLUDED

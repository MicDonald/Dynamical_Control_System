#ifndef POTENTIAL_QUARTICANGLE_H_INCLUDED
#define POTENTIAL_QUARTICANGLE_H_INCLUDED

#include "PotentialAngle.h"


class Potential_QuarticAngle : public PotentialAngle {
public:
	Potential_QuarticAngle (
		double theta0 = 0.65,
		double k2 = 2.,
		double k3 = 1.1,
		double k4 = 0.5
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
	double m_theta0;
	double m_k2;
	double m_k3;
	double m_k4;
};

#endif // POTENTIAL_QUARTICANGLE_H_INCLUDED

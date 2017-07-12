#ifndef POTENTIAL_COMPASS_H_INCLUDED
#define POTENTIAL_COMPASS_H_INCLUDED

#include "PotentialAngle.h"

class Potential_Compass : public PotentialAngle {
public:
	Potential_Compass (
		double A = 0.2,
		double B = 0.3,
		double C = 0.25,
		double theta0 = 0.33,
		double rij0 = 1.1,
		double rik0 = 1.5
	) noexcept;


	double
	ObjectiveFunction (
		double cosTheta,
		double rij,
		double rik
	) const noexcept;


	_1stDerivative_t
	_1stDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept;


	_2ndDerivative_t
	_2ndDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept;

private:
	double m_A;
	double m_B;
	double m_C;
	double m_theta0;
	double m_rij0;
	double m_rik0;
};

#endif // POTENTIAL_COMPASS_H_INCLUDED

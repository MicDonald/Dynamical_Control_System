#ifndef POTENTIAL_SCREENEDHARMONIC_H_INCLUDED
#define POTENTIAL_SCREENEDHARMONIC_H_INCLUDED

#include "PotentialAngle.h"

class Potential_ScreenedHarmonic : public PotentialAngle {
public:
	Potential_ScreenedHarmonic (
		double k = 0.5,
		double theta0 = 0.3,
		double rho1 = 0.17,
		double rho2 = 0.13
	) noexcept;

private:
	double
	ObjectiveFunction (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;


	_1stDerivative_t
	_1stDerivative (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;


	_2ndDerivative_t
	_2ndDerivative (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;

private:
	double m_k;
	double m_theta0;
	double m_inverseRho1;
	double m_inverseRho2;
};

#endif // POTENTIAL_SCREENEDHARMONIC_H_INCLUDED

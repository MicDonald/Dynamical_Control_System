#ifndef POTENTIAL_TERSOFF_H_INCLUDED
#define POTENTIAL_TERSOFF_H_INCLUDED

#include "Potential.h"
#include <cmath>
#include <memory>


class Potential_Tersoff : public Potential {
public:
	Potential_Tersoff (
		double R = 1.95,
		double D = 0.15,
		double A = 1393.6,
		double B = 346.74,
		double c = 38049.,
		double d = 4.3484,
		double lambda1 = 3.4879,
		double lambda2 = 2.2119,
		double lambda3 = 0.,
		double beta = 1.5724e-7,
		double gamma = 1.,
		double m = 3.,
		double n = 0.72751,
		double cosTheta0 = -0.57058
	) noexcept;

private:
	double
	ObjectiveFunction (
		double rij_
	) const noexcept override;


	double
	ObjectiveFunction (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept override;


	Potential::Pair_1stDerivative_t
	_1stDerivative (
		double rij_
	) const noexcept override;


	Potential::ManyBody_1stDerivative_t
	_1stDerivative (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept override;


	Potential::Pair_2ndDerivative_t
	_2ndDerivative (
		double rij_
	) const noexcept override;


	Potential::ManyBody_2ndDerivative_t
	_2ndDerivative (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept override;

public:
	double
	fc (
		double rij_
	) const noexcept;


	double
	_1stDfc (
		double rij_
	) const noexcept;


	double
	_2ndDfc (
		double rij_
	) const noexcept;

	//-------------------------------------//

	double
	fR (
		double rij_
	) const noexcept;


	double
	_1stDfR (
		double rij_
	) const noexcept;


	double
	_2ndDfR (
		double rij_
	) const noexcept;

	//-------------------------------------//

	double
	fA (
		double rij_
	) const noexcept;


	double
	_1stDfA (
		double rij_
	) const noexcept;


	double
	_2ndDfA (
		double rij_
	) const noexcept;

	//-------------------------------------//

	double
	Zeta (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;


	Angle_Namespace::_1stDerivative_t
	_1stDZeta (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;


	Angle_Namespace::_2ndDerivative_t
	_2ndDZeta (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;

	//-------------------------------------//

	double
	g (
		double cosTheta
	) const noexcept;


	double
	_1stDg (
		double cosTheta
	) const noexcept;


	double
	_2ndDg (
		double cosTheta
	) const noexcept;

	//-------------------------------------//

	double
	b (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept;


	ManyBody_Namespace::_1stDerivative_t
	_1stDb (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept;


	ManyBody_Namespace::_2ndDerivative_t
	_2ndDb (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept;

private:
	double m_R;
	double m_D;
	double m_A;
	double m_B;
	double m_cSq;
	double m_dSq;
	double m_lambda1;
	double m_lambda2;
	double m_lambda3;
	double m_beta;
	double m_gamma;
	double m_m;
	double m_n;
	double m_cosTheta0;

	double m_lowerBoundR;
	double m_upperBoundR;
};

#endif // POTENTIAL_TERSOFF_H_INCLUDED

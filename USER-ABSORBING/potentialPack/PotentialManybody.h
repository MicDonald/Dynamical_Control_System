#ifndef POTENTIALMANYBODY_H_INCLUDED
#define POTENTIALMANYBODY_H_INCLUDED

#include "Potential.h"


class PotentialManybody : public Potential {
public:
	using _1stDerivative_t = Manybody_Namespace::_1stDerivative_t;
	using _2ndDerivative_t = Manybody_Namespace::_2ndDerivative_t;
	using _1stD_Func = std::function< _1stDerivative_t(double, const std::vector<double>&, const std::vector<double>&) >;
	using _2ndD_Func = std::function< _2ndDerivative_t(double, const std::vector<double>&, const std::vector<double>&) >;


	double
	Energy (
		const array3d& rij,
		const std::vector<array3d>& rik
	) const noexcept override;

	//-------------------------------------//

	std::vector<array3d> // fij, fik[]
	Force (
		const array3d& rij,
		const std::vector<array3d>& rik
	) const noexcept override;


	std::vector<array3d>
	Force (
		const array3d& rij,
		const std::vector<array3d>& rik,
		FiniteDifference_t
	) const noexcept override;

	//-------------------------------------//

	std::vector<matrix3d> // Hij, Hik[], Hjk[], Hk[]k[]
	Hessian (
		const array3d& rij,
		const std::vector<array3d>& rik
	) const noexcept override;


	std::vector<matrix3d>
	Hessian (
		const array3d& rij,
		const std::vector<array3d>& rik,
		FiniteDifference_t
	) const noexcept override;

private:
	std::vector<array3d>
	ForceImp (
		const array3d& rij,
		const std::vector<array3d>& rik,
		_1stD_Func _1stDeriFunc
	) const noexcept;


	std::vector<matrix3d>
	HessianImp (
		const array3d& rij,
		const std::vector<array3d>& rik,
		_1stD_Func _1stDeriFunc,
		_2ndD_Func _2ndDeriFunc
	) const noexcept;

	//-------------------------------------//
	//-------------------------------------//
	//-------------------------------------//

	virtual double
	ObjectiveFunction (
		double /*rij_*/,
		const std::vector<double>& /*rik_*/,
		const std::vector<double>& /*cosTheta*/
	) const noexcept = 0;


	virtual _1stDerivative_t // {drij, dcosTheta}, {drik, dcosTheta}[]
	_1stDerivative (
		double rij_,
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept;


	virtual _2ndDerivative_t // drij_drij, drij_drik[], drij_dcosTheta[]
	_2ndDerivative (		// drik[]_drik[], drik[]_dcosTheta[]
		double rij_,		// dcosTheta[]_dcosTheta[]
		const std::vector<double>& rik_,
		const std::vector<double>& cosTheta
	) const noexcept;
};

#endif // POTENTIALMANYBODY_H_INCLUDED

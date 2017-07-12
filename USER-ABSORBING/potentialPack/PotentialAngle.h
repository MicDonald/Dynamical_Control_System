#ifndef POTENTIALANGLE_H_INCLUDED
#define POTENTIALANGLE_H_INCLUDED

#include "Potential.h"


class PotentialAngle : public Potential {
public:
	using _1stDerivative_t = Angle_Namespace::_1stDerivative_t;
	using _2ndDerivative_t = Angle_Namespace::_2ndDerivative_t;
	using _1stD_Func = std::function< _1stDerivative_t(double, double, double) >;
	using _2ndD_Func = std::function< _2ndDerivative_t(double, double, double) >;


	double
	Energy (
		const array3d& rij,
		const array3d& rik
	) const noexcept override;

	//-------------------------------------//

	std::array<array3d, 2> // fij, fik
	Force (
		const array3d& rij,
		const array3d& rik
	) const noexcept override;


	std::array<array3d, 2>
	Force (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept override;

	//-------------------------------------//

	std::array<matrix3d, 3> // Hij, Hik, Hjk
	Hessian (
		const array3d& rij,
		const array3d& rik
	) const noexcept override;


	std::array<matrix3d, 3>
	Hessian (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept override;


	void
	_1stDebug(
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	void
	_2ndDebug(
		const array3d& rij,
		const array3d& rik
	) const noexcept;

private:
	std::array<array3d, 2>
	ForceImp (
		const array3d& rij,
		const array3d& rik,
		_1stD_Func _1stDeriFunc
	) const noexcept;


	std::array<matrix3d, 3>
	HessianImp (
		const array3d& rij,
		const array3d& rik,
		_1stD_Func _1stDeriFunc,
		_2ndD_Func _2ndDeriFunc
	) const noexcept;

	//-------------------------------------//
	//-------------------------------------//
	//-------------------------------------//

	virtual double
	ObjectiveFunction (
		double /*cosTheta*/,
		double /*rij_*/,
		double /*rik_*/
	) const noexcept = 0;


	virtual _1stDerivative_t // dcosTheta, drij, drik
	_1stDerivative (
		double cosTheta,
		double rij_,
		double rik_
	) const noexcept;


	virtual _2ndDerivative_t	// dcosTheta_dcosTheta, dcosTheta_drij, dcosTheta_drik
	_2ndDerivative (		// drij_drij, drij_drik
		double cosTheta,	// drik_drik
		double rij_,
		double rik_
	) const noexcept;
};

#endif // POTENTIALANGLE_H_INCLUDED

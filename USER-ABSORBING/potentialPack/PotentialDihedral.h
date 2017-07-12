#ifndef POTENTIALDIHEDRAL_H_INCLUDED
#define POTENTIALDIHEDRAL_H_INCLUDED

#include "Potential.h"


class PotentialDihedral : public Potential {
public:
	using _1stDerivative_t = Dihedral_Namespace::_1stDerivative_t;
	using _2ndDerivative_t = Dihedral_Namespace::_2ndDerivative_t;
	using _1stD_Func = std::function< _1stDerivative_t(double, double, double, double, double, double) >;
	using _2ndD_Func = std::function< _2ndDerivative_t(double, double, double, double, double, double) >;


	double
	Energy (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept override;

	//-------------------------------------//

	std::array<array3d, 3> // fji, fjk, fkl
	Force (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept override;


	std::array<array3d, 3>
	Force (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		FiniteDifference_t
	) const noexcept override;

	//-------------------------------------//

	std::array<matrix3d, 6> // Hij, Hik, Hil, Hjk, Hjl, Hkl
	Hessian (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept override;


	std::array<matrix3d, 6>
	Hessian (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		FiniteDifference_t
	) const noexcept override;

private:
	std::array<array3d, 3>
	ForceImp (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		_1stD_Func _1stDeriFunc
	) const noexcept;


	std::array<matrix3d, 6>
	HessianImp (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		_1stD_Func _1stDeriFunc,
		_2ndD_Func _2ndDeriFunc
	) const noexcept;

	//-------------------------------------//
	//-------------------------------------//
	//-------------------------------------//

	virtual double
	ObjectiveFunction (
		double /*cosPhi*/,
		double /*rij_*/,
		double /*rjk_*/,
		double /*rkl_*/,
		double /*cosThetaIJK*/,
		double /*cosThetaJKL*/
	) const noexcept = 0;


	virtual _1stDerivative_t // dcosPhi, drij, drik, drkl, dcosThetaIJK, dcosThetaJKL
	_1stDerivative (
		double cosPhi,
		double rij_,
		double rjk_,
		double rkl_,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;


	virtual _2ndDerivative_t
	_2ndDerivative (
		double cosPhi,
		double rij_,
		double rjk_,
		double rkl_,
		double cosThetaIJK,
		double cosThetaJKL
	) const noexcept;
};

#endif // POTENTIALDIHEDRAL_H_INCLUDED

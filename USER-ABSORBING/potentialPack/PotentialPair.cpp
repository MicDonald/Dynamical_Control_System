#include "PotentialPair.h"

using namespace std;

constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


double
PotentialPair::Energy (
	const array3d& rij
) const noexcept
{
	return ObjectiveFunction( Norm(rij) );
}

//---------------------------------------------------------------------------//

array3d
PotentialPair::Force (
	const array3d& rij
) const noexcept
{
	auto _1stDeriFunc = [this](double rij) noexcept {
		return this->_1stDerivative(rij);
	};
	return ForceImp(rij, _1stDeriFunc);
}


array3d
PotentialPair::Force (
	const array3d& rij,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_) noexcept {
		return this->PotentialPair::_1stDerivative(rij_);
	};
	return ForceImp(rij, _1stDeriFunc);
}


//---------------------------------------------------------------------------//

matrix3d
PotentialPair::Hessian (
	const array3d& rij
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_) noexcept {
		return this->_1stDerivative(rij_);
	};
	auto _2ndDeriFunc = [this](double rij_) noexcept {
		return this->_2ndDerivative(rij_);
	};
	return HessianImp(rij, _1stDeriFunc, _2ndDeriFunc);
}


matrix3d
PotentialPair::Hessian (
	const array3d& rij,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_) noexcept {
		return this->PotentialPair::_1stDerivative(rij_);
	};
	auto _2ndDeriFunc = [this](double rij_) noexcept {
		return this->PotentialPair::_2ndDerivative(rij_);
	};
	return HessianImp(rij, _1stDeriFunc, _2ndDeriFunc);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

array3d
PotentialPair::ForceImp (
	const array3d& rij,
	_1stD_Func _1stDeriFunc
) const noexcept
{
	auto rij_ = Norm(rij);
	auto _1stDeri = _1stDeriFunc(rij_);
	return Scale( UnitVector(rij, 1./rij_), _1stDeri );
}

//---------------------------------------------------------------------------//

matrix3d
PotentialPair::HessianImp (
	const array3d& rij,
	_1stD_Func _1stDeriFunc,
	_2ndD_Func _2ndDeriFunc
) const noexcept
{
	auto rij_ = Norm(rij);
	auto inverseRij = 1./rij_;
	auto nij = UnitVector(rij, inverseRij);
	auto _1stDeri = _1stDeriFunc(rij_);
	auto _2ndDeri = _2ndDeriFunc(rij_);
	auto nijnij = OuterDot(nij, nij);
	return Addition (
		Scale(nijnij, -_2ndDeri),
		Scale(Minus(nijnij, Identity()), _1stDeri*inverseRij)
	);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

PotentialPair::_1stDerivative_t
PotentialPair::_1stDerivative (
	double rij_
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(rij_);
	auto incrementObjFun = ObjectiveFunction(rij_+FiniteDistance);
	return (incrementObjFun-currentObjFun) * InverseFiniteDistance;
}


PotentialPair::_2ndDerivative_t
PotentialPair::_2ndDerivative (
	double rij_
) const noexcept
{
	auto current1stDeri = _1stDerivative(rij_);
	auto increment1stDeri = _1stDerivative(rij_+FiniteDistance);
	return (increment1stDeri-current1stDeri) * InverseFiniteDistance;
}


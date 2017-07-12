#include "PotentialAngle.h"
#include <iostream>

using namespace std;

constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


double
PotentialAngle::Energy (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	auto rij_ = Norm(rij);
	auto rik_ = Norm(rik);
	auto cosTheta = CosTheta (
		UnitVector(rij, 1./rij_),
		UnitVector(rik, 1./rik_)
	);
	return ObjectiveFunction( cosTheta, rij_, rik_ );	
}

//---------------------------------------------------------------------------//

std::array<array3d, 2>
PotentialAngle::Force (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) noexcept {
		return this->_1stDerivative(cosTheta, rij_, rik_);
	};
	return ForceImp(rij, rik, _1stDeriFunc);
}


std::array<array3d, 2>
PotentialAngle::Force (
	const array3d& rij,
	const array3d& rik,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) noexcept {
		return this->PotentialAngle::_1stDerivative(cosTheta, rij_, rik_);
	};
	return ForceImp(rij, rik, _1stDeriFunc);
}

//---------------------------------------------------------------------------//

std::array<matrix3d, 3>
PotentialAngle::Hessian (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) noexcept {
		return this->_1stDerivative(cosTheta, rij_, rik_);
	};
	auto _2ndDeriFunc = [this](double cosTheta, double rij_, double rik_) noexcept {
		return this->_2ndDerivative(cosTheta, rij_, rik_);
	};
	return HessianImp(rij, rik, _1stDeriFunc, _2ndDeriFunc);
}


std::array<matrix3d, 3>
PotentialAngle::Hessian (
	const array3d& rij,
	const array3d& rik,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosTheta, double rij_, double rik_) noexcept {
		return this->PotentialAngle::_1stDerivative(cosTheta, rij_, rik_);
	};
	auto _2ndDeriFunc = [this](double cosTheta, double rij_, double rik_) noexcept {
		return this->PotentialAngle::_2ndDerivative(cosTheta, rij_, rik_);
	};
	return HessianImp(rij, rik, _1stDeriFunc, _2ndDeriFunc);
}


void
PotentialAngle::_1stDebug (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	using namespace Angle_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij = 1./rij_;
	auto rik_ = Norm(rik);
	auto inverseRik = 1./rik_;
	auto nij = UnitVector(rij, inverseRij);
	auto nik = UnitVector(rik, inverseRik);
	auto cosTheta = CosTheta(nij, nik);
	auto _1stDeri = _1stDerivative(cosTheta, rij_, rik_);
	std::cout << _1stDeri[DCosTheta] << " "
		<< _1stDeri[DRij] << " "
		<< _1stDeri[DRik] << std::endl;

	_1stDeri = PotentialAngle::_1stDerivative(cosTheta, rij_, rik_);
	std::cout << _1stDeri[DCosTheta] << " "
		<< _1stDeri[DRij] << " "
		<< _1stDeri[DRik] << std::endl;
}


void
PotentialAngle::_2ndDebug (
	const array3d& rij,
	const array3d& rik
) const noexcept
{
	using namespace Angle_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij = 1./rij_;
	auto rik_ = Norm(rik);
	auto inverseRik = 1./rik_;
	auto nij = UnitVector(rij, inverseRij);
	auto nik = UnitVector(rik, inverseRik);
	auto cosTheta = CosTheta(nij, nik);
	auto _2ndDeri = _2ndDerivative(cosTheta, rij_, rik_);
	std::cout << _2ndDeri[DCosTheta_DCosTheta] << " "
		<< _2ndDeri[DCosTheta_DRij] << " "
		<< _2ndDeri[DCosTheta_DRik] << " "
		<< _2ndDeri[DRij_DRij] << " "
		<< _2ndDeri[DRij_DRik] << " "
		<< _2ndDeri[DRik_DRik] << std::endl;

	_2ndDeri = PotentialAngle::_2ndDerivative(cosTheta, rij_, rik_);
	std::cout << _2ndDeri[DCosTheta_DCosTheta] << " "
		<< _2ndDeri[DCosTheta_DRij] << " "
		<< _2ndDeri[DCosTheta_DRik] << " "
		<< _2ndDeri[DRij_DRij] << " "
		<< _2ndDeri[DRij_DRik] << " "
		<< _2ndDeri[DRik_DRik] << std::endl;
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

std::array<array3d, 2>
PotentialAngle::ForceImp (
	const array3d& rij,
	const array3d& rik,
	_1stD_Func _1stDeriFunc
) const noexcept
{
	using namespace Angle_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij = 1./rij_;
	auto rik_ = Norm(rik);
	auto inverseRik = 1./rik_;
	auto nij = UnitVector(rij, inverseRij);
	auto nik = UnitVector(rik, inverseRik);
	auto cosTheta = CosTheta(nij, nik);
	auto _1stDeri = _1stDeriFunc(cosTheta, rij_, rik_);

	auto termIJ_1 = _1stDeri[DCosTheta] * inverseRij;
	auto termIJ_2 = _1stDeri[DRij] - termIJ_1*cosTheta;
	auto termIK_1 = _1stDeri[DCosTheta] * inverseRik;
	auto termIK_2 = _1stDeri[DRik] - termIK_1*cosTheta;
	auto forceJI = Addition (
		Scale(nik, -termIJ_1),
		Scale(nij, -termIJ_2)
	);
	auto forceKI = Addition (
		Scale(nij, -termIK_1),
		Scale(nik, -termIK_2)
	);
	return {{ move(forceJI), move(forceKI) }};
}

//---------------------------------------------------------------------------//

std::array<matrix3d, 3>
PotentialAngle::HessianImp (
	const array3d& rij,
	const array3d& rik,
	_1stD_Func _1stDeriFunc,
	_2ndD_Func _2ndDeriFunc
) const noexcept
{
	using namespace Angle_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij = 1. / rij_;
	auto rik_ = Norm(rik);
	auto inverseRik = 1. / rik_;
	auto nij = UnitVector(rij, inverseRij);
	auto nik = UnitVector(rik, inverseRik);
	auto cosTheta = CosTheta(nij, nik);
	auto nij_nij = OuterDot(nij, nij);
	auto nik_nik = OuterDot(nik, nik);

	auto cosTheta_drj = Scale( Addition(nik, Scale(nij, -cosTheta)), inverseRij );
	auto cosTheta_drk = Scale( Addition(nij, Scale(nik, -cosTheta)), inverseRik );
	auto cosTheta_dri = Scale( Addition( cosTheta_drj, cosTheta_drk ), -1. );

	auto cosTheta_dri_drj = Scale( Addition (
		OuterDot( nij, cosTheta_drj ),
		OuterDot( cosTheta_dri, -nij ),
		Scale( Minus(nik_nik, Identity()), inverseRik ),
		Scale( Minus(nij_nij, Identity()), -cosTheta*inverseRij )
	), inverseRij );
	auto cosTheta_drj_drk = Scale( Addition (
		Scale( Minus(Identity(), nij_nij), inverseRij ),
		OuterDot( cosTheta_drj, -nik )
	), inverseRik );
	auto cosTheta_dri_drk = Scale( Addition (
		OuterDot( nik, cosTheta_drk ),
		OuterDot( cosTheta_dri, -nik ),
		Scale( Minus(nij_nij, Identity()), inverseRij ),
		Scale( Minus(nik_nik, Identity()), -cosTheta*inverseRik )
	), inverseRik );

	auto rij_dri_drj = Scale( Minus(nij_nij, Identity()), inverseRij );
	auto rik_dri_drk = Scale( Minus(nik_nik, Identity()), inverseRik );

	//-------------------------------------//

	auto _1stDeri = _1stDeriFunc(cosTheta, rij_, rik_);
	auto _2ndDeri = _2ndDeriFunc(cosTheta, rij_, rik_);

	auto hessianIJ = Addition (
		OuterDot( Addition (
			Scale(nij, -_2ndDeri[DRij_DRij]),
			Scale(nik, -_2ndDeri[DRij_DRik]),
			Scale(cosTheta_dri, _2ndDeri[DCosTheta_DRij])
		), nij),
		Scale(rij_dri_drj, _1stDeri[DRij]),
		OuterDot( Addition (
			Scale(nij, -_2ndDeri[DCosTheta_DRij]),
			Scale(nik, -_2ndDeri[DCosTheta_DRik]),
			Scale(cosTheta_dri, _2ndDeri[DCosTheta_DCosTheta])
		), cosTheta_drj),
		Scale(cosTheta_dri_drj, _1stDeri[DCosTheta])
	);
	auto hessianIK = Addition (
		OuterDot( Addition (
			Scale(nij, -_2ndDeri[DRij_DRik]),
			Scale(nik, -_2ndDeri[DRik_DRik]),
			Scale(cosTheta_dri, _2ndDeri[DCosTheta_DRik])
		), nik),
		Scale(rik_dri_drk, _1stDeri[DRik]),
		OuterDot ( Addition (
			Scale(cosTheta_dri, _2ndDeri[DCosTheta_DCosTheta]),
			Scale(nij, -_2ndDeri[DCosTheta_DRij]),
			Scale(nik, -_2ndDeri[DCosTheta_DRik])
		), cosTheta_drk),
		Scale( cosTheta_dri_drk, _1stDeri[DCosTheta] )
	);
	auto hessianJK = Addition (
		OuterDot( Addition (
			Scale(nij, _2ndDeri[DRij_DRik]),
			Scale(cosTheta_drj, _2ndDeri[DCosTheta_DRik])
		), nik),
		OuterDot( Addition (
			Scale(cosTheta_drj, _2ndDeri[DCosTheta_DCosTheta]),
			Scale(nij, _2ndDeri[DCosTheta_DRij])
		), cosTheta_drk),
		Scale(cosTheta_drj_drk, _1stDeri[DCosTheta])
	);

	return { move(hessianIJ), move(hessianIK), move(hessianJK) };
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

PotentialAngle::_1stDerivative_t
PotentialAngle::_1stDerivative (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(cosTheta, rij_, rik_);
	auto incrementObjFun_CosTheta = ObjectiveFunction(cosTheta+FiniteDistance, rij_, rik_);
	auto incrementObjFun_Rij = ObjectiveFunction(cosTheta, rij_+FiniteDistance, rik_);
	auto incrementObjFun_Rik = ObjectiveFunction(cosTheta, rij_, rik_+FiniteDistance);
	return {
		(incrementObjFun_CosTheta-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rij-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rik-currentObjFun) * InverseFiniteDistance
	};
}


PotentialAngle::_2ndDerivative_t
PotentialAngle::_2ndDerivative (
	double cosTheta,
	double rij_,
	double rik_
) const noexcept
{
	auto current1stDeri = _1stDerivative(cosTheta, rij_, rik_);
	auto increment1stDeri_CosTheta = _1stDerivative(cosTheta+FiniteDistance, rij_, rik_);
	auto increment1stDeri_Rij = _1stDerivative(cosTheta, rij_+FiniteDistance, rik_);
	auto increment1stDeri_Rik = _1stDerivative(cosTheta, rij_, rik_+FiniteDistance);
	return {
		(increment1stDeri_CosTheta[0]-current1stDeri[0]) * InverseFiniteDistance,
		(increment1stDeri_CosTheta[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_CosTheta[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rij[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_Rij[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rik[2]-current1stDeri[2]) * InverseFiniteDistance
	};
}


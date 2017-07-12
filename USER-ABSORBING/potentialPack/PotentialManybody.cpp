#include "PotentialManybody.h"

using namespace std;

constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


double
PotentialManybody::Energy (
	const array3d& rij,
	const std::vector<array3d>& rik
) const noexcept
{
	auto rij_ = Norm(rij);
	auto nij = UnitVector(rij, 1./rij_);
	std::vector<double> rik_;
	std::vector<double> cosTheta;
	for (const auto& iRik : rik)
	{
		auto iRik_ = Norm(iRik);
		auto nik = UnitVector(iRik, 1./iRik_);
		rik_.emplace_back(iRik_);
		cosTheta.emplace_back( CosTheta(nij, nik) );
	}
	return ObjectiveFunction(rij_, rik_, cosTheta);
}

//---------------------------------------------------------------------------//

std::vector<array3d>
PotentialManybody::Force (
	const array3d& rij,
	const std::vector<array3d>& rik
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_, const std::vector<double>& rik_, const std::vector<double>& cosTheta) noexcept {
		return this->_1stDerivative(rij_, rik_, cosTheta);
	};
	return ForceImp(rij, rik, _1stDeriFunc);
}


std::vector<array3d>
PotentialManybody::Force (
	const array3d& rij,
	const std::vector<array3d>& rik,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_, const std::vector<double>& rik_, const std::vector<double>& cosTheta) noexcept {
		return this->PotentialManybody::_1stDerivative(rij_, rik_, cosTheta);
	};
	return ForceImp(rij, rik, _1stDeriFunc);
}

//---------------------------------------------------------------------------//

std::vector<matrix3d>
PotentialManybody::Hessian (
	const array3d& rij,
	const std::vector<array3d>& rik
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_, const std::vector<double>& rik_, const std::vector<double>& cosTheta) noexcept {
		return this->_1stDerivative(rij_, rik_, cosTheta);
	};
	auto _2ndDeriFunc = [this](double rij_, const std::vector<double>& rik_, const std::vector<double>& cosTheta) noexcept {
		return this->_2ndDerivative(rij_, rik_, cosTheta);
	};
	return HessianImp(rij, rik, _1stDeriFunc, _2ndDeriFunc);
}


std::vector<matrix3d>
PotentialManybody::Hessian (
	const array3d& rij,
	const std::vector<array3d>& rik,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double rij_, const std::vector<double>& rik_, const std::vector<double>& cosTheta) noexcept {
		return this->PotentialManybody::_1stDerivative(rij_, rik_, cosTheta);
	};
	auto _2ndDeriFunc = [this](double rij_, const std::vector<double>& rik_, const std::vector<double>& cosTheta) noexcept {
		return this->PotentialManybody::_2ndDerivative(rij_, rik_, cosTheta);
	};
	return HessianImp(rij, rik, _1stDeriFunc, _2ndDeriFunc);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

std::vector<array3d>
PotentialManybody::ForceImp (
	const array3d& rij,
	const std::vector<array3d>& rik,
	_1stD_Func _1stDeriFunc
) const noexcept
{
	auto rij_ = Norm(rij);
	auto _rij_ = 1./rij_;
	auto nij = UnitVector(rij, _rij_);
	std::vector<double> rik_;
	std::vector<array3d> nik;
	std::vector<double> cosTheta;
	for (const auto& iRik : rik)
	{
		auto iRik_ = Norm(iRik);
		auto iNik = UnitVector(iRik, 1./iRik_);
		nik.emplace_back(iNik);
		rik_.emplace_back(iRik_);
		cosTheta.emplace_back( CosTheta(nij, iNik) );
	}
	auto _1stDeri = _1stDeriFunc(rij_, rik_, cosTheta);

	std::vector<array3d> forces;
	auto _1stDeriIter = _1stDeri.begin();
	forces.emplace_back( Scale(nij, _1stDeriIter->first) );
	++_1stDeriIter;
	auto nikIter = nik.begin();
	for (	auto rik_Iter=rik_.begin(), cosThetaIter=cosTheta.begin();
		rik_Iter!=rik_.end();
		++nikIter, ++rik_Iter, ++cosThetaIter, ++_1stDeriIter
	)
	{
		auto factorIJ = _rij_ * _1stDeriIter->second;
		Increment( forces.front(),
			Scale(nij, -*cosThetaIter * factorIJ),
			Scale(*nikIter, factorIJ)
		);
		auto factorIK = _1stDeriIter->second / *rik_Iter;
		forces.emplace_back( Addition (
			Scale(nij, factorIK),
			Scale(*nikIter, _1stDeriIter->first - factorIK * *cosThetaIter)
		) );
	}
	return forces;
}

//---------------------------------------------------------------------------//

std::vector<matrix3d>
PotentialManybody::HessianImp (
	const array3d& rij,
	const std::vector<array3d>& rik,
	_1stD_Func _1stDeriFunc,
	_2ndD_Func _2ndDeriFunc
) const noexcept
{
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

PotentialManybody::_1stDerivative_t
PotentialManybody::_1stDerivative (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(rij_, rik_, cosTheta);
	_1stDerivative_t _1stDeri;
	_1stDeri.emplace_back(
		(ObjectiveFunction(rij_+FiniteDistance, rik_, cosTheta)-currentObjFun) * InverseFiniteDistance,
		0.
	);

	auto copyRik_(rik_);
	auto copyCosTheta(cosTheta);
	for ( auto rik_Iter=copyRik_.begin(), cosThetaIter=copyCosTheta.begin(); rik_Iter!=copyRik_.end(); ++rik_Iter, ++cosThetaIter )
	{
		*rik_Iter += FiniteDistance;
		*cosThetaIter += FiniteDistance;
		_1stDeri.emplace_back( make_pair (
			(ObjectiveFunction(rij_, copyRik_, cosTheta)-currentObjFun) * InverseFiniteDistance,
			(ObjectiveFunction(rij_, rik_, copyCosTheta)-currentObjFun) * InverseFiniteDistance
		));
		*rik_Iter -= FiniteDistance;
		*cosThetaIter -= FiniteDistance;
	}
	return _1stDeri;
}


PotentialManybody::_2ndDerivative_t
PotentialManybody::_2ndDerivative (
	double rij_,
	const std::vector<double>& rik_,
	const std::vector<double>& cosTheta
) const noexcept
{
	_2ndDerivative_t _2ndDeri(rik_.size());
	auto current1stDeri = _1stDerivative(rij_, rik_, cosTheta);
	auto increment1stDeri_Rij = _1stDerivative(rij_+FiniteDistance, rik_, cosTheta);
	_2ndDeri.drij_drij = (increment1stDeri_Rij.front().first-current1stDeri.front().first) * InverseFiniteDistance;
	auto drij_drikIter = _2ndDeri.drij_drik.begin();
	auto drij_dcosTheta = _2ndDeri.drij_dcosTheta.begin();
	for (	auto currentIter=current1stDeri.begin()+1, incrementIter=increment1stDeri_Rij.begin()+1;
		currentIter!=current1stDeri.end();
		++currentIter, ++incrementIter, ++drij_drikIter, ++drij_dcosTheta
	)
	{
		*drij_drikIter = (incrementIter->first-currentIter->first) * InverseFiniteDistance;
		*drij_dcosTheta = (incrementIter->second-currentIter->second) * InverseFiniteDistance;
	}

	auto copyRik_(rik_);
	auto drik_drikOuterIter = _2ndDeri.drik_drik.begin();
	auto drik_dcosTOuterIter = _2ndDeri.drik_dcosTheta.begin();
	for ( auto& iRik : copyRik_ )
	{
		iRik += FiniteDistance;
		auto increment1st = _1stDerivative(rij_, copyRik_, cosTheta);
		auto drik_drikInnerIter = drik_drikOuterIter->begin();
		auto drik_dcosTInnerIter = drik_dcosTOuterIter->begin();
		for ( auto currentIter=current1stDeri.begin()+1, incrementIter=increment1st.begin()+1;
			currentIter!=current1stDeri.end();
			++currentIter, ++incrementIter, ++drik_drikInnerIter, ++drik_dcosTInnerIter
		)
		{
			*drik_drikInnerIter = (incrementIter->first-currentIter->first) * InverseFiniteDistance;
			*drik_dcosTInnerIter = (incrementIter->second-currentIter->second) * InverseFiniteDistance;
		}
		iRik -= FiniteDistance;
		++drik_drikOuterIter;
		++drik_dcosTOuterIter;
	}

	auto copyCosTheta(cosTheta);
	auto dcosT_dcosTOuterIter = _2ndDeri.dcosTheta_dcosTheta.begin();
	for ( auto& icosTheta : copyCosTheta )
	{
		icosTheta += FiniteDistance;
		auto increment1st = _1stDerivative(rij_, rik_, copyCosTheta);
		auto dcosT_dcosTInnerIter = dcosT_dcosTOuterIter->begin();
		for ( auto currentIter=current1stDeri.begin()+1, incrementIter=increment1st.begin()+1;
			currentIter!=current1stDeri.end();
			++currentIter, ++incrementIter, ++dcosT_dcosTInnerIter
		)
			*dcosT_dcosTInnerIter = (incrementIter->second-currentIter->second) * InverseFiniteDistance;
		icosTheta -= FiniteDistance;
		++dcosT_dcosTOuterIter;
	}
	return _2ndDeri;
}


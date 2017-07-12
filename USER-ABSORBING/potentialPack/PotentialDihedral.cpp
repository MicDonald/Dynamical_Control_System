#include "PotentialDihedral.h"

constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


double
PotentialDihedral::Energy (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto rij_ = Norm(rij);
	auto rjk_ = Norm(rjk);
	auto rkl_ = Norm(rkl);
	auto nij = UnitVector(rij, 1./rij_);
	auto njk = UnitVector(rjk, 1./rjk_);
	auto nkl = UnitVector(rkl, 1./rkl_);
	auto cosThetaIJK = -CosTheta(nij, njk);
	auto cosThetaJKL = -CosTheta(njk, nkl);
	auto cosPhi = CosPhi(rij, rjk, rkl);
	return ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
}

//---------------------------------------------------------------------------//

std::array<array3d, 3>
PotentialDihedral::Force (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return ForceImp(rij, rjk, rkl, _1stDeriFunc);
}


std::array<array3d, 3>
PotentialDihedral::Force (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->PotentialDihedral::_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return ForceImp(rij, rjk, rkl, _1stDeriFunc);
}

//---------------------------------------------------------------------------//

std::array<matrix3d, 6>
PotentialDihedral::Hessian (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	auto _2ndDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->_2ndDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return HessianImp(rij, rjk, rkl, _1stDeriFunc, _2ndDeriFunc);
}


std::array<matrix3d, 6>
PotentialDihedral::Hessian (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	FiniteDifference_t
) const noexcept
{
	auto _1stDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->PotentialDihedral::_1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	auto _2ndDeriFunc = [this](double cosPhi, double rij_, double rjk_, double rkl_, double cosThetaIJK, double cosThetaJKL) noexcept {
		return this->PotentialDihedral::_2ndDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	};
	return HessianImp(rij, rjk, rkl, _1stDeriFunc, _2ndDeriFunc);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

std::array<array3d, 3>
PotentialDihedral::ForceImp (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	_1stD_Func _1stDeriFunc
) const noexcept
{
	using namespace Dihedral_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij_ = 1./rij_;
	auto nij = UnitVector(rij, inverseRij_);
	auto rjk_ = Norm(rjk);
	auto inverseRjk_ = 1./rjk_;
	auto njk = UnitVector(rjk, inverseRjk_);
	auto rkl_ = Norm(rkl);
	auto inverseRkl_ = 1./rkl_;
	auto nkl = UnitVector(rkl, inverseRkl_);
	auto cosThetaIJK = -CosTheta(nij, njk);
	auto cosThetaJKL = -CosTheta(njk, nkl);
	auto nijk = Cross(rij, rjk);
	auto njkl = Cross(rjk, rkl);
	auto inverseNijk_ = 1./Norm(nijk);
	auto inverseNjkl_ = 1./Norm(njkl);
	auto cosPhi = CosTheta(UnitVector(nijk, inverseNijk_), UnitVector(njkl, inverseNjkl_));
	auto _1stDeri = _1stDeriFunc(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);

	auto inverseNijk_Sq = inverseNijk_*inverseNijk_;
	auto inverseNjkl_Sq = inverseNjkl_*inverseNjkl_;
	auto inverseNijk_Njkl_ = inverseNijk_*inverseNjkl_;
	auto rijCnijk = Cross(rij, nijk);
	auto rijCnjkl = Cross(rij, njkl);
	auto rjkCnijk = Cross(rjk, nijk);
	auto rjkCnjkl = Cross(rjk, njkl);
	auto rklCnijk = Cross(rkl, nijk);
	auto rklCnjkl = Cross(rkl, njkl);

	auto forceIJ = Addition (
		Scale(nij, _1stDeri[DRij] -
			_1stDeri[DCosThetaIJK]*cosThetaIJK*inverseRij_),
		Scale(njk, -_1stDeri[DCosThetaIJK]*inverseRij_),
		Scale(rjkCnijk, -_1stDeri[DCosPhi]*cosPhi*inverseNijk_Sq),
		Scale(rjkCnjkl, _1stDeri[DCosPhi]*inverseNijk_Njkl_)
	);
	auto forceKJ = Addition (
		Scale(nij, _1stDeri[DCosThetaIJK]*inverseRjk_),
		Scale(njk, -_1stDeri[DRjk] +
			(_1stDeri[DCosThetaIJK]*cosThetaIJK +
			_1stDeri[DCosThetaJKL]*cosThetaJKL) * inverseRjk_),
		Scale(nkl, _1stDeri[DCosThetaJKL]*inverseRjk_),
		Scale(rijCnijk, -_1stDeri[DCosPhi]*cosPhi*inverseNijk_Sq),
		Scale(rijCnjkl, _1stDeri[DCosPhi]*inverseNijk_Njkl_),
		Scale(rklCnijk, -_1stDeri[DCosPhi]*inverseNijk_Njkl_),
		Scale(rklCnjkl, _1stDeri[DCosPhi]*cosPhi*inverseNjkl_Sq)
	);
	auto forceLK = Addition (
		Scale(njk, _1stDeri[DCosThetaJKL]*inverseRkl_),
		Scale(nkl, -_1stDeri[DRkl] +
			_1stDeri[DCosThetaJKL]*cosThetaJKL*inverseRkl_),
		Scale(rjkCnijk, _1stDeri[DCosPhi]*inverseNijk_Njkl_),
		Scale(rjkCnjkl, -_1stDeri[DCosPhi]*cosPhi*inverseNjkl_Sq)
	);
	return {{ move(forceIJ), move(forceKJ), move(forceLK) }};
}

//---------------------------------------------------------------------------//

std::array<matrix3d, 6>
PotentialDihedral::HessianImp (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl,
	_1stD_Func _1stDeriFunc,
	_2ndD_Func _2ndDeriFunc
) const noexcept
{
	using namespace Dihedral_Namespace;

	auto rij_ = Norm(rij);
	auto inverseRij_ = 1./rij_;
	auto nij = UnitVector(rij, inverseRij_);
	auto rjk_ = Norm(rjk);
	auto rjk_Sq = rjk_*rjk_;
	auto inverseRjk_ = 1./rjk_;
	auto njk = UnitVector(rjk, inverseRjk_);
	auto rkl_ = Norm(rkl);
	auto rkl_Sq = rkl_*rkl_;
	auto inverseRkl_ = 1./rkl_;
	auto nkl = UnitVector(rkl, inverseRkl_);
	auto cosThetaIJK = -CosTheta(nij, njk);
	auto cosThetaJKL = -CosTheta(njk, nkl);
	auto nijk = Cross(rij, rjk);
	auto njkl = Cross(rjk, rkl);
	auto inverseNijk_ = 1./Norm(nijk);
	auto inverseNjkl_ = 1./Norm(njkl);
	auto cosPhi = CosTheta(UnitVector(nijk, inverseNijk_), UnitVector(njkl, inverseNjkl_));
	auto _1stDeri = _1stDeriFunc(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto _2ndDeri = _2ndDeriFunc(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);

	auto rik = Addition(rij, rjk);
	auto rjl = Addition(rjk, rkl);
	auto rijCnijk = Cross(rij, nijk);
	auto rjkCnijk = Cross(rjk, nijk);
	auto rklCnijk = Cross(rkl, nijk);
	auto rijCnjkl = Cross(rij, njkl);
	auto rjkCnjkl = Cross(rjk, njkl);
	auto rklCnjkl = Cross(rkl, njkl);
	auto inverseNijk_Njkl_ = inverseNijk_*inverseNjkl_;
	auto inverseNijk_Sq = inverseNijk_*inverseNijk_;
	auto inverseNjkl_Sq = inverseNjkl_*inverseNjkl_;
	auto cosPhi_dri = Cross (
		rjk,
		Minus (
			Scale(nijk, cosPhi*inverseNijk_Sq),
			Scale(njkl, inverseNijk_Njkl_)
		)
	);
	auto cosPhi_drkj = Addition (
		Scale( Minus (
			rklCnijk,
			rijCnjkl
		), inverseNijk_Njkl_),
		Scale(rijCnijk, cosPhi*inverseNijk_Sq),
		Scale(rklCnjkl, -cosPhi*inverseNjkl_Sq)
	);
	auto cosPhi_drj = Scale(Addition(cosPhi_dri, cosPhi_drkj), -1.);
	auto cosPhi_drl = Cross (
		rjk,
		Minus (
			Scale(njkl, cosPhi*inverseNjkl_Sq),
			Scale(nijk, inverseNijk_Njkl_)
		)
	);
	auto cosPhi_drk = Addition(Scale(cosPhi_drl, -1.), cosPhi_drkj);
	auto cosPhi_dri_drj = Addition (
		OuterDot (
			Scale(rjkCnijk, inverseNijk_Sq*inverseNijk_Njkl_),
			Minus ( Addition( rijCnjkl, rjkCnjkl),
				rklCnijk
			)
		),
		Scale( Addition (
			OuterDot(rjk, rkl),
			Scale(Identity(), Dot(rkl, rjk)),
			OuterDot(rkl, Scale(rjk, -2.))
		), inverseNijk_Njkl_),
		OuterDot (
			cosPhi_dri,
			Scale(rklCnjkl, inverseNjkl_Sq)
		),
		OuterDot (
			Addition (
				cosPhi_dri,
				Scale(rjkCnijk, 2.*cosPhi*inverseNijk_Sq)
			),
			Scale( Addition (
				rijCnijk, rjkCnijk
			), -inverseNijk_Sq)
		),
		Scale( Minus (
			Addition (
				OuterDot(rjk, rij),
				Scale(Identity(), Dot(rij, rjk)+rjk_Sq)
			),
			Addition (
				OuterDot(rij, Scale(rjk, 2.)),
				OuterDot(rjk, rjk)
			)
		), cosPhi*inverseNijk_Sq)
	);
	auto cosPhi_dri_drk = Addition (
		Scale( OuterDot (
			rjkCnijk,
			Minus( Addition (
				rklCnijk,
				rjkCnijk
			), rijCnjkl)
		), inverseNijk_Sq*inverseNijk_Njkl_),
		Scale( Minus (
			Addition (
				OuterDot(rjl, rjk),
				OuterDot(rkl, rjk)
			),
			Addition (
				OuterDot(rjk, rkl),
				Scale(Identity(), Dot(rjk,rkl)+rjk_Sq)
			)
		), inverseNijk_Njkl_),
		OuterDot (
			cosPhi_dri,
			Scale(rijCnijk, inverseNijk_Sq)
		),
		Scale( Minus (
			OuterDot(rij, Scale(rjk, 2.)),
			Addition (
				OuterDot(rjk, rij),
				Scale(Identity(), Dot(rij, rjk))
			)
		), cosPhi*inverseNijk_Sq),
		Scale( OuterDot(rjkCnijk, rijCnijk), 2.*cosPhi*inverseNijk_Sq*inverseNijk_Sq),
		Scale( OuterDot (
			cosPhi_dri,
			Addition(rjkCnjkl, rklCnjkl)
		), -inverseNjkl_Sq)
	);
	auto cosPhi_dri_drl = Addition (
		Scale( Minus (
			Scale(Identity(), rjk_Sq),
			OuterDot(rjk, rjk)
		), inverseNijk_Njkl_),
		OuterDot(
			Scale(rjkCnijk, -inverseNijk_Sq*inverseNijk_Njkl_),
		 	rjkCnijk
		),
		OuterDot (
			cosPhi_dri,
			Scale( rjkCnjkl, inverseNjkl_Sq )
		)
	);
	auto cosPhi_drj_drl = Addition (
		Scale( Minus (
			OuterDot(rjk, Addition(rij, rik)),
			Addition (
				OuterDot(rij, rjk),
				Scale(Identity(), Dot(rik, rjk))
			)
		), inverseNijk_Njkl_),
		Scale( OuterDot (
			Addition(rijCnijk, rjkCnijk),
			rjkCnijk
		), inverseNijk_Sq*inverseNijk_Njkl_),
		Scale( OuterDot(
			rklCnjkl,
			rjkCnijk
		), -inverseNijk_Njkl_*inverseNjkl_Sq),
		OuterDot(
			cosPhi_drj,
			Scale(rjkCnjkl, inverseNjkl_Sq)
		),
		Scale( Minus (
			OuterDot(rjk, Scale(rkl, 2.)),
			Addition (
				OuterDot(rkl, rjk),
				Scale(Identity(), Dot(rjk, rkl))
			)
		), cosPhi*inverseNjkl_Sq),
		Scale( OuterDot (
			rklCnjkl,
			rjkCnjkl
		), 2.*cosPhi*inverseNjkl_Sq*inverseNjkl_Sq)
	);
	auto cosPhi_drj_drk = Addition (
		Scale( Minus (
			Addition (
				OuterDot(rjk, rkl),
				Scale(Identity(), Dot(rkl, rik)+Dot(rkl, rij))
			),
			Addition (
				OuterDot(rkl, Addition(Scale(rjk, 2.), rij)),
				OuterDot(rij, rkl)
			)
		), inverseNijk_Njkl_),
		Scale( OuterDot (
			Minus (
				Scale(rklCnjkl, inverseNjkl_Sq),
				Scale( Addition(rijCnijk, rjkCnijk), inverseNijk_Sq )
			),
			Minus(rklCnijk, rijCnjkl)
		), inverseNijk_Njkl_),
		OuterDot( Minus (
			cosPhi_drj,
			Scale( Addition(rijCnijk, rjkCnijk), 2.*cosPhi*inverseNijk_Sq)
		), Scale(rijCnijk, inverseNijk_Sq)),
		Scale( Minus (
			Addition (
				OuterDot(rjk, rij),
				Scale(Identity(), Dot(rij, rik))
			),
			OuterDot(rij, Addition(Scale(rjk, 2.), rij))
		), cosPhi*inverseNijk_Sq),
		OuterDot( Addition (
			cosPhi_drj,
			Scale(rklCnjkl, 2.*cosPhi*inverseNjkl_Sq)
		), Scale(rklCnjkl, -inverseNjkl_Sq)),
		Scale( Minus (
			Scale(Identity(), rkl_Sq),
			OuterDot(rkl, rkl)
		), cosPhi*inverseNjkl_Sq),
		Scale(cosPhi_drj_drl, -1.)
	);
	auto cosPhi_drk_drl = Addition (
		Scale( Minus (
			Addition (
				OuterDot(rij, rjk),
				Scale(Identity(), Dot(rij, rjk))
			),
			OuterDot(rjk, Scale(rij, 2.))
		), inverseNijk_Njkl_),
		Scale( OuterDot(
			rijCnijk,
			rjkCnijk
		), -inverseNijk_Sq*inverseNijk_Njkl_),
		Scale( OuterDot (
			Addition(rjkCnjkl, rklCnjkl),
			rjkCnijk
		), inverseNijk_Njkl_*inverseNjkl_Sq),
		Scale( OuterDot (
			cosPhi_drk,
			rjkCnjkl
		), inverseNjkl_Sq),
		Scale( OuterDot (
			Addition (rjkCnjkl, rklCnjkl),
			rjkCnjkl
		), -2.*cosPhi*inverseNjkl_Sq*inverseNjkl_Sq),
		Scale( Minus (
			Addition (
				OuterDot(rkl, rjk),
				Scale(Identity(), Dot(rjk, rjl))
			),
			Addition (
				OuterDot(rjk, rjl),
				OuterDot(rjk, rkl)
			)
		), cosPhi*inverseNjkl_Sq)
	);
	auto cosThetaIJK_dri = Addition (
		Scale(nij, cosThetaIJK*inverseRij_),
		Scale(njk, -inverseRij_)
	);
	auto cosThetaIJK_drj = Addition (
		Scale(nij, inverseRjk_-cosThetaIJK*inverseRij_),
		Scale(njk, cosThetaIJK*inverseRjk_-inverseRij_)
	);
	auto cosThetaIJK_drk = Addition (
		Scale(nij, -inverseRjk_),
		Scale(njk, -cosThetaIJK*inverseRjk_)
	);
	auto cosThetaJKL_drj = Addition (
		Scale(nkl, inverseRjk_),
		Scale(njk, cosThetaJKL*inverseRjk_)
	);
	auto cosThetaJKL_drk = Addition (
		Scale(njk, -cosThetaJKL*inverseRjk_+inverseRkl_),
		Scale(nkl, -inverseRjk_+cosThetaJKL*inverseRkl_)
	);
	auto cosThetaJKL_drl = Addition (
		Scale(njk, -inverseRkl_),
		Scale(nkl, -cosThetaJKL*inverseRkl_)
	);
	auto nijOnij = OuterDot(nij, nij);
	auto njkOnjk = OuterDot(njk, njk);
	auto nklOnkl = OuterDot(nkl, nkl);
	auto nij_dri = Scale( Minus(nijOnij, Identity()), inverseRij_);
	auto njk_drj = Scale( Minus(njkOnjk, Identity()), inverseRjk_);
	auto nkl_drk = Scale( Minus(nklOnkl, Identity()), inverseRkl_);
	auto cosThetaIJK_dri_drj = Scale( Addition (
		OuterDot(cosThetaIJK_dri, nij),
		nij_dri,
		Scale( OuterDot(nij, cosThetaIJK_drj), -1.)
	), -inverseRij_);
	auto cosThetaIJK_dri_drk = Scale( Minus (
		nij_dri,
		OuterDot(cosThetaIJK_drj, njk)
	), inverseRjk_);
	auto cosThetaIJK_drj_drk = Scale( Minus (
		Addition (
			nij_dri,
			OuterDot(njk, cosThetaIJK_drk)
		),
		Addition (
			OuterDot(cosThetaIJK_drj, njk),
			Scale(njk_drj, cosThetaIJK)
		)
	), inverseRjk_);
	auto cosThetaJKL_drj_drk = Addition (
		Scale( Addition (
			OuterDot(njk, cosThetaJKL_drj),
			OuterDot(cosThetaJKL_drj, njk),
			Scale(njk_drj, cosThetaJKL)
		), -inverseRjk_ ),
		Scale( Addition (
			njk_drj,
			OuterDot(cosThetaJKL_drj, nkl)
		), inverseRkl_ )
	);
	auto cosThetaJKL_drj_drl = Scale( Addition (
		njk_drj,
		OuterDot(cosThetaJKL_drj, nkl)
	), -inverseRkl_);
	auto cosThetaJKL_drk_drl = Scale( Addition (
		OuterDot(nkl, cosThetaJKL_drl),
		OuterDot(cosThetaJKL_drl, nkl),
		Scale(nkl_drk, -cosThetaJKL)
	), -inverseRkl_);

	matrix3d hessianIJ = Addition (
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, -_2ndDeri[DCosPhi_DRij]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosPhi_DCosThetaIJK])
		), cosPhi_drj),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DRij]),
			Scale(nij, -_2ndDeri[DRij_DRij]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRij_DCosThetaIJK])
		), nij),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DRjk]),
			Scale(nij, -_2ndDeri[DRij_DRjk]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRjk_DCosThetaIJK])
		), Scale(njk, -1.)),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosThetaIJK]),
			Scale(nij, -_2ndDeri[DRij_DCosThetaIJK]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaIJK])
		), cosThetaIJK_drj),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(nij, -_2ndDeri[DRij_DCosThetaJKL]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaJKL])
		), cosThetaJKL_drj),
		Scale(cosPhi_dri_drj, _1stDeri[DCosPhi]),
		Scale(nij_dri, _1stDeri[DRij]),
		Scale(cosThetaIJK_dri_drj, _1stDeri[DCosThetaIJK])
	);
	matrix3d hessianIK = Addition (
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, -_2ndDeri[DCosPhi_DRij]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosPhi_DCosThetaIJK])
		), cosPhi_drk),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DRjk]),
			Scale(nij, -_2ndDeri[DRij_DRjk]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRjk_DCosThetaIJK])
		), njk),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DRkl]),
			Scale(nij, -_2ndDeri[DRij_DRkl]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRkl_DCosThetaIJK])
		), Scale(nkl, -1.)),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosThetaIJK]),
			Scale(nij, -_2ndDeri[DRij_DCosThetaIJK]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaIJK])
		), cosThetaIJK_drk),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(nij, -_2ndDeri[DRij_DCosThetaJKL]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaJKL])
		), cosThetaJKL_drk),
		Scale(cosPhi_dri_drk, _1stDeri[DCosPhi]),
		Scale(cosThetaIJK_dri_drk, _1stDeri[DCosThetaIJK])
	);
	matrix3d hessianIL = Addition (
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, -_2ndDeri[DCosPhi_DRij]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosPhi_DCosThetaIJK])
		), cosPhi_drl),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DRkl]),
			Scale(nij, -_2ndDeri[DRij_DRkl]),
			Scale(cosThetaIJK_dri, _2ndDeri[DRkl_DCosThetaIJK])
		), nkl),
		OuterDot( Addition (
			Scale(cosPhi_dri, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(nij, -_2ndDeri[DRij_DCosThetaJKL]),
			Scale(cosThetaIJK_dri, _2ndDeri[DCosThetaIJK_DCosThetaJKL])
		), cosThetaJKL_drl),
		Scale(cosPhi_dri_drl, _1stDeri[DCosPhi])
	);
	matrix3d hessianJK = Addition (
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, _2ndDeri[DCosPhi_DRij]),
			Scale(njk, -_2ndDeri[DCosPhi_DRjk]),
			Scale(cosThetaIJK_drj, _2ndDeri[DCosPhi_DCosThetaIJK]),
			Scale(cosThetaJKL_drj, _2ndDeri[DCosPhi_DCosThetaJKL])
		), cosPhi_drk),
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DRjk]),
			Scale(nij, _2ndDeri[DRij_DRjk]),
			Scale(njk, -_2ndDeri[DRjk_DRjk]),
			Scale(cosThetaIJK_drj, _2ndDeri[DRjk_DCosThetaIJK]),
			Scale(cosThetaJKL_drj, _2ndDeri[DRjk_DCosThetaJKL])
		), njk),
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DRkl]),
			Scale(nij, _2ndDeri[DRij_DRkl]),
			Scale(njk, -_2ndDeri[DRjk_DRkl]),
			Scale(cosThetaIJK_drj, _2ndDeri[DRkl_DCosThetaIJK]),
			Scale(cosThetaJKL_drj, _2ndDeri[DRkl_DCosThetaJKL])
		), Scale(nkl, -1.)),
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DCosThetaIJK]),
			Scale(nij, _2ndDeri[DRij_DCosThetaIJK]),
			Scale(njk, -_2ndDeri[DRjk_DCosThetaIJK]),
			Scale(cosThetaIJK_drj, _2ndDeri[DCosThetaIJK_DCosThetaIJK]),
			Scale(cosThetaJKL_drj, _2ndDeri[DCosThetaIJK_DCosThetaJKL])
		), cosThetaIJK_drk),
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(nij, _2ndDeri[DRij_DCosThetaJKL]),
			Scale(njk, -_2ndDeri[DRjk_DCosThetaJKL]),
			Scale(cosThetaIJK_drj, _2ndDeri[DCosThetaIJK_DCosThetaJKL]),
			Scale(cosThetaJKL_drj, _2ndDeri[DCosThetaJKL_DCosThetaJKL])
		), cosThetaJKL_drk),
		Scale(cosPhi_drj_drk, _1stDeri[DCosPhi]),
		Scale(njk_drj, _1stDeri[DRjk]),
		Scale(cosThetaIJK_drj_drk, _1stDeri[DCosThetaIJK]),
		Scale(cosThetaJKL_drj_drk, _1stDeri[DCosThetaJKL])
	);
	matrix3d hessianJL = Addition (
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(nij, _2ndDeri[DCosPhi_DRij]),
			Scale(njk, -_2ndDeri[DCosPhi_DRjk]),
			Scale(cosThetaIJK_drj, _2ndDeri[DCosPhi_DCosThetaIJK]),
			Scale(cosThetaJKL_drj, _2ndDeri[DCosPhi_DCosThetaJKL])
		), cosPhi_drl),
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DRkl]),
			Scale(nij, _2ndDeri[DRij_DRkl]),
			Scale(njk, -_2ndDeri[DRjk_DRkl]),
			Scale(cosThetaIJK_drj, _2ndDeri[DRkl_DCosThetaIJK]),
			Scale(cosThetaJKL_drj, _2ndDeri[DRkl_DCosThetaJKL])
		), nkl),
		OuterDot( Addition (
			Scale(cosPhi_drj, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(nij, _2ndDeri[DRij_DCosThetaJKL]),
			Scale(njk, -_2ndDeri[DRjk_DCosThetaJKL]),
			Scale(cosThetaIJK_drj, _2ndDeri[DCosThetaIJK_DCosThetaJKL]),
			Scale(cosThetaJKL_drj, _2ndDeri[DCosThetaJKL_DCosThetaJKL])
		), cosThetaJKL_drl),
		Scale(cosPhi_drj_drl, _1stDeri[DCosPhi]),
		Scale(cosThetaJKL_drj_drl, _1stDeri[DCosThetaJKL])
	);
	matrix3d hessianKL = Addition (
		OuterDot( Addition (
			Scale(cosPhi_drk, _2ndDeri[DCosPhi_DCosPhi]),
			Scale(njk, _2ndDeri[DCosPhi_DRjk]),
			Scale(nkl, -_2ndDeri[DCosPhi_DRkl]),
			Scale(cosThetaIJK_drk, _2ndDeri[DCosPhi_DCosThetaIJK]),
			Scale(cosThetaJKL_drk, _2ndDeri[DCosPhi_DCosThetaJKL])
		), cosPhi_drl),
		OuterDot( Addition (
			Scale(cosPhi_drk, _2ndDeri[DCosPhi_DRkl]),
			Scale(njk, _2ndDeri[DRjk_DRkl]),
			Scale(nkl, -_2ndDeri[DRkl_DRkl]),
			Scale(cosThetaIJK_drk, _2ndDeri[DRkl_DCosThetaIJK]),
			Scale(cosThetaJKL_drk, _2ndDeri[DRkl_DCosThetaJKL])
		), nkl),
		OuterDot( Addition (
			Scale(cosPhi_drk, _2ndDeri[DCosPhi_DCosThetaJKL]),
			Scale(njk, _2ndDeri[DRjk_DCosThetaJKL]),
			Scale(nkl, -_2ndDeri[DRkl_DCosThetaJKL]),
			Scale(cosThetaIJK_drk, _2ndDeri[DCosThetaIJK_DCosThetaJKL]),
			Scale(cosThetaJKL_drk, _2ndDeri[DCosThetaJKL_DCosThetaJKL])
		), cosThetaJKL_drl),
		Scale(cosPhi_drk_drl, _1stDeri[DCosPhi]),
		Scale(nkl_drk, _1stDeri[DRkl]),
		Scale(cosThetaJKL_drk_drl, _1stDeri[DCosThetaJKL])
	);
	return {{
		move(hessianIJ), move(hessianIK), move(hessianIL),
		move(hessianJK), move(hessianJL),
		move(hessianKL)
	}};
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

PotentialDihedral::_1stDerivative_t
PotentialDihedral::_1stDerivative (
	double cosPhi,
	double rij_,
	double rjk_,
	double rkl_,
	double cosThetaIJK,
	double cosThetaJKL
) const noexcept
{
	auto currentObjFun = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_CosPhi = ObjectiveFunction(cosPhi+FiniteDistance, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rij = ObjectiveFunction(cosPhi, rij_+FiniteDistance, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rjk = ObjectiveFunction(cosPhi, rij_, rjk_+FiniteDistance, rkl_, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_Rkl = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_+FiniteDistance, cosThetaIJK, cosThetaJKL);
	auto incrementObjFun_CosThetaIJK = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK+FiniteDistance, cosThetaJKL);
	auto incrementObjFun_CosThetaJKL = ObjectiveFunction(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL+FiniteDistance);
	return {
		(incrementObjFun_CosPhi-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rij-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rjk-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_Rkl-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_CosThetaIJK-currentObjFun) * InverseFiniteDistance,
		(incrementObjFun_CosThetaJKL-currentObjFun) * InverseFiniteDistance,
	};
}

PotentialDihedral::_2ndDerivative_t
PotentialDihedral::_2ndDerivative (
	double cosPhi,
	double rij_,
	double rjk_,
	double rkl_,
	double cosThetaIJK,
	double cosThetaJKL
) const noexcept
{
	auto current1stDeri = _1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_CosPhi = _1stDerivative(cosPhi+FiniteDistance, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rij = _1stDerivative(cosPhi, rij_+FiniteDistance, rjk_, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rjk = _1stDerivative(cosPhi, rij_, rjk_+FiniteDistance, rkl_, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_Rkl = _1stDerivative(cosPhi, rij_, rjk_, rkl_+FiniteDistance, cosThetaIJK, cosThetaJKL);
	auto increment1stDeri_CosThetaIJK = _1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK+FiniteDistance, cosThetaJKL);
	auto increment1stDeri_CosThetaJKL = _1stDerivative(cosPhi, rij_, rjk_, rkl_, cosThetaIJK, cosThetaJKL+FiniteDistance);
	return {
		(increment1stDeri_CosPhi[0]-current1stDeri[0]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_CosPhi[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_Rij[1]-current1stDeri[1]) * InverseFiniteDistance,
		(increment1stDeri_Rij[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rij[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_Rij[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_Rij[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[2]-current1stDeri[2]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_Rjk[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_Rkl[3]-current1stDeri[3]) * InverseFiniteDistance,
		(increment1stDeri_Rkl[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_Rkl[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_CosThetaIJK[4]-current1stDeri[4]) * InverseFiniteDistance,
		(increment1stDeri_CosThetaIJK[5]-current1stDeri[5]) * InverseFiniteDistance,
		(increment1stDeri_CosThetaJKL[5]-current1stDeri[5]) * InverseFiniteDistance,
	};
}


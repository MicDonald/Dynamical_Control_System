#ifndef POTENTIALFLAGS_H_INCLUDED
#define POTENTIALFLAGS_H_INCLUDED

#include <array>
#include <vector>
#include <map>

constexpr double PI = 3.14159265358979323846;


struct FiniteDifference_t {};


namespace Angle_Namespace {
	enum {DCosTheta, DRij, DRik};
	enum {DCosTheta_DCosTheta, DCosTheta_DRij, DCosTheta_DRik, DRij_DRij, DRij_DRik, DRik_DRik};
	using _1stDerivative_t = std::array<double, 3>;
	using _2ndDerivative_t = std::array<double, 6>;
};


namespace Dihedral_Namespace {
	enum {DCosPhi, DRij, DRjk, DRkl, DCosThetaIJK, DCosThetaJKL};
	enum {	DCosPhi_DCosPhi, DCosPhi_DRij, DCosPhi_DRjk, DCosPhi_DRkl, DCosPhi_DCosThetaIJK, DCosPhi_DCosThetaJKL,
		DRij_DRij, DRij_DRjk, DRij_DRkl, DRij_DCosThetaIJK, DRij_DCosThetaJKL,
		DRjk_DRjk, DRjk_DRkl, DRjk_DCosThetaIJK, DRjk_DCosThetaJKL,
		DRkl_DRkl, DRkl_DCosThetaIJK, DRkl_DCosThetaJKL,
		DCosThetaIJK_DCosThetaIJK, DCosThetaIJK_DCosThetaJKL,
		DCosThetaJKL_DCosThetaJKL
	};
	using _1stDerivative_t = std::array<double, 6>;
	using _2ndDerivative_t = std::array<double, 21>;
};


namespace Manybody_Namespace {
	using _1stDerivative_t = std::vector<std::pair<double,double>>;
	struct _2ndDerivative_t {
		double drij_drij;
		std::vector<double> drij_drik;
		std::vector<double> drij_dcosTheta;
		std::vector<std::vector<double>> drik_drik;
		std::vector<std::vector<double>> drik_dcosTheta;
		std::vector<std::vector<double>> dcosTheta_dcosTheta;

		_2ndDerivative_t (
			unsigned n
		) noexcept :
			drij_drik(n),
			drij_dcosTheta(n),
			drik_drik(n),
			drik_dcosTheta(n),
			dcosTheta_dcosTheta(n)
		{
			for ( auto& ik_ik : drik_drik )
				ik_ik.resize(n);
			for ( auto& ik_cosT : drik_dcosTheta )
				ik_cosT.resize(n);
			for ( auto& cosT_cosT : dcosTheta_dcosTheta )
				cosT_cosT.resize(n);
		}
	};
};

#endif // POTENTIALFLAGS_H_INCLUDED

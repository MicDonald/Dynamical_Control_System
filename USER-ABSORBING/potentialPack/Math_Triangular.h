#include <utility>

struct Math_Triangular {
	//notice, all the derivatives are about cosTheta

	double
	_1stDTheta (
		double cosTheta
	) const noexcept;


	double
	_2ndDTheta (
		double cosTheta
	) const noexcept;

	//-------------------------------------//

	double
	CosNTheta (
		unsigned n,
		double cosTheta
	) const noexcept;


	double
	_1stDCosNTheta (
		unsigned n,
		double cosTheta
	) const noexcept;


	double
	_2ndDCosNTheta (
		unsigned n,
		double cosTheta
	) const noexcept;

	//-------------------------------------//

	double
	CosNTheta_delta (
		unsigned n,
		double cosTheta,
		double delta
	) const noexcept;


	double
	_1stDCosNTheta_delta (
		unsigned n,
		double cosTheta,
		double delta
	) const noexcept;


	double
	_2ndDCosNTheta_delta (
		unsigned n,
		double cosTheta,
		double delta
	) const noexcept;

	//-------------------------------------//

	std::pair<double, double> // cosNTheta, sinNTheta
	TriNTheta (
		unsigned n,
		double cosTheta
	) const noexcept;


	std::pair<double, double>
	_1stDTriNTheta (
		unsigned n,
		double cosTheta
	) const noexcept;


	std::pair<double, double>
	_2ndDTriNTheta (
		unsigned n,
		double cosTheta
	) const noexcept;
};

//---------------------------------------------------------------------------//

inline double
_1stDTheta (
	double cosTheta
) noexcept
{
	return Math_Triangular()._1stDTheta(cosTheta);
}


inline double
_2ndDTheta (
	double cosTheta
) noexcept
{
	return Math_Triangular()._2ndDTheta(cosTheta);
}


inline double
CosNTheta (
	unsigned n,
	double cosTheta
) noexcept
{
	return Math_Triangular().CosNTheta(n, cosTheta);
}


inline double
_1stDCosNTheta (
	unsigned n,
	double cosTheta
) noexcept
{
	return Math_Triangular()._1stDCosNTheta(n, cosTheta);
}


inline double
_2ndDCosNTheta (
	unsigned n,
	double cosTheta
) noexcept
{
	return Math_Triangular()._2ndDCosNTheta(n, cosTheta);
}


inline double
CosNTheta_delta (
	unsigned n,
	double cosTheta,
	double delta
) noexcept
{
	return Math_Triangular().CosNTheta_delta(n, cosTheta, delta);
}


inline double
_1stDCosNTheta_delta (
	unsigned n,
	double cosTheta,
	double delta
) noexcept
{
	return Math_Triangular()._1stDCosNTheta_delta(n, cosTheta, delta);
}


inline double
_2ndDCosNTheta_delta (
	unsigned n,
	double cosTheta,
	double delta
) noexcept
{
	return Math_Triangular()._2ndDCosNTheta_delta(n, cosTheta, delta);
}


#ifndef DATA3D_H_INCLUDED
#define DATA3D_H_INCLUDED

#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <functional>


using array3d = std::array<double, 3>;
using matrix3d = std::array<array3d, 3>;
namespace CartesianCoordinate_Namespace {
	enum {X, Y, Z};
};


inline void
Increment (
	array3d& a,
	const array3d& b
) noexcept
{
	std::transform (
		a.begin(), a.end(),
		b.begin(),
		a.begin(),
		std::plus<double>()
	);
}


inline array3d
Minus (
	const array3d& to,
	const array3d& from
) noexcept
{
	array3d diff;
	std::transform (
		to.begin(), to.end(),
		from.begin(),
		diff.begin(),
		std::minus<double>()
	);
	return diff;
}


inline array3d
Scale (
	const array3d& data,
	double factor
) noexcept
{
	auto copyData(data);
	for ( auto& copy : copyData )
		copy *= factor;
	return copyData;
}


inline double
Dot (
	const array3d& data1,
	const array3d& data2
) noexcept
{
	return std::inner_product (
		data1.begin(), data1.end(),
		data2.begin(),
		0.
	);
}


inline matrix3d
OuterDot (
	const array3d& a,
	const array3d& b
) noexcept
{
	using namespace CartesianCoordinate_Namespace;

	return {{
		{a[X]*b[X], a[X]*b[Y], a[X]*b[Z]},
		{a[Y]*b[X], a[Y]*b[Y], a[Y]*b[Z]},
		{a[Z]*b[X], a[Z]*b[Y], a[Z]*b[Z]}
	}};
}


inline array3d
Cross (
	const array3d& a,
	const array3d& b
) noexcept
{
	using namespace CartesianCoordinate_Namespace;

	return {a[Y]*b[Z]-a[Z]*b[Y], a[Z]*b[X]-a[X]*b[Z], a[X]*b[Y]-a[Y]*b[X]};
}


inline double
TripleDot ( // a . ( b x c )
	const array3d& a,
	const array3d& b,
	const array3d& c
) noexcept
{
	return Dot(a, Cross(b, c));
}


inline array3d
VectorTripleDot ( // a x ( b x c )
	const array3d& a,
	const array3d& b,
	const array3d& c
) noexcept
{
	return Minus(
		Scale(b, Dot(a, c)),
		Scale(c, Dot(a, b))
	);
}


inline double
Square (
	const array3d& data
) noexcept
{
	return Dot(data, data);
}


inline array3d
operator- (
	const array3d& a
) noexcept
{
	return Scale(a, -1.);
}


inline double
Norm (
	const array3d& rij
) noexcept
{
	return sqrt( Square(rij) );
}


inline array3d
UnitVector (
	const array3d& rij,
	double inverseRij
) noexcept
{
	return Scale(rij, inverseRij);
}


inline array3d
UnitVector (
	const array3d& rij
) noexcept
{
	auto rij_ = Norm(rij);
	return UnitVector(rij, 1./rij_);
}


inline double
CosTheta (
	const array3d& nij,
	const array3d& nik
) noexcept
{
	return Dot(nij, nik);
}


inline double
CosPhi (
	const array3d& rij,
	const array3d& rjk,
	const array3d& rkl
) noexcept
{
	auto normalIJK = UnitVector(Cross(rij, rjk));
	auto normalJKL = UnitVector(Cross(rjk, rkl));
	return CosTheta(normalIJK, normalJKL);
}

//---------------------------------------------------------------------------//

inline matrix3d
Identity () noexcept
{
	return {{
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	}};
}


inline void
Increment (
	matrix3d& a,
	const matrix3d& b
) noexcept
{
	auto aIter = a.begin();
	auto bIter = b.begin();
	for(; aIter!=a.end(); ++aIter, ++bIter)
		Increment(*aIter, *bIter);
}


inline matrix3d
Minus (
	const matrix3d& a,
	const matrix3d& b
) noexcept
{
	return {{
		{a[0][0]-b[0][0], a[0][1]-b[0][1], a[0][2]-b[0][2]},
		{a[1][0]-b[1][0], a[1][1]-b[1][1], a[1][2]-b[1][2]},
		{a[2][0]-b[2][0], a[2][1]-b[2][1], a[2][2]-b[2][2]}
	}};
}


inline matrix3d
Scale (
	const matrix3d& m,
	double factor
) noexcept
{
	auto copyMatrix(m);
	for ( auto& array : copyMatrix )
		for ( auto& d : array )
			d *= factor;
	return copyMatrix;
}

//---------------------------------------------//

template < typename T, typename... Ts >
void Increment (
	T& a,
	const T& b,
	const Ts&... ts
) noexcept
{
	Increment(a, b);
	Increment(a, ts...);
}


template < typename T >
T
Addition (
	const T& a,
	const T& b
) noexcept
{
	auto result(a);
	Increment(result, b);
	return result;
}


template < typename T, typename... Ts >
T
Addition (
	const T& a,
	const T& b,
	const Ts&... ts
) noexcept
{
	auto c = Addition(a, b);
	Increment(c, ts...);
	return c;
}

#endif // DATA3D_H_INCLUDED

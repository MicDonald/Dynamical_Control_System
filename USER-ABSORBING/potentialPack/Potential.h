#ifndef POTENTIAL_H_INCLUDED
#define POTENTIAL_H_INCLUDED

/* definition:
	rij = rj - ri
	fij = force at i due to j
	Hij = potential second derivative about ri and rj, dri( drj(E) )

	pair: i-j
	angle: j-i-k
	dihedral: i-j-k-l
*/

#include "array3d.h"
#include "PotentialFlags.h"
#include <functional>


class Potential {
public:
	// pair
	virtual double
	Energy (
		const array3d& rij
	) const noexcept;


	// angle
	virtual double
	Energy (
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	// dihedral
	virtual double
	Energy (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept;


	// many body
	virtual double
	Energy (
		const array3d& rij,
		const std::vector<array3d>& rik
	) const noexcept;

	//-------------------------------------//

	// pair
	virtual array3d // fij
	Force (
		const array3d& rij
	) const noexcept;


	// pair, finite difference
	virtual array3d
	Force (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept;


	// angle
	virtual std::array<array3d, 2> // fji, fki
	Force (
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	// angle, finite difference
	virtual std::array<array3d, 2>
	Force (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept;


	// dihedral
	virtual std::array<array3d, 3> // fij, fkj, flk
	Force (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept;


	// dihedral, finite difference
	virtual std::array<array3d, 3>
	Force (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		FiniteDifference_t
	) const noexcept;


	// many body
	virtual std::vector<array3d> // fji, fki[]
	Force (
		const array3d& rij,
		const std::vector<array3d>& rik
	) const noexcept;


	// many body, finite difference
	virtual std::vector<array3d>
	Force (
		const array3d& rij,
		const std::vector<array3d>& rik,
		FiniteDifference_t
	) const noexcept;

	//-------------------------------------//

	// pair
	virtual matrix3d // Hij
	Hessian (
		const array3d& rij
	) const noexcept;


	// pair, finite difference
	virtual matrix3d
	Hessian (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept;


	// angle
	virtual std::array<matrix3d, 3> // Hij, Hik, Hjk
	Hessian (
		const array3d& rij,
		const array3d& rik
	) const noexcept;


	// angle, finite difference
	virtual std::array<matrix3d, 3>
	Hessian (
		const array3d& rij,
		const array3d& rik,
		FiniteDifference_t
	) const noexcept;


	// dihedral
	virtual std::array<matrix3d, 6> // Hij, Hik, Hil, Hjk, Hjl, Hkl
	Hessian (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl
	) const noexcept;


	// dihedral, finite difference
	virtual std::array<matrix3d, 6>
	Hessian (
		const array3d& rij,
		const array3d& rjk,
		const array3d& rkl,
		FiniteDifference_t
	) const noexcept;


	// many body
	virtual std::vector<matrix3d> // Hij, Hik[], Hjk[], Hk[]k[]
	Hessian (
		const array3d& rij,
		const std::vector<array3d>& rik
	) const noexcept;


	// many body, finite difference
	virtual std::vector<matrix3d>
	Hessian (
		const array3d& rij,
		const std::vector<array3d>& rik,
		FiniteDifference_t
	) const noexcept;

protected:
	constexpr array3d
	ZeroForce () const noexcept;


	constexpr matrix3d
	ZeroHessian () const noexcept;
};

#endif // POTENTIAL_H_INCLUDED

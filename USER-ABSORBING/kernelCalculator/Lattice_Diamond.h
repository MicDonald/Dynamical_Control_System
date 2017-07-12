#ifndef LATTICE_DIAMOND_H_INCLUDED
#define LATTICE_DAIMOND_H_INCLUDED

#include "array3d.h"

struct Lattice_Diamond {
	enum { DIM=3, N=8 };

	Lattice_Diamond (
		double spacing = 5.4309 // Si
	) noexcept;


	array3d bases[N];
	matrix3d primitive;
};

#endif // LATTICE_DIAMOND_H_INCLUDED

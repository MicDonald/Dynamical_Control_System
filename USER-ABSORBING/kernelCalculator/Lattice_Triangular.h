#ifndef LATTICE_TRIANGULAR_H_INCLUDED
#define LATTICE_TRIANGULAR_H_INCLUDED

/*
   ┌────────┐
   │        │
   │     O  │ <--  id 1
   │        │
   │        │
   │ O      │ <--  id 0
   └────────┘

   0: (0., 0.)
   1: (0.5, 0.866)
*/

#include "array3d.h"

class Lattice_Triangular {
public:
	enum {DIM=2, N=2};


	Lattice_Triangular (
		double spacing = 1.
	) noexcept;

	array3d bases[N];
	array3d primitive[DIM];
};

#endif // LATTICE_TRIANGULAR_H_INCLUDED

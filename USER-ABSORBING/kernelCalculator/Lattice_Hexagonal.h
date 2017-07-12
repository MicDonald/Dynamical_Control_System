/*
   ┌────────┐
   │ O      │ <--  id 0
   │     O  │ <--  id 1
   │        │
   │     O  │ <--  id 2
   │ O      │ <--  id 3
   └────────┘

   0: (0., 2.)
   1: (0.866, 1.5)
   2: (0.866, 0.5)
   3: (0., 0.)

   Hexagonal: spacing = 2.4612/sqrt(3) = 1.42097448252951
*/

#include <cmath>

struct Lattice_Hexagonal {
	enum {DIM=2, N=4};


	double bases[N][DIM];
	double primitive[DIM][DIM];


	Lattice_Hexagonal (
		double spacing = 1.42097448252951
	) noexcept
	{
		bases[0][0] = 0.;
		bases[0][1] = 2.*spacing;

		bases[1][0] = sqrt(3.)*0.5*spacing;
		bases[1][1] = 1.5*spacing;

		bases[2][0] = bases[1][0];
		bases[2][1] = 0.5*spacing;

		bases[3][0] = 0.;
		bases[3][1] = 0.;


		primitive[0][0] = bases[1][0] * 2.;
		primitive[0][1] = 0.;

		primitive[1][0] = 0.;
		primitive[1][1] = 3.*spacing;
	}
};


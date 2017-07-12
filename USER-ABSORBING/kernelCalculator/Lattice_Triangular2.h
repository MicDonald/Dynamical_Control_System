/*
   ┌────────┐
   │        │
   │ O      │ <--  id 1
   │        │
   │        │
   │     O  │ <--  id 0
   └────────┘

   0: (0.5, 0.)
   1: (0, 0.866)
*/

#include <cmath>

struct Lattice_Triangular2 {
	enum {DIM=2, N=2};


	double bases[N][DIM];
	double primitive[DIM][DIM];


	Lattice_Triangular2 (
		double spacing = 1.
	) noexcept
	{
		auto ha = 0.5*spacing;
		bases[0][0] = ha;
		bases[0][1] = 0.;

		bases[1][0] = 0.;
		bases[1][1] = ha*sqrt(3.);


		primitive[0][0] = spacing;
		primitive[0][1] = 0.;

		primitive[1][0] = 0.;
		primitive[1][1] = sqrt(3.)*spacing;
	}
};

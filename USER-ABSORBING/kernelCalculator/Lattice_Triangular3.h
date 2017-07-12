/*
   ┌──────────┐
   │          │
   │       O  │ <--  id 1
   │          │
   │ O        │ <--  id 0
   └──────────┘

   0: (0., 0.)
   1: (0.866, 0.5)
*/

#include <cmath>

struct Lattice_Triangular3 {
	enum {DIM=2, N=2};


	double bases[N][DIM];
	double primitive[DIM][DIM];


	Lattice_Triangular3 (
		double spacing = 1.
	) noexcept
	{
		auto ha = 0.5*spacing;
		bases[0][0] = 0.;
		bases[0][1] = 0.;

		bases[1][0] = ha*sqrt(3.);
		bases[1][1] = ha;


		primitive[0][0] = sqrt(3.)*spacing;
		primitive[0][1] = 0.;

		primitive[1][0] = 0.;
		primitive[1][1] = spacing;
	}
};

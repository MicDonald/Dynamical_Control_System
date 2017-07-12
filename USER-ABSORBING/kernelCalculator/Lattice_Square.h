/*
    
    |
    O-
   0: (0., 0.)
   
*/

#include <cmath>

struct Lattice_Square {
	enum {DIM=2, N=1};

	double bases[N][DIM];
	double primitive[DIM][DIM];


	Lattice_Square (
		double spacing = 1.
	) noexcept
	{		
		bases[0][0] = 0.;
		bases[0][1] = 0.;
		
		primitive[0][0] = spacing;
		primitive[0][1] = 0.;
		
		primitive[1][0] = 0.;
		primitive[1][1] = spacing;
		
	}
};

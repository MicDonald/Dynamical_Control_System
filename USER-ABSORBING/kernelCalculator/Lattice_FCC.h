/*
   ┌────────┐
   │        │                y 
   │ O   o  │ <--  id 3 2    |
   │        │                |
   │        │                z --- x
   │ o   O  │ <--  id 0 1
   └────────┘

   0: (0., 0., 0.)
   1: (0.5, 0., 0.5)
   2: (0.5, 0.5, 0.)
   3: (0., 0.5, 0.5)
*/

struct Lattice_FCC {
	enum {DIM=3, N=4};


	double bases[N][DIM];
	double primitive[DIM][DIM];


	Lattice_FCC (
		double spacing = 1.
	) noexcept
	{
		double ha = 0.5 * spacing;

		bases[1][0] = ha;
		bases[1][1] = 0.;
		bases[1][2] = ha;

		bases[2][0] = ha;
		bases[2][1] = ha;
		bases[2][2] = 0.;

		bases[3][0] = 0.;
		bases[3][1] = ha;
		bases[3][2] = ha;


		primitive[0][0] = spacing;
		primitive[0][1] = 0.;
		primitive[0][2] = 0.;

		primitive[1][0] = 0.;
		primitive[1][1] = spacing;
		primitive[1][2] = 0.;

		primitive[2][0] = 0.;
		primitive[2][1] = 0.;
		primitive[2][2] = spacing;
	}
};


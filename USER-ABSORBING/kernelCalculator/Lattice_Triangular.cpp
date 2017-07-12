#include "Lattice_Triangular.h"
#include <cmath>

Lattice_Triangular::Lattice_Triangular (
	double spacing
) noexcept
{
	auto ha = 0.5*spacing;
	bases[0] = {0., 0.};
	bases[1] = {ha, ha*sqrt(3.)};

	primitive[0] = {spacing, 0.};
	primitive[1] = {0., sqrt(3.)*spacing};
}


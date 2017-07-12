#include "Lattice_Diamond.h"

Lattice_Diamond::Lattice_Diamond (
	double spacing
) noexcept
{
	auto ha = 0.5 * spacing;
	auto qa = 0.25 * spacing;
	auto qa3 = 0.75 * spacing;

	bases[0] = {0., 0., 0.};
	bases[1] = {ha, ha, 0.};
	bases[2] = {ha, 0., ha};
	bases[3] = {0., ha, ha};
	bases[4] = {qa, qa, qa};
	bases[5] = {qa3, qa3, qa};
	bases[6] = {qa, qa3, qa3};
	bases[7] = {qa3, qa, qa3};

	primitive[0] = {spacing, 0., 0.};
	primitive[1] = {0., spacing, 0.};
	primitive[2] = {0., 0., spacing};
}


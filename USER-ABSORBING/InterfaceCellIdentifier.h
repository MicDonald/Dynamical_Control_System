#ifndef INTERFACECELLIDENTIFIER_H_INCLUDED
#define INTERFACECELLIDENTIFIER_H_INCLUDED

#include "InterfaceCellIdentifier1D.h"
#include "InterfaceCellIdentifier2D.h"
#include "InterfaceCellIdentifier3D.h"
#include <type_traits>


template < typename LATTICE >
using InterfaceCellIdentifier = std::conditional_t <
	LATTICE::DIM == 1,
	InterfaceCellIdentifier1D< LATTICE >,
	std::conditional_t <
		LATTICE::DIM == 2,
		InterfaceCellIdentifier2D< LATTICE >,
		InterfaceCellIdentifier3D< LATTICE >
	>
>;

#endif // INTERFACECELLIDENTIFIER_H_INCLUDED

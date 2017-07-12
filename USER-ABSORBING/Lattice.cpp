#include "Lattice.h"

constexpr std::array<data3d, 2> Lattice_Triangular::primitive;
constexpr std::array<data3d, 2> Lattice_Triangular::shift;
constexpr std::array<data3d, 2> Lattice_Triangular::basis;

//---------------------------------------------//

constexpr std::array<data3d, 2> Lattice_Hexagonal::primitive;
constexpr std::array<data3d, 2> Lattice_Hexagonal::shift;
constexpr std::array<data3d, 4> Lattice_Hexagonal::basis;

//---------------------------------------------------------------------------//

constexpr std::array<data3d, 3> Lattice_FCC::primitive;
constexpr std::array<data3d, 3> Lattice_FCC::shift;
constexpr std::array<data3d, 4> Lattice_FCC::basis;

//---------------------------------------------------------------------------//

constexpr std::array<data3d, 3> Lattice_Diamond::primitive;
constexpr std::array<data3d, 3> Lattice_Diamond::shift;
constexpr std::array<data3d, 8> Lattice_Diamond::basis;

//---------------------------------------------------------------------------//

constexpr std::array<data3d, 2> Lattice_Square::primitive;
constexpr std::array<data3d, 2> Lattice_Square::shift;
constexpr std::array<data3d, 1> Lattice_Square::basis;



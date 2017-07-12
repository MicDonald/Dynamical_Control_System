#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED

#include "data3d.h"
#include <cmath>
#include <vector>
#include <array>
#include <string>


enum LatticeType {TRIANGULAR, HEXAGONAL, FCC, DIAMOND, Si, SQUARE,NOTYPE};

auto String2LatticeType = []( const std::string& type )
	{
		if ( type.compare("Triangular") == 0 ) return TRIANGULAR;
		else if ( type.compare("Hexagonal") == 0 ) return HEXAGONAL;
		else if ( type.compare("FCC") == 0 ) return FCC;
		else if ( type.compare("Square") == 0 ) return SQUARE;
		throw 0;
		return NOTYPE;
	};

//---------------------------------------------------------------------------//

struct ShiftLattice {
	/*template < unsigned DIM, template T >
	T
	operator() (
		const std::array<data3d, DIM>& primitive,
		const std::array<data3d, DIM>& shift,
		T basis
	) const noexcept
	{
		for ( auto& iBasis : basis )
		{
			iBasis[0] += shift[0][0];
			if ( iBasis[0] > primitive[0][0] )
				iBasis[0] -= primitive[0][0];
		}
		return basis;
	}*/


	template < unsigned long DIM, unsigned long N >
	std::array<data3d, N>
	operator() (
		const std::array<data3d, DIM>& primitive,
		const std::array<data3d, DIM>& shift,
		const std::array<data3d, N>& basis
	) const noexcept
	{
		auto shiftBasis(basis);
		for ( auto& iBasis : shiftBasis )
		{
			iBasis[0] += shift[0][0];
			if ( iBasis[0] >= primitive[0][0] )
				iBasis[0] -= primitive[0][0];
		}
		return shiftBasis;
	}


	template < unsigned long DIM >
	std::vector<data3d>
	operator() (
		const std::array<data3d, DIM>& primitive,
		const std::array<data3d, DIM>& shift,
		const std::vector<data3d>& basis
	) const noexcept
	{
		auto shiftBasis(basis);
		for ( auto& iBasis : shiftBasis )
		{
			iBasis[0] += shift[0][0];
			if ( iBasis[0] >= primitive[0][0] )
				iBasis[0] -= primitive[0][0];
		}
		return shiftBasis;
	}
};

//---------------------------------------------//

struct Lattice_Triangular {
	static constexpr auto kernelFile = "_2D_Tri";
	static constexpr unsigned DIM = 2;
	static constexpr unsigned N = 2;
	static constexpr std::array<data3d, DIM> primitive = {{
		{1., 0.},
		{0., sqrt(3.)}
	}};
	static constexpr std::array<data3d, DIM> shift = {{
		{0.5, 0.},
		{0., 0.5*sqrt(3.)}
	}};
	static constexpr std::array<data3d, N> basis = {{
		{0., 0.},
		{0.5, 0.5*sqrt(3.)}
	}};


	static bool
	Rotate ( // return orientation shifted after rotate
		unsigned& localID,
		char axis,
		int angle
	) noexcept
	{
		if ( axis == 'z' && angle == 180 )
			localID = 1 - localID;
		return false;
	}
};

//---------------------------------------------//

struct Lattice_Hexagonal {
	static constexpr auto kernelFile = "_2D_Hex";
	static constexpr unsigned DIM = 2;
	static constexpr unsigned N = 4;
	static constexpr std::array<data3d, DIM> primitive = {{
		{1., 0.},
		{0., sqrt(3.)}
	}};
	static constexpr std::array<data3d, DIM> shift = {{
		{0.5, 0.},
		{0., 0.5*sqrt(3.)}
	}};
	static constexpr std::array<data3d, N> basis = {{
		{0., 0.},
		{0.5, 0.5/sqrt(3.)},
		{0.5, 1.5/sqrt(3.)},
		{0., 2./sqrt(3.)}
	}};


	static bool
	Rotate ( // return orientation shifted after rotate
		unsigned& localID,
		char axis,
		int angle
	) noexcept
	{
		if ( axis == 'z' && angle == 180 )
			localID = 3 - localID;
		return true;
	}
};

using Lattice_Honeycomb = Lattice_Hexagonal;

//---------------------------------------------------------------------------//

struct Lattice_FCC {
	static constexpr auto kernelFile = "_3D_FCC";
	static constexpr unsigned DIM = 3;
	static constexpr unsigned N = 4;
	static constexpr std::array<data3d, DIM> primitive = {{
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	}};
	static constexpr std::array<data3d, DIM> shift = {{
		{0.5, 0., 0.},
		{0., 0.5, 0.},
		{0., 0., 0.5}
	}};
	static constexpr std::array<data3d, N> basis = {{
		{0., 0., 0.},
		{0., 0.5, 0.5},
		{0.5, 0., 0.5},
		{0.5, 0.5, 0.} 
	}};


	static bool
	Rotate (
		unsigned& localID,
		char axis,
		int angle
	) noexcept
	{
		if ( angle == 180 )
		{
			if ( axis == 'x' )
				localID = localID < 2 ? 1-localID : 5-localID;
			else if ( axis == 'y' )
				localID = localID%2==0 ? 2-localID : 4-localID;
			else if ( axis == 'z' )
				localID = 3 - localID;
			return false;
		}
		else if ( angle == 90 )
		{
			if ( axis == 'x' )
				localID = localID<2 ? 3-localID : localID-2;
			else if ( axis == 'y' )
				localID = localID%2==0 ? 2-localID : 4-localID;
			else if ( axis == 'z' )
			{
				if ( localID>0 && localID<3 )
					localID = 3 - localID;
			}
			return true;
		}
		return false;
	}
};

struct Lattice_Ar : public Lattice_FCC {
	static constexpr auto kernelFile = "_3D_Ar";
};	

//---------------------------------------------//

struct Lattice_Diamond {
	static constexpr unsigned DIM = 3;
	static constexpr unsigned N = 8;
	static constexpr std::array<data3d, DIM> primitive = {{
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	}};
	static constexpr std::array<data3d, DIM> shift = {{
		{0.25, 0., 0.},
		{0., 0.25, 0.},
		{0., 0., 0.25}
	}};
	static constexpr std::array<data3d, N> basis = {{
		{0., 0., 0.},
		{0.5, 0.5, 0.},
		{0.5, 0., 0.5},
		{0., 0.5, 0.5},
		{0.25, 0.25, 0.25},
		{0.75, 0.75, 0.25},
		{0.25, 0.75, 0.75},
		{0.75, 0.25, 0.75}
	}};


	static bool
	Rotate (
		unsigned& localID,
		char axis,
		int angle
	) noexcept
	{
/*		if ( axis == 'x' && angle == 180 )
			switch ( localID )
			{
			default: ;
			}
		else if ( axis == 'y' && angle == 180 )
			switch ( localID )
			{
			default: ;
			}
		else if ( axis == 'z' && angle == 180 )
			;*/
		return false;
	}
};


struct Lattice_Si : public Lattice_Diamond {
	static constexpr auto kernelFile = "_3D_Si";
};

struct Lattice_Square {
        static constexpr auto kernelFile = "_2D_Squ";
        static constexpr unsigned DIM = 2;
        static constexpr unsigned N = 1;
        static constexpr std::array<data3d, DIM> primitive = {{
                {1., 0.},
                {0., 1.}
        }};
        static constexpr std::array<data3d, DIM> shift = {{
                {0.25, 0.},
                {0., 0.25}
        }};
        static constexpr std::array<data3d, N> basis = {{
                {0.5, 0.5},
        }};


        static bool
        Rotate ( // return orientation shifted after rotate
                unsigned& localID,
                char axis,
                int angle
        ) noexcept
        {
		localID = localID;
        	return true;
              
        }
};



#endif // LATTICE_H_INCLUDED

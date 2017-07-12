/*
  <---  width  --->
  ┌───┬───┬───┬───┐
  │ m │ n │ o │ p │
  ├───┼───┼───┼───┤
  │ i │ j │ k │ l │
  ├───┼───┼───┼───┤
  │ e │ f │ g │ h │
  ├───┼───┼───┼───┤
  │ a │ b │ c │ d │
  └───┴───┴───┴───┘

   interface 
  ┌───┐¦┌───┬───┬───┬───┐
  │   │|│ M │ N │ O │ P │          y
  ├───┤|├───┼───┼───┼───┤          ↑
  │ w │|│ I │ J │ K │ L │          │
  ├───┤|├───┼───┼───┼───┤         z└───→ x
  │ v │|│ E │ F │ G │ H │
  ├───┤|├───┼───┼───┼───┤ 
  │   │|│ A │ B │ C │ D │
  │ u │|└───┴───┴───┴───┘
  │   │|----------------- 
  ├───┼─────┬───┬───┬───┐
  │ q │  r  │ s │ t │   │
  └───┴─────┴───┴───┴───┘

  F = -KU - K'U'
  U = [a b c d ...]^T	(alone x->y) -> alone z -> next layer
  U' = [A B C D ...]^T	alone y -> alone z
    ┌                           ┐
    │  Kaa  Kba  Kca  Kda  ...  │
    │  Kab  Kbb  Kcb  Kdb  ...  │
  K=│  Kac  Kbc  Kcc  Kdc  ...  │
    │  Kad  Kbd  Kcd  Kdd  ...  │
    │   ¦    ¦    ¦    ¦        │
    └                           ┘
*/
#include <vector>


struct Topology_3D_Corner {
	MKL_INT nUnits;
	MKL_INT width;
	MKL_INT planeUnits;
	MKL_INT lineUnits;
	enum {DIM = 3};


	static constexpr MKL_INT
	numOfUnits_Of (
		MKL_INT width
	) noexcept
	{
		return 3*width*(width-1)+1;
	}


	Topology_3D_Corner (
		MKL_INT nUnits
	) noexcept :
		nUnits( nUnits )
	{
		width = (3.+sqrt(12*nUnits-3.))/6. + 0.5;
		planeUnits = width*width;
		lineUnits = width*2 - 1;
	}


	template < typename STIFF >
	void
	checkDIM () const noexcept
	{
		static_assert(unsigned(DIM)==unsigned(STIFF::DIM), "diff DIM");
	}


	template < typename STIFF >
	decltype(auto) 
	operator() (
		MKL_INT& i,
		MKL_INT& j,
		const STIFF&
	) const noexcept
	{
		checkDIM<STIFF>();

		constexpr MKL_INT nDOFPerUnit = STIFF::DOF;
		auto nDOFPerLayer = nDOFPerUnit * nUnits;
		auto nDOFPerPlane = nDOFPerUnit * planeUnits;
		auto nDOFPerLine = nDOFPerUnit * lineUnits;

		auto iLayer = i/nDOFPerLayer; // x
		i %= nDOFPerLayer;
		MKL_INT locationI[DIM] = {};
		if( i<nDOFPerPlane )
		{
			auto ii = i/nDOFPerUnit;
			locationI[1] = ii/width;
			locationI[0] = ii%width;
			i %= nDOFPerUnit;
		}
		else
		{
			i -= nDOFPerPlane;
			locationI[2] = i/nDOFPerLine + 1;
			i %= nDOFPerLine;
			auto iUnit = i/nDOFPerUnit;
			locationI[1] = iUnit >= width ? iUnit-width+1 : 0;
			locationI[0] = iUnit >= width ? 0 : iUnit;
			i %= nDOFPerUnit;
		}
		locationI[0] -= iLayer;
		locationI[1] -= iLayer;
		locationI[2] -= iLayer;

		auto jLayer = j/nDOFPerLayer;
		j %= nDOFPerLayer;
		MKL_INT locationJ[DIM] = {};
		if( j<nDOFPerPlane )
		{
			auto jj = j/nDOFPerUnit;
			locationJ[1] = jj/width;
			locationJ[0] = jj%width;
			j %= nDOFPerUnit;
		}
		else
		{
			j -= nDOFPerPlane;
			locationJ[2] = j/nDOFPerLine + 1;
			j %= nDOFPerLine;
			auto jUnit = j/nDOFPerUnit;
			locationJ[1] = jUnit >= width ? jUnit-width+1 : 0;
			locationJ[0] = jUnit >= width ? 0 : jUnit;
			j %= nDOFPerUnit;
		}
		locationJ[0] -= jLayer;
		locationJ[1] -= jLayer;
		locationJ[2] -= jLayer;

		if(locationI[2]==locationJ[2])
		{
			if(locationI[0]==locationJ[0])
			{
				if(locationI[1]==locationJ[1]) return &STIFF::_000;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::_010;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::_0_10;
			}
			else if(locationI[0]==locationJ[0]+1)
			{
				if(locationI[1]==locationJ[1]) return &STIFF::_100;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::_110;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::_1_10;
			}
			else if(locationI[0]==locationJ[0]-1)
			{
				if(locationI[1]==locationJ[1]) return &STIFF::__100;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::__110;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::__1_10;
			}
		}
		else if(locationI[2]==locationJ[2]+1)
		{
			if(locationI[0]==locationJ[0])
			{
				if(locationI[1]==locationJ[1]) return &STIFF::_001;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::_011;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::_0_11;
			}
			else if(locationI[0]==locationJ[0]+1)
			{
				if(locationI[1]==locationJ[1]) return &STIFF::_101;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::_111;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::_1_11;
			}
			else if(locationI[0]==locationJ[0]-1)
			{
				if(locationI[1]==locationJ[1]) return &STIFF::__101;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::__111;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::__1_11;
			}
		}
		else if(locationI[2]==locationJ[2]-1)
		{
			if(locationI[0]==locationJ[0])
			{
				if(locationI[1]==locationJ[1]) return &STIFF::_00_1;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::_01_1;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::_0_1_1;
			}
			else if(locationI[0]==locationJ[0]+1)
			{
				if(locationI[1]==locationJ[1]) return &STIFF::_10_1;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::_11_1;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::_1_1_1;
			}
			else if(locationI[0]==locationJ[0]-1)
			{
				if(locationI[1]==locationJ[1]) return &STIFF::__10_1;
				else if(locationI[1]==locationJ[1]+1) return &STIFF::__11_1;
				else if(locationI[1]==locationJ[1]-1) return &STIFF::__1_1_1;
			}
		}
		return &STIFF::zero;
	}


	std::vector<MKL_INT>
	range (	MKL_INT width,
		MKL_INT dofPerUnit
	) const noexcept
	{
		width = width < this->width-1 ? width : this->width-1;
		auto startP = 0;
		auto endP = (width+1) * dofPerUnit;
		auto startL1 = 0;
		auto endL1 = width * dofPerUnit;
		auto startL2 = this->width*dofPerUnit;
		auto endL2 = startL2 + width*dofPerUnit;

		auto dofPerWidth = this->width * dofPerUnit;
		auto dofPerLine = lineUnits * dofPerUnit;
		auto dofPerPlane = planeUnits * dofPerUnit;

		std::vector<MKL_INT> ids;
		for(MKL_INT j=0; j<=width; ++j)
			for(MKL_INT i=startP; i<endP; ++i)
				ids.push_back(i+j*dofPerWidth);

		for(MKL_INT j=0; j<width; ++j)
		{
			for(MKL_INT i=startL1; i<endL1; ++i)
				ids.push_back(dofPerPlane + j*dofPerLine + i);
			for(MKL_INT i=startL2; i<endL2; ++i)
				ids.push_back(dofPerPlane + j*dofPerLine + i);
		}
		return ids;
	}
};


/*
               interface 
  ┌───┬───┬───┬───┐¦┌───┐
  │   │   │   │ d │¦│ D │          y
  ├───┼───┼───┼───┤¦├───┤          ↑
  │   │   │   │ c │¦│ C │          │
  ├───┼───┼───┼───┤¦├───┤         z└───→ x
  │   │   │   │ b │¦│ B │
  ├───┼───┼───┼───┤¦├───┤ 
  │   │   │   │ a │¦│ A │
  └───┴───┴───┴───┘¦└───┘

  F = -KU - K'U'
  U = [a b c d ...]^T	alone y -> alone z -> alone -x
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


struct Topology_3D_Face {
	MKL_INT nUnits;
	MKL_INT yUnits;
	MKL_INT zUnits;
	enum {DIM = 3};


	static constexpr unsigned
	numOfUnits_Of (
		MKL_INT width
	) noexcept
	{
		return width*width;
	}


	Topology_3D_Face (
		MKL_INT nUnits
	) noexcept :
		nUnits( nUnits )
	{
		yUnits = sqrt(nUnits);
		zUnits = nUnits / yUnits;
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

		auto nDOFPerUnit = STIFF::DOF;
		auto nDOFPerSide = nDOFPerUnit * yUnits;
		auto nDOFPerLayer = nDOFPerUnit * nUnits;

		auto iLayer = i/nDOFPerLayer; // x
		i %= nDOFPerLayer;
		auto izUnit = i/nDOFPerSide;
		i %= nDOFPerSide;
		auto iyUnit = i/nDOFPerUnit;
		i %= nDOFPerUnit;

		auto jLayer = j/nDOFPerLayer;
		j %= nDOFPerLayer;
		auto jzUnit = j/nDOFPerSide;
		j %= nDOFPerSide;
		auto jyUnit = j/nDOFPerUnit;
		j %= nDOFPerUnit;
		if(iLayer == jLayer)
		{
			if (izUnit==jzUnit && iyUnit==jyUnit) return &STIFF::_000;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit) return &STIFF::_00_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit) return &STIFF::_001;

			else if(izUnit==jzUnit && iyUnit==jyUnit-1) return &STIFF::_0_10;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit-1) return &STIFF::_0_1_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit-1) return &STIFF::_0_11;

			else if(izUnit==jzUnit && iyUnit==jyUnit+1) return &STIFF::_010;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit+1) return &STIFF::_01_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit+1) return &STIFF::_011;
		}
		else if(iLayer == jLayer+1)
		{
			if(izUnit==jzUnit && iyUnit==jyUnit) return &STIFF::__100;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit) return &STIFF::__10_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit) return &STIFF::__101;

			else if(izUnit==jzUnit && iyUnit==jyUnit-1) return &STIFF::__1_10;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit-1) return &STIFF::__1_1_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit-1) return &STIFF::__1_11;

			else if(izUnit==jzUnit && iyUnit==jyUnit+1) return &STIFF::__110;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit+1) return &STIFF::__11_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit+1) return &STIFF::__111;
		}
		else if(iLayer == jLayer-1)
		{
			if(izUnit==jzUnit && iyUnit==jyUnit) return &STIFF::_100;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit) return &STIFF::_10_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit) return &STIFF::_101;

			else if(izUnit==jzUnit && iyUnit==jyUnit-1) return &STIFF::_1_10;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit-1) return &STIFF::_1_1_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit-1) return &STIFF::_1_11;

			else if(izUnit==jzUnit && iyUnit==jyUnit+1) return &STIFF::_110;
			else if(izUnit==jzUnit-1 && iyUnit==jyUnit+1) return &STIFF::_11_1;
			else if(izUnit==jzUnit+1 && iyUnit==jyUnit+1) return &STIFF::_111;
		}
		return &STIFF::zero;
	}


	std::vector<MKL_INT>
	range (	MKL_INT width,
		MKL_INT dofPerUnit
	) const noexcept
	{
		auto yWidth = width < yUnits/2 ? width : yUnits/2;
		auto zWidth = width < zUnits/2 ? width : zUnits/2;
		auto yStart = (yUnits/2 - yWidth) * dofPerUnit;
		auto yEnd = (yUnits/2 + yWidth + 1) * dofPerUnit;
		auto zStart = zUnits/2 - zWidth;
		auto zEnd = zUnits/2 + zWidth + 1;
//std::cout<<yStart<<" "<<yEnd<<" "<<zStart<<" "<<zEnd<<std::endl;
		std::vector<MKL_INT> ids;
		for(MKL_INT j=zStart; j<zEnd; ++j)
			for(MKL_INT i=yStart; i<yEnd; ++i)
				ids.push_back( j*yUnits*dofPerUnit+i );

		return ids;
	}
};


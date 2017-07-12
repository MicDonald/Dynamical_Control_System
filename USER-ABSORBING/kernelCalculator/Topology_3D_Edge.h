/*
   interface 
  ┌───┐¦┌───┬───┬───┬───┐
  │   │|│ A │   │   │   │          y
  ├───┤|├───┼───┼───┼───┤          ↑
  │ a │|│ B │   │   │   │          │
  ├───┤|├───┼───┼───┼───┤         z└───→ x
  │ b │|│ C │   │   │   │
  ├───┤|├───┼───┼───┼───┤ 
  │   │|│ D │ E │ F │ G │
  │ c │|└───┴───┴───┴───┘
  │   │|----------------- 
  ├───┼─────┬───┬───┬───┐
  │ d │  e  │ f │ g │   │
  └───┴─────┴───┴───┴───┘

  F = -KU - K'U'
  U = [a b c d ...]^T	(alone -y->x) -> alone z -> next layer
  U' = [A B C D ...]^T	alone -y->x -> alone z
    ┌                           ┐
    │  Kaa  Kba  Kca  Kda  ...  │
    │  Kab  Kbb  Kcb  Kdb  ...  │
  K=│  Kac  Kbc  Kcc  Kdc  ...  │
    │  Kad  Kbd  Kcd  Kdd  ...  │
    │   ¦    ¦    ¦    ¦        │
    └                           ┘
*/
#include <vector>


struct Topology_3D_Edge {
	MKL_INT nUnits;
	MKL_INT lineUnits;
	MKL_INT halfLineUnits;
	MKL_INT zUnits;
	enum {DIM = 3};


	static constexpr MKL_INT
	numOfUnits_Of (
		MKL_INT width
	) noexcept
	{
		return width*width;
	}


	Topology_3D_Edge (
		MKL_INT nUnits
	) noexcept :
		nUnits( nUnits )
	{
		zUnits = sqrt(nUnits);
		lineUnits = nUnits/zUnits;
		halfLineUnits = lineUnits/2;
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
		auto nDOFPerLine = nDOFPerUnit * lineUnits;
		auto nDOFPerLayer = nDOFPerUnit * nUnits;

		auto iLayer = i/nDOFPerLayer; // x
		i %= nDOFPerLayer;
		auto izUnit = i/nDOFPerLine;
		i %= nDOFPerLine;
		auto iUnit = i/nDOFPerUnit;
		i %= nDOFPerUnit;

		auto jLayer = j/nDOFPerLayer;
		j %= nDOFPerLayer;
		auto jzUnit = j/nDOFPerLine;
		j %= nDOFPerLine;
		auto jUnit = j/nDOFPerUnit;
		j %= nDOFPerUnit;

		if(iLayer == jLayer)
		{
			if(izUnit==jzUnit)
			{
				if(iUnit==jUnit) return &STIFF::_000;
				else if(iUnit==jUnit+1)
					if(iUnit<=halfLineUnits) return &STIFF::_0_10;
					else return &STIFF::_100;
				else if(iUnit==jUnit-1)
					if(iUnit<halfLineUnits) return &STIFF::_010;
					else return &STIFF::__100;
				else if(iUnit==jUnit+2 && iUnit==halfLineUnits+1) return &STIFF::_1_10;
				else if(iUnit==jUnit-2 && iUnit==halfLineUnits-1) return &STIFF::__110;
			}
			else if(izUnit==jzUnit+1)
			{
				if(iUnit==jUnit) return &STIFF::_001;
				else if(iUnit==jUnit+1)
					if(iUnit<=halfLineUnits) return &STIFF::_0_11;
					else return &STIFF::_101;
				else if(iUnit==jUnit-1)
					if(iUnit<halfLineUnits) return &STIFF::_011;
					else return &STIFF::__101;
				else if(iUnit==jUnit+2 && iUnit==halfLineUnits+1) return &STIFF::_1_11;
				else if(iUnit==jUnit-2 && iUnit==halfLineUnits-1) return &STIFF::__111;
			}
			else if(izUnit==jzUnit-1)
			{
				if(iUnit==jUnit) return &STIFF::_00_1;
				else if(iUnit==jUnit+1)
					if(iUnit<=halfLineUnits) return &STIFF::_0_1_1;
					else return &STIFF::_10_1;
				else if(iUnit==jUnit-1)
					if(iUnit<halfLineUnits) return &STIFF::_01_1;
					else return &STIFF::__10_1;
				else if(iUnit==jUnit+2 && iUnit==halfLineUnits+1) return &STIFF::_1_1_1;
				else if(iUnit==jUnit-2 && iUnit==halfLineUnits-1) return &STIFF::__11_1;
			}
		}
		else if(iLayer == jLayer+1)
		{
			if(izUnit==jzUnit)
			{
				if(iUnit==jUnit) return &STIFF::__1_10;
				else if(iUnit==jUnit+1 && jUnit>=halfLineUnits) return &STIFF::_0_10;
				else if(iUnit==jUnit+2 && jUnit>=halfLineUnits) return &STIFF::_1_10;
				else if(iUnit==jUnit-1 && jUnit<=halfLineUnits) return &STIFF::__100;
				else if(iUnit==jUnit-2 && jUnit<=halfLineUnits) return &STIFF::__110;
			}
			else if(izUnit==jzUnit+1)
			{
				if(iUnit==jUnit) return &STIFF::__1_11;
				else if(iUnit==jUnit+1 && jUnit>=halfLineUnits) return &STIFF::_0_11;
				else if(iUnit==jUnit+2 && jUnit>=halfLineUnits) return &STIFF::_1_11;
				else if(iUnit==jUnit-1 && jUnit<=halfLineUnits) return &STIFF::__101;
				else if(iUnit==jUnit-2 && jUnit<=halfLineUnits) return &STIFF::__111;
			}
			else if(izUnit==jzUnit-1)
			{
				if(iUnit==jUnit) return &STIFF::__1_1_1;
				else if(iUnit==jUnit+1 && jUnit>=halfLineUnits) return &STIFF::_0_1_1;
				else if(iUnit==jUnit+2 && jUnit>=halfLineUnits) return &STIFF::_1_1_1;
				else if(iUnit==jUnit-1 && jUnit<=halfLineUnits) return &STIFF::__10_1;
				else if(iUnit==jUnit-2 && jUnit<=halfLineUnits) return &STIFF::__11_1;
			}
		}
		else if(iLayer == jLayer-1)
		{
			if(izUnit==jzUnit)
			{
				if(iUnit==jUnit) return &STIFF::_110;
				else if(iUnit==jUnit+1 && iUnit<=halfLineUnits) return &STIFF::_100;
				else if(iUnit==jUnit+2 && iUnit<=halfLineUnits) return &STIFF::_1_10;
				else if(iUnit==jUnit-1 && iUnit>=halfLineUnits) return &STIFF::_010;
				else if(iUnit==jUnit-2 && iUnit>=halfLineUnits) return &STIFF::__110;
			}
			else if(izUnit==jzUnit+1)
			{
				if(iUnit==jUnit) return &STIFF::__111;
				else if(iUnit==jUnit+1 && iUnit<=halfLineUnits) return &STIFF::_101;
				else if(iUnit==jUnit+2 && iUnit<=halfLineUnits) return &STIFF::_1_11;
				else if(iUnit==jUnit-1 && iUnit>=halfLineUnits) return &STIFF::_011;
				else if(iUnit==jUnit-2 && iUnit>=halfLineUnits) return &STIFF::__111;
			}
			else if(izUnit==jzUnit-1)
			{
				if(iUnit==jUnit) return &STIFF::__11_1;
				else if(iUnit==jUnit+1 && iUnit<=halfLineUnits) return &STIFF::_10_1;
				else if(iUnit==jUnit+2 && iUnit<=halfLineUnits) return &STIFF::_1_1_1;
				else if(iUnit==jUnit-1 && iUnit>=halfLineUnits) return &STIFF::_01_1;
				else if(iUnit==jUnit-2 && iUnit>=halfLineUnits) return &STIFF::__11_1;
			}
		}
		return &STIFF::zero;
	}


	std::vector<MKL_INT>
	range (	MKL_INT width,
		MKL_INT dofPerUnit
	) const noexcept
	{
		auto lineWidth = width < halfLineUnits ? width : halfLineUnits;
		auto zWidth = width < zUnits/2 ? width : zUnits/2;
		auto lineStart = (halfLineUnits - lineWidth) * dofPerUnit;
		auto lineEnd = (halfLineUnits + lineWidth + 1) * dofPerUnit;
		auto zStart = zUnits/2 - zWidth;
		auto zEnd = zUnits/2 + zWidth + 1;

		std::vector<MKL_INT> ids;
		for(MKL_INT j=zStart; j<zEnd; ++j)
			for(MKL_INT i=lineStart; i<lineEnd; ++i)
				ids.push_back( j*lineUnits*dofPerUnit + i );
		return ids;
	}
};


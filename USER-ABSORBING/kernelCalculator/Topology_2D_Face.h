/*
               interface 
  ┌───┬───┬───┬───┐¦┌───┐
  │ m │ i │ e │ a │¦│ A │
  ├───┼───┼───┼───┤¦├───┤          y
  │ n │ j │ f │ b │¦│ B │          │
  ├───┼───┼───┼───┤¦├───┤          └─── x
  │ o │ k │ g │ c │¦│ C │<-half
  ├───┼───┼───┼───┤¦├───┤ 
  │ p │ l │ h │ d │¦│ D │
  └───┴───┴───┴───┘¦└───┘

  F = -KU - K'U'
  U = [a b c d e f g ...]^T   alone -y, alone -x
  U' = [A B C D ...]^T
    ┌                           ┐
    │  Kaa  Kba  Kca  Kda  ...  │
    │  Kab  Kbb  Kcb  Kdb  ...  │
  K=│  Kac  Kbc  Kcc  Kdc  ...  │
    │  Kad  Kbd  Kcd  Kdd  ...  │
    │   ¦    ¦    ¦    ¦        │
    └                           ┘
     ┌                         ┐
     │ KAa  KBa  KCa  KDa  ... │
     │ KAb  KBb  KCb  KDb  ... │
  K'=│ KAc  KBc  KCc  KDc  ... │
     │ KAd  KBd  KCd  KDd  ... │
     │  ¦    ¦    ¦    ¦       │
     └                         ┘
*/
#include <vector>

struct Topology_2D_Face {
	MKL_INT nUnits;
	enum {DIM = 2};


	static constexpr MKL_INT
	numOfUnits_Of (
		MKL_INT width
	) noexcept
	{
		return width;
	}


	Topology_2D_Face (
		MKL_INT nUnits
	) noexcept :
		nUnits( nUnits )
	{}


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
		auto nDOFPerLayer = nDOFPerUnit * nUnits;
		auto iLayer = i/nDOFPerLayer;
		i %= nDOFPerLayer;
		auto iUnit = i/nDOFPerUnit;
		i %= nDOFPerUnit;
		auto jLayer = j/nDOFPerLayer;
		j %= nDOFPerLayer;
		auto jUnit = j/nDOFPerUnit;
		j %= nDOFPerUnit;
		if(iLayer == jLayer)
		{
			if(iUnit==jUnit) return &STIFF::_00;
			else if(iUnit==jUnit+1) return &STIFF::_0_1;
			else if(iUnit==jUnit-1) return &STIFF::_01;
		}
		else if(iLayer == jLayer+1)
		{
			if(iUnit==jUnit) return &STIFF::__10;
			else if(iUnit==jUnit+1) return &STIFF::__1_1;
			else if(iUnit==jUnit-1) return &STIFF::__11;
		}
		else if(iLayer == jLayer-1)
		{
			if(iUnit==jUnit) return &STIFF::_10;
			else if(iUnit==jUnit+1) return &STIFF::_1_1;
			else if(iUnit==jUnit-1) return &STIFF::_11;
		}
		return &STIFF::zero;
	}


	std::vector<MKL_INT>
	range (	MKL_INT width,
		MKL_INT dofPerUnit
	) const noexcept
	{
		auto halfNUnits = nUnits/2;
		width = width<halfNUnits ? width : halfNUnits;
		auto start = (halfNUnits-width) * dofPerUnit;
		auto end = (halfNUnits+width+1) * dofPerUnit;
		std::vector<MKL_INT> ids;
		ids.reserve(end-start);
		for(auto i=start; i!=end; ++i)
			ids.push_back(i);
		return ids;
	}
};


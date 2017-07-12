#ifndef UNITCELLSTIFF_H_INCLUDED
#define UNITCELLSTIFF_H_INCLUDED

#include "array3d.h"
#include <type_traits>
#include <cmath>

template < typename LATTICE,
	typename POTENTIAL >
class UnitCellStiff_2DImp {
public:
	enum {DIM = 2, DOF = LATTICE::N*DIM};


	UnitCellStiff_2DImp (
		const LATTICE& lattice,
		const POTENTIAL& potential,
		double cutoff = 1.01
	) noexcept;


	void
	smallKMatrix (
		int i,
		int j,
		double kMatrix[DOF][DOF]
	) noexcept;


	double
	_11 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m11[i][j];	}


	double
	_10 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m10[i][j];	}


	double
	_1_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m1_1[i][j];	}


	double
	_01 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m01[i][j];	}


	double
	_00 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m00[i][j];	}


	double
	_0_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _01(j,i);	}


	double
	__11 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _1_1(j,i);	}


	double
	__10 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _10(j,i);	}


	double
	__1_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _11(j,i);	}

	//-------------------------------------//

	double
	zero (	unsigned,
		unsigned
	) const noexcept
	{	return 0.;	}

private:
	array3d
	Base (
		int i,
		int j
	) const noexcept;


	array3d
	Diff (
		unsigned i,
		unsigned j,
		const array3d& base
	) const noexcept;

private:
	const LATTICE& lattice;
	const POTENTIAL& potential;
	double cutoff;

	double m11[DOF][DOF];
	double m10[DOF][DOF];
	double m1_1[DOF][DOF];
	double m01[DOF][DOF];
	double m00[DOF][DOF];
};


template < typename LATTICE,
	typename POTENTIAL >
UnitCellStiff_2DImp<LATTICE, POTENTIAL>::UnitCellStiff_2DImp (
	const LATTICE& lattice,
	const POTENTIAL& potential,
	double cutoff
) noexcept :
	lattice( lattice ),
	potential( potential ),
	cutoff( cutoff )
{
	smallKMatrix(1, 1, m11);
	smallKMatrix(1, 0, m10);
	smallKMatrix(1, -1, m1_1);
	smallKMatrix(0, 1, m01);
	smallKMatrix(0, 0, m00);
	constexpr unsigned N = LATTICE::N;
	for(unsigned i=0; i<N; ++i)
	{
		m00[ i*DIM ][ i*DIM ] = 0.;
		m00[ i*DIM ][ i*DIM+1 ] = 0.;
		m00[ i*DIM+1 ][ i*DIM ] = 0.;
		m00[ i*DIM+1 ][ i*DIM+1 ] = 0.;
	}

	for( unsigned i=0; i<DOF; ++i )
	{
		double sum = 0.;
		for( unsigned j=0; j<DOF; ++j )
			sum += _11(i,j) + _10(i,j) + _1_1(i,j) +
				_01(i,j) + _00(i,j) + _0_1(i,j) +
				__11(i,j) + __10(i,j) + __1_1(i,j);
		m00[i][i] = -sum;
	}
}


template < typename LATTICE,
	typename POTENTIAL >
array3d
UnitCellStiff_2DImp<LATTICE, POTENTIAL>::Base (
	int i,
	int j
) const noexcept
{
	return array3d {
		i*lattice.primitive[0][0] + j*lattice.primitive[1][0],
		i*lattice.primitive[0][1] + j*lattice.primitive[1][1],
		0.
	};
}


template < typename LATTICE,
	typename POTENTIAL >
array3d
UnitCellStiff_2DImp<LATTICE, POTENTIAL>::Diff (
	unsigned i, // in cell 0,0
	unsigned j, // in cell base
	const array3d& base
) const noexcept
{
	return array3d {
		-lattice.bases[i][0] + lattice.bases[j][0] + base[0],
		-lattice.bases[i][1] + lattice.bases[j][1] + base[1],
		0.
	};
}


template < typename LATTICE,
	typename POTENTIAL >
void
UnitCellStiff_2DImp<LATTICE, POTENTIAL>::smallKMatrix (
	int iCell,
	int jCell,
	double kMatrix[DOF][DOF]
) noexcept
{
	for(unsigned i=0; i<DOF; ++i)
		for(unsigned j=0; j<DOF; ++j)
			kMatrix[i][j] = 0.;

	constexpr unsigned N = LATTICE::N;
	array3d base = Base( iCell, jCell );

	// pair
	for(unsigned iN=0; iN<N; ++iN)
	for( unsigned jN=0; jN<N; ++jN )
	{
		if(iCell==0 && jCell==0 && iN==jN) continue; // ignore self

		array3d rij = Diff( iN, jN, base );
		if ( Norm(rij) > cutoff ) continue;
		auto k = potential.Hessian(rij);
		kMatrix[ jN*DIM ][ iN*DIM ] += k[0][0];
		kMatrix[ jN*DIM+1 ][ iN*DIM ] += k[0][1];
		kMatrix[ jN*DIM ][ iN*DIM+1 ] += k[1][0];
		kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += k[1][1];
		std::cout<<"k00:"<<k[0][0]<<"\n"<<"k01:"<<k[0][1]<<"\n"<<"k10:"<<k[1][0]<<"\n"<<"k11:"<<k[1][1]<<"\n";
	}

	// angle
	for (unsigned iN=0; iN<N; ++iN)
	for (unsigned jN=0; jN<N; ++jN)
	{
		if(iCell==0 && jCell==0 && iN==jN) continue; // ignore self
		array3d rij = Diff( iN, jN, base );
		auto rij_ = Norm(rij);
		for ( int ii=-1; ii<=1; ++ii )
		for ( int jj=-1; jj<=1; ++jj )
		{
			auto baseK = Base(ii, jj);
			auto baseJK = Base(ii-iCell, jj-jCell);
		for ( unsigned kN=0; kN<N; ++kN )
		{
			array3d rik = Diff( iN, kN, baseK );
			array3d rjk = Diff( jN, kN, baseJK );
			auto rik_ = Norm(rik);
			auto rjk_ = Norm(rjk);
			if ( rik_ < 1.e-5 || rjk_ < 1.e-5 ) continue;

			if ( rij_ < cutoff && rjk_ < cutoff ) // i-j-k
			{
				auto kIK = potential.Hessian(Scale(rij, -1.), rjk);
				kMatrix[ jN*DIM ][ iN*DIM ] += kIK[0][0][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIK[0][0][1];
				kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIK[0][1][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIK[0][1][1];
			}
			else if ( rij_ < cutoff && rik_ < cutoff ) // j-i-k
			{
				auto kIK = potential.Hessian(rij, rik);
				kMatrix[ jN*DIM ][ iN*DIM ] += kIK[0][0][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIK[0][0][1];
				kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIK[0][1][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIK[0][1][1];
			}
			else if ( rik_ < cutoff && rjk_ < cutoff ) // i-k-j
			{
				auto kIK = potential.Hessian(Scale(rik, -1.), Scale(rjk, -1.));
				kMatrix[ jN*DIM ][ iN*DIM ] += kIK[2][0][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIK[2][0][1];
				kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIK[2][1][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIK[2][1][1];
			}
		}
		}
	}
}

//---------------------------------------------------------------------------//

template < typename LATTICE,
	typename POTENTIAL >
class UnitCellStiff_3DImp {
public:
	enum {DIM = 3, DOF=LATTICE::N*DIM};

	UnitCellStiff_3DImp (
		const LATTICE& lattice,
		const POTENTIAL& potential,
		double cutoff = 1.01
	) noexcept;


	void
	smallKMatrix (
		int i,
		int j,
		int k,
		double kMatrix[DOF][DOF]
	) noexcept;


	double
	_111 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m111[i][j];	}


	double
	_110 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m110[i][j];	}


	double
	_11_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m11_1[i][j];	}


	double
	_101 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m101[i][j];	}


	double
	_100 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m100[i][j];	}


	double
	_10_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m10_1[i][j];	}


	double
	_1_11 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m1_11[i][j];	}


	double
	_1_10 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m1_10[i][j];	}


	double
	_1_1_1 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return m1_1_1[i][j];	}


	double
	_011 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m011[i][j];	}


	double
	_010 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m010[i][j];	}


	double
	_01_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m01_1[i][j];	}


	double
	_001 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m001[i][j];	}


	double
	_000 (	unsigned i,
		unsigned j
	) const noexcept
	{	return m000[i][j];	}


	double
	_00_1 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _001(j,i);	}


	double
	_0_11 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _01_1(j,i);	}


	double
	_0_10 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _010(j,i);	}


	double
	_0_1_1 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return _011(j,i);	}


	double
	__111 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _1_1_1(j,i);	}


	double
	__110 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _1_10(j,i);	}


	double
	__11_1 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return _1_11(j,i);	}


	double
	__101 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _10_1(j,i);	}


	double
	__100 (	unsigned i,
		unsigned j
	) const noexcept
	{	return _100(j, i);	}


	double
	__10_1 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return _101(j, i);	}


	double
	__1_11 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return _11_1(j, i);	}


	double
	__1_10 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return _110(j, i);	}


	double
	__1_1_1 (
		unsigned i,
		unsigned j
	) const noexcept
	{	return _111(j, i);	}

	//-------------------------------------//

	double
	zero (	unsigned,
		unsigned
	) const noexcept
	{	return 0.;	}

private:
	array3d
	Base (
		int i,
		int j,
		int k
	) const noexcept;


	array3d
	Diff (
		unsigned i,
		unsigned j,
		const array3d& base
	) const noexcept;


	const LATTICE& lattice;
	const POTENTIAL& potential;
	double cutoff;

	double m111[DOF][DOF];
	double m110[DOF][DOF];
	double m11_1[DOF][DOF];
	double m101[DOF][DOF];
	double m100[DOF][DOF];
	double m10_1[DOF][DOF];
	double m1_11[DOF][DOF];
	double m1_10[DOF][DOF];
	double m1_1_1[DOF][DOF];

	double m011[DOF][DOF];
	double m010[DOF][DOF];
	double m01_1[DOF][DOF];
	double m001[DOF][DOF];
	double m000[DOF][DOF];
};


template < typename LATTICE,
	typename POTENTIAL >
UnitCellStiff_3DImp<LATTICE, POTENTIAL>::UnitCellStiff_3DImp (
	const LATTICE& lattice,
	const POTENTIAL& potential,
	double cutoff
) noexcept :
	lattice( lattice ),
	potential( potential ),
	cutoff( cutoff )
{
	smallKMatrix( 1, 1, 1, m111 );
	smallKMatrix( 1, 1, 0, m110 );
	smallKMatrix( 1, 1, -1, m11_1 );
	smallKMatrix( 1, 0, 1, m101 );
	smallKMatrix( 1, 0, 0, m100 );
	smallKMatrix( 1, 0, -1, m10_1 );
	smallKMatrix( 1, -1, 1, m1_11 );
	smallKMatrix( 1, -1, 0, m1_10 );
	smallKMatrix( 1, -1, -1, m1_1_1 );
	smallKMatrix( 0, 1, 1, m011 );
	smallKMatrix( 0, 1, 0, m010 );
	smallKMatrix( 0, 1, -1, m01_1 );
	smallKMatrix( 0, 0, 1, m001 );
	smallKMatrix( 0, 0, 0, m000 );
	constexpr unsigned N = LATTICE::N;
	for(unsigned i=0; i<N; ++i)
	{
		m000[ i*DIM ][ i*DIM ] = 0.;
		m000[ i*DIM ][ i*DIM+1 ] = 0.;
		m000[ i*DIM ][ i*DIM+2 ] = 0.;
		m000[ i*DIM+1 ][ i*DIM ] = 0.;
		m000[ i*DIM+1 ][ i*DIM+1 ] = 0.;
		m000[ i*DIM+1 ][ i*DIM+2 ] = 0.;
		m000[ i*DIM+2 ][ i*DIM ] = 0.;
		m000[ i*DIM+2 ][ i*DIM+1 ] = 0.;
		m000[ i*DIM+2 ][ i*DIM+2 ] = 0.;
	}
	for( unsigned i=0; i<DOF; ++i )
	{
		double sum = 0.;
		for( unsigned j=0; j<DOF; ++j )
			sum +=	_111(i,j) + _110(i,j) + _11_1(i,j) +
				_101(i,j) + _100(i,j) + _10_1(i,j) +
				_1_11(i,j) + _1_10(i,j) + _1_1_1(i,j) +
				_011(i,j) + _010(i,j) + _01_1(i,j) +
				_001(i,j) + _000(i,j) + _00_1(i,j) +
				_0_11(i,j) + _0_10(i,j) + _0_1_1(i,j) +
				__111(i,j) + __110(i,j) + __11_1(i,j) +
				__101(i,j) + __100(i,j) + __10_1(i,j) +
				__1_11(i,j) + __1_10(i,j) + __1_1_1(i,j);
		m000[i][i] = -sum;
	}
}


template < typename LATTICE,
	typename POTENTIAL >
void
UnitCellStiff_3DImp<LATTICE, POTENTIAL>::smallKMatrix (
	int iCell,
	int jCell,
	int kCell,
	double kMatrix[DOF][DOF]
) noexcept
{
	for(unsigned i=0; i<DOF; ++i)
		for(unsigned j=0; j<DOF; ++j)
			kMatrix[i][j] = 0.;

	constexpr unsigned N = LATTICE::N;
	auto base = Base( iCell, jCell, kCell );
	// pair
	for(unsigned iN=0; iN<N; ++iN)
	for( unsigned jN=0; jN<N; ++jN )
	{
		auto rij = Diff( iN, jN, base );
		if ( Norm(rij) > cutoff ) continue;
		auto kIJ = potential.Hessian( rij );
		kMatrix[ jN*DIM ][ iN*DIM ] += kIJ[0][0];
		kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIJ[0][1];
		kMatrix[ jN*DIM+2 ][ iN*DIM ] += kIJ[0][2];
		kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIJ[1][0];
		kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIJ[1][1];
		kMatrix[ jN*DIM+2 ][ iN*DIM+1 ] += kIJ[1][2];
		kMatrix[ jN*DIM ][ iN*DIM+2 ] += kIJ[2][0];
		kMatrix[ jN*DIM+1 ][ iN*DIM+2 ] += kIJ[2][1];
		kMatrix[ jN*DIM+2 ][ iN*DIM+2 ] += kIJ[2][2];
	}

	// angle
/*	for ( unsigned iN=0; iN<N; ++iN )
	for ( unsigned jN=0; jN<N; ++jN )
	{
		if(iCell==0 && jCell==0 && iN==jN) continue;
		auto rij = Diff( iN, jN, base );
		auto rij_ = Norm(rij);
		for ( int ii=-1; ii<=1; ++ii )
		for ( int jj=-1; jj<=1; ++jj )
		for ( int kk=-1; kk<=1; ++kk )
		{
			auto baseJK = Base( ii-iCell, jj-jCell, kk-kCell );
			auto baseK = Base( ii, jj, kk );
		for ( unsigned kN=0; kN<N; ++kN )
		{
			auto rik = Diff( iN, kN, baseK );
			auto rjk = Diff( jN, kN, baseJK );
			auto rik_ = Norm(rik);
			auto rjk_ = Norm(rjk);
			if ( rik_ < 1.e-5 || rjk_ < 1.e-5 ) continue;

			if ( rij_ < cutoff && rjk_ < cutoff ) // i-j-k
			{
				auto kIK = potential.Hessian(Scale(rij, -1.), rjk);
				kMatrix[ jN*DIM ][ iN*DIM ] += kIK[0][0][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIK[0][0][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM ] += kIK[0][0][2];
				kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIK[0][1][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIK[0][1][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM+1 ] += kIK[0][1][2];
				kMatrix[ jN*DIM ][ iN*DIM+2 ] += kIK[0][2][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+2 ] += kIK[0][2][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM+2 ] += kIK[0][2][2];
			}
			else if ( rij_ < cutoff && rik_ < cutoff ) // j-i-k
			{
				auto kIK = potential.Hessian(rij, rik);
				kMatrix[ jN*DIM ][ iN*DIM ] += kIK[0][0][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIK[0][0][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM ] += kIK[0][0][2];
				kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIK[0][1][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIK[0][1][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM+1 ] += kIK[0][1][2];
				kMatrix[ jN*DIM ][ iN*DIM+2 ] += kIK[0][2][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+2 ] += kIK[0][2][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM+2 ] += kIK[0][2][2];
			}
			else if ( rjk_ < cutoff && rik_ < cutoff ) // i-k-j
			{
				auto kIK = potential.Hessian(Scale(rik, -1.), Scale(rjk, -1.));
				kMatrix[ jN*DIM ][ iN*DIM ] += kIK[2][0][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM ] += kIK[2][0][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM ] += kIK[2][0][2];
				kMatrix[ jN*DIM ][ iN*DIM+1 ] += kIK[2][1][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+1 ] += kIK[2][1][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM+1 ] += kIK[2][1][2];
				kMatrix[ jN*DIM ][ iN*DIM+2 ] += kIK[2][2][0];
				kMatrix[ jN*DIM+1 ][ iN*DIM+2 ] += kIK[2][2][1];
				kMatrix[ jN*DIM+2 ][ iN*DIM+2 ] += kIK[2][2][2];
			}
		}
		}
	}*/
}

template < typename LATTICE,
	typename POTENTIAL >
array3d
UnitCellStiff_3DImp<LATTICE, POTENTIAL>::Base (
	int i,
	int j,
	int k
) const noexcept
{
	return {
		i*lattice.primitive[0][0] + j*lattice.primitive[1][0] + k*lattice.primitive[2][0],
		i*lattice.primitive[0][1] + j*lattice.primitive[1][1] + k*lattice.primitive[2][1],
		i*lattice.primitive[0][2] + j*lattice.primitive[1][2] + k*lattice.primitive[2][2]
	};
}


template < typename LATTICE,
	typename POTENTIAL >
array3d
UnitCellStiff_3DImp<LATTICE, POTENTIAL>::Diff (
	unsigned i, // in cell 0,0,0
	unsigned j, // in cell base
	const array3d& base
) const noexcept
{
	return {
		-lattice.bases[i][0] + lattice.bases[j][0] + base[0],
		-lattice.bases[i][1] + lattice.bases[j][1] + base[1],
		-lattice.bases[i][2] + lattice.bases[j][2] + base[2]
	};
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

template < typename LATTICE,
	typename POTENTIAL >
using UnitCellStiff = std::conditional_t <
	LATTICE::DIM == 2,
	UnitCellStiff_2DImp<LATTICE, POTENTIAL>,
	UnitCellStiff_3DImp<LATTICE, POTENTIAL>
>;

#endif // UNITCELLSTIFF_H_INCLUDED

#include <fstream>
#include <string>
#include <cmath>
#include <memory>
#include "Math.h"
#include <iostream>

struct Potential_EAM {
	double cutoff;
	unsigned nRho;
	double dRho;
	unsigned nR;
	double dR;
	std::unique_ptr<double[]> _rho;
	std::unique_ptr<double[]> _f;
	std::unique_ptr<double[]> _phi;
	double equiliRho;


	Potential_EAM (
		const char* filename = "eam/Au_u3.eam"
	) noexcept
	{
		_rho = nullptr;
		_f = nullptr;
		_phi = nullptr;
		auto lattice = readEAMFile(filename);
		calculateEquiliRho(lattice.second, lattice.first);
	}


	std::pair<std::string, double>
	readEAMFile (
		const char* filename
	) noexcept
	{
		std::ifstream fin(filename);
		std::string sline;
		double dtemp;
		double spacing;
		std::string lattice;
		std::getline(fin, sline);
		fin >> dtemp >> dtemp >> spacing >> lattice;
		std::getline(fin, sline);

		fin >> nRho >> dRho >> nR >> dR >> cutoff;
		allocate();
		for(unsigned i=0; i<nRho; ++i)
			fin >> _f[i];
		for(unsigned i=0; i<nR; ++i)
			fin >> _phi[i];
		for(unsigned i=0; i<nR; ++i)
			fin >> _rho[i];

		effectiveCharge2Phi();
		return std::make_pair(lattice, spacing);
	}


	void
	effectiveCharge2Phi () noexcept
	{
		for( unsigned i=1; i<nR; ++i )
			_phi[i] *= 0.52917721092*27.21138602 * _phi[i] / (dR*i);
	}


	void
	calculateEquiliRho (
		double spacing,
		std::string lattice
	) noexcept
	{
		equiliRho = 0.;
		double cutoffSq = cutoff*cutoff / (spacing*spacing);
		if(lattice.compare("FCC")==0 || lattice.compare("fcc")==0)
		{
			constexpr unsigned nNeighbors[] = {0, 12, 6, 24, 8, 24, 8};
			double rijSq = 0.5;
			unsigned times = 1;
			while(rijSq < cutoffSq)
			{
				double rij = sqrt( rijSq )*spacing;
				equiliRho += rho( rij ) * nNeighbors[times];

				++times;
				rijSq += 0.5;
			}
		}
		else if(lattice.compare("BCC") == 0 || lattice.compare("bcc")==0)
		{
			constexpr unsigned nNeighbors[] = {0, 8, 6, 12, 24, 8};
			constexpr double rijSq[] = {0., 0.75, 1., 2., 2.75, 3.};
			unsigned times = 1;
			while(rijSq[times] < cutoffSq)
			{
				double rij = sqrt( rijSq[times] ) * spacing;
				equiliRho += rho( rij ) * nNeighbors[times];
				++times;
			}
		}
	}


	double
	rho (	double rij ) const noexcept
	{
		return interpolate(rij, _rho.get(), dR, nR);
	}


	void
	allocate () noexcept
	{
		_rho.reset( new double[nR] );
		_f.reset( new double[nRho] );
		_phi.reset( new double[nR] );
	}

	//-------------------------------------------------------------------//

	double
	interpolate (
		double x,
		const double data[],
		double spacing,
		unsigned nMax
	) const noexcept
	{
		Interpolation inter(data, spacing, nMax);
		return inter.interpolate(x);
	}


	double
	diff (	double x,
		const double data[],
		double spacing,
		unsigned nMax
	) const noexcept
	{
		Interpolation inter(data, spacing, nMax);
		return inter.diff(x);
	}


	double
	diffSq (
		double x,
		const double data[],
		double spacing,
		unsigned nMax
	) const noexcept
	{
		Interpolation inter(data, spacing, nMax);
		return inter.diffSq(x);
	}

	//-------------------------------------------------------------------//

	bool
	_2ndDerivative (
		double rij[3],
		double k[3][3]
	) const noexcept
	{
		double rSq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
		double r = sqrt(rSq);
		rij[0] /= r;
		rij[1] /= r;
		rij[2] /= r;
		auto K = _2ndDerivative(r);
		k[0][0] = -K*rij[0]*rij[0];
		k[0][1] = -K*rij[0]*rij[1];
		k[0][2] = -K*rij[0]*rij[2];
		k[1][1] = -K*rij[1]*rij[1];
		k[1][2] = -K*rij[1]*rij[2];
		k[2][2] = -K*rij[2]*rij[2];
		k[1][0] = k[0][1];
		k[2][0] = k[0][2];
		k[2][1] = k[1][2];
		return true;
	}


	double
	_2ndDerivative (
		double rij
	) const noexcept
	{
		if(rij>cutoff) return 0.;
		if(rij<1.e-5) return 0.;
		auto _dSqF = diffSq(equiliRho, _f.get(), dRho, nRho);
		auto _dF = diff(equiliRho, _f.get(), dRho, nRho);
		auto _dSqRho = diffSq(rij, _rho.get(), dR, nR);
		auto _dRho = diff(rij, _rho.get(), dR, nR);
		auto _dSqPhi = diffSq(rij, _phi.get(), dR, nR);

		return _dSqF*_dRho*_dRho + _dF*_dSqRho + _dSqPhi;
	}


	bool
	_2ndDerivative (
		double [3],
		double [3],
		double [3][3]
	) const noexcept
	{
		return false;
	}
};


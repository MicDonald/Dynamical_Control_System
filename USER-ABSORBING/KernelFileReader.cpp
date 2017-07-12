#include "KernelFileReader.h"
#include <cmath>
#include <iostream>

using namespace std;


KernelFileReaderImp::KernelFileReaderImp (
	const char* filename
) noexcept
{
	fin.open( filename, ios::in | ios::binary );
}


void
KernelFileReaderImp::Read (
	double tCut,
	double timeStep,
	KernelFunctionsMap& kernelFunctionsMap
) noexcept
{
	double dt, tAll;
	fin.read( (char*)&dt, sizeof(double) );
	fin.read( (char*)&tAll, sizeof(double) );
	auto nData = unsigned(tAll/dt + 0.01);
	auto nDataSize = nData * sizeof(double);
	unique_ptr<double[]> kernelFunction( new double[nData] );
	unsigned nUsed = round( (tCut<tAll ?  tCut : tAll)/timeStep );
//cout<<"Kernel File: \n";
	while ( !fin.eof() )
	{
		InterfacePairConfiguration config;
		char orientation;
		fin.read( &orientation, sizeof(char) );
		if ( fin.eof() ) break;
		config.orientation = orientation == '0' ? REGULAR : SHIFT;
		fin.read( (char*)config.outerCellID.data(), sizeof(int)*3 );
		fin.read( (char*)&config.outerLocal, sizeof(unsigned) );
		fin.read( (char*)config.innerCellID.data(), sizeof(int)*3 );
		fin.read( (char*)&config.innerLocal, sizeof(unsigned) );
//cout<<config.orientation<<": "
//<<config.outerCellID[0]<<"_"<<config.outerCellID[1]<<"_"<<config.outerCellID[2]<<"_"<<config.outerLocal<<", "
//<<config.innerCellID[0]<<"_"<<config.innerCellID[1]<<"_"<<config.innerCellID[2]<<"_"<<config.innerLocal<<", : (";

		shared_ptr<KernelFunctions> functionsPtr( new KernelFunctions );
		unsigned nCurves;
		fin.read( (char*)&nCurves, sizeof(unsigned) );
		for ( unsigned i=0; i<nCurves; ++i )
		{
			char type;// xx xy xz yx yy yz zx zy zz
			fin.read( &type, sizeof(char) );
			fin.read( (char*)kernelFunction.get(), nDataSize );
			switch( type )
			{
			case '0': // xx
//cout<<"xx, ";
				functionsPtr->xx.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->xx.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->xx.get(), dt, nUsed );
				break;
			case '1': // xy
//cout<<"xy, ";
				functionsPtr->xy.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->xy.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->xy.get(), dt, nUsed );
				break;
			case '2': // xz
//cout<<"xz, ";
				functionsPtr->xz.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->xz.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->xz.get(), dt, nUsed );
				break;
			case '3': // yx
//cout<<"yx, ";
				functionsPtr->yx.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->yx.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->yx.get(), dt, nUsed );
				break;
			case '4': // yy
//cout<<"yy, ";
				functionsPtr->yy.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->yy.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->yy.get(), dt, nUsed );
				break;
			case '5': // yz
//cout<<"yz, ";
				functionsPtr->yz.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->yz.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->yz.get(), dt, nUsed );
				break;
			case '6': // zx
//cout<<"zx, ";
				functionsPtr->zx.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->zx.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->zx.get(), dt, nUsed );
				break;
			case '7': // zy
//cout<<"zy, ";
				functionsPtr->zy.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->zy.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->zy.get(), dt, nUsed );
				break;
			case '8': // zz
//cout<<"zz, ";
				functionsPtr->zz.reset( new double[nUsed] );
				Interpolate( kernelFunction.get(), nData, dt, functionsPtr->zz.get(), nUsed, timeStep );
				NumericalIntegralScaling( functionsPtr->zz.get(), dt, nUsed );
				break;
			}
		}
//cout<<")"<<endl;

		kernelFunctionsMap.insert( make_pair (
			move(config), functionsPtr
		) );
	}
}


void
KernelFileReaderImp::NumericalIntegralScaling (
	double kernelFunction[],
	double scalingFactor,
	unsigned nData
) const noexcept
{
	for ( unsigned i=0; i<nData; ++i )
		kernelFunction[i] *= scalingFactor;
}


void
KernelFileReaderImp::Interpolate (
	const double function[],
	unsigned nData,
	double dt,
	double interpolateFunction[],
	unsigned nUsed,
	double timeStep
) const noexcept
{
	for ( unsigned i=1; i<=nUsed; ++i )
		interpolateFunction[i-1] = Interpolate (
			function, nData, dt, i*timeStep
		);
}


double
KernelFileReaderImp::Interpolate (
	const double function[],
	unsigned nData,
	double dx,
	double x
) const noexcept
{
	if ( x > dx*nData ) return 0.;

	double inverseDx = 1./dx;
	unsigned leftPoint = x * inverseDx;
	double extDx = x - leftPoint*dx;
	double leftFx = leftPoint==0 ? 0. : function[leftPoint-1];
	double rightFx = leftPoint>=nData ? 0. : function[leftPoint];
	double dF = rightFx - leftFx;
	return leftFx + dF * extDx * inverseDx;
}


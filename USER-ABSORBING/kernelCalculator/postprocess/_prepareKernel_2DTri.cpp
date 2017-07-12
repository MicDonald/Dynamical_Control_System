#include <iostream>
#include <fstream>
#include <string>
#include <ios>
#include <iomanip>
#include <vector>

using namespace std;

enum TYPE {XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };

void ReadAndWrite (
	ofstream& fout,
	char inFilename[],
	char type,
	int nData,
	bool negative
)
{
	ifstream fin( inFilename );
	if(!fin.is_open()) {
		cout<<"file "<<inFilename<<" not existed."<<endl;
		return;
	}
	fout.write(&type, sizeof(char));
	double kernel;
	for(int k=0; k<nData; ++k) {
		fin >> kernel;
		if( negative ) kernel = -kernel;
		fout.write((char*)&kernel, sizeof(double));
	}
}


int main()
{
	double dt=0.01, tall=50.;
	int nsides = 40;
	int nsidesV = 25;
	int ncorners_ = -6;
	int center = 35;
	int ncorners = 10;
	string outFilename = "_2D_Tri";

	ofstream fout(outFilename.c_str(), ios::out | ios::binary);
	fout.write((char*)&dt, sizeof(double));
	fout.write((char*)&tall, sizeof(double));

	char inFilename[100];
	int nData = int(tall/dt+0.01);
	char regular = '0', shift = '1';
	int one = 1;
	int izero = 0;
	int pbc = -1;
	bool positive = false, negative = true;
	//side -
//	for(int i=-nsides; i<=nsides; ++i){
	for(int i=0; i<nsides; ++i){
		unsigned bcLocal = 1, innerLocal = 0;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&pbc, sizeof(int));///	-1_0_-1_1
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&i, sizeof(int));///		i_1_0_0
		fout.write((char*)&one, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		int nCurves = 4;
		fout.write((char*)&nCurves, sizeof(int));

		sprintf(inFilename, "tri_side/50_1y-%d_0y", 49+i);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_side/50_1y-%d_0x", 49+i);//xy
		ReadAndWrite(fout, inFilename, '1', nData, true);
		sprintf(inFilename, "tri_side/50_1x-%d_0y", 49+i);//yx
		ReadAndWrite(fout, inFilename, '3', nData, true);
		sprintf(inFilename, "tri_side/50_1x-%d_0x", 49+i);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);
	}
	//side |
//	for(int i=-nsidesV; i<=nsidesV; ++i){
	for(int i=0; i<nsidesV; ++i){
		unsigned bcLocal = 0, innerLocal = 0;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&one, sizeof(int));///
		fout.write((char*)&i, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		int nCurves = 4;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_sideV/50_0x-%d_0x", 50-i);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_sideV/50_0x-%d_0y", 50-i);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_sideV/50_0y-%d_0x", 50-i);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "tri_sideV/50_0y-%d_0y", 50-i);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);

		//===============================================

		bcLocal = 0;
		innerLocal = 1;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&one, sizeof(int));///
		fout.write((char*)&i, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		nCurves = 2;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_sideV/50_0x-%d_1x", 50-i);
		ReadAndWrite(fout, inFilename, '0', nData, false);
//		sprintf(inFilename, "tri_sideV/50_0x-%d_1y", 50-i);//no
//		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_sideV/50_0y-%d_1x", 50-i);
		ReadAndWrite(fout, inFilename, '3', nData, false);
//		sprintf(inFilename, "tri_sideV/50_0y-%d_1y", 50-i);//no
//		ReadAndWrite(fout, inFilename, '4', nData, false);

		//===============================================

		bcLocal = 1;
		innerLocal = 0;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&one, sizeof(int));///
		fout.write((char*)&i, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		nCurves = 4;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_sideV/50_1x-%d_0x", 50-i);
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_sideV/50_1x-%d_0y", 50-i);
		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_sideV/50_1y-%d_0x", 50-i);
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "tri_sideV/50_1y-%d_0y", 50-i);
		ReadAndWrite(fout, inFilename, '4', nData, false);

		//===============================================

		bcLocal = 1;
		innerLocal = 1;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&one, sizeof(int));///
		fout.write((char*)&i, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		nCurves = 2;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_sideV/50_1x-%d_1x", 50-i);
		ReadAndWrite(fout, inFilename, '0', nData, false);
//		sprintf(inFilename, "tri_sideV/50_1x-%d_1y", 50-i);//no
//		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_sideV/50_1y-%d_1x", 50-i);
		ReadAndWrite(fout, inFilename, '3', nData, false);
//		sprintf(inFilename, "tri_sideV/50_1y-%d_1y", 50-i);//no
//		ReadAndWrite(fout, inFilename, '4', nData, false);
	}

	// corner L
	for(int i=ncorners_; i<=ncorners; ++i)
	for(int j=ncorners_+1; j<=ncorners-1; ++j) {
		int outerX = i < 0 ? 0 : i;
		int outerY = i < 0 ? -i : 0;
		int innerX = j < 0 ? 1 : j+1;
		int innerY = j < 0 ? -j+1 : 1;

		unsigned bcLocal = 0, innerLocal = 0;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_0
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_0
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		int nCurves = 4;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerL/%d_0x-%d_0x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_0x-%d_0y", center+i, center+j);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_0y-%d_0x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_0y-%d_0y", center+i, center+j);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);

		//===============================================

		bcLocal = 0;
		innerLocal = 1;
		nCurves = j<=0 ? 2 : 0;
if ( nCurves > 0 )
{
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_0
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_1
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerL/%d_0x-%d_1x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
//		sprintf(inFilename, "tri_cornerL/%d_0x-%d_1y", center+i, center+j);//xy
//		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_0y-%d_1x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
//		sprintf(inFilename, "tri_cornerL/%d_0y-%d_1y", center+i, center+j);//yy
//		ReadAndWrite(fout, inFilename, '4', nData, false);
}
		//===============================================

		bcLocal = 1;
		innerLocal = 0;
		nCurves = 4;
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_1
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_0
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerL/%d_1x-%d_0x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_1x-%d_0y", center+i, center+j);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_1y-%d_0x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_1y-%d_0y", center+i, center+j);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);

		//===============================================

		bcLocal = 1;
		innerLocal = 1;
		nCurves = j<=0 ? 2 : 0;
if ( nCurves > 0 )
{
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_1
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_1
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerL/%d_1x-%d_1x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
//		sprintf(inFilename, "tri_cornerL/%d_1x-%d_1y", center+i, center+j);//xy
//		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_cornerL/%d_1y-%d_1x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
//		sprintf(inFilename, "tri_cornerL/%d_1y-%d_1y", center+i, center+j);//yy
//		ReadAndWrite(fout, inFilename, '4', nData, false);
}
	}

	// corner O
	for(int i=ncorners_; i<=ncorners; ++i)
	for(int j=ncorners_+1; j<=ncorners-1; ++j) {
		int outerX = i < 0 ? 0 : i;
		int outerY = i < 0 ? -i : 0;
		int innerX = j < 0 ? 1 : j+1;
		int innerY = j < 0 ? -j+1 : 1;

		unsigned bcLocal = 0, innerLocal = 0;
		fout.write((char*)&shift, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_0
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_0
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		int nCurves = i<0 ? 2 : 4;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerO/%d_0x-%d_0x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
if ( i>=0 ) {
		sprintf(inFilename, "tri_cornerO/%d_0x-%d_0y", center+i, center+j);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
}
		sprintf(inFilename, "tri_cornerO/%d_0y-%d_0x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
if ( i>=0 ) {
		sprintf(inFilename, "tri_cornerO/%d_0y-%d_0y", center+i, center+j);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);
}

		//===============================================

		bcLocal = 0; innerLocal = 1;
		nCurves = i<=0 ? 4 : 0;
if ( nCurves > 0 )
{
		fout.write((char*)&shift, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_0
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_1
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerO/%d_0x-%d_1x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_cornerO/%d_0x-%d_1y", center+i, center+j);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_cornerO/%d_0y-%d_1x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "tri_cornerO/%d_0y-%d_1y", center+i, center+j);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);
}

		//===============================================

		bcLocal = 1; innerLocal = 0;
		fout.write((char*)&shift, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_1
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_0
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		nCurves = i<0 ? 2 : 4;
		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerO/%d_1x-%d_0x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
if ( i>=0 ) {
		sprintf(inFilename, "tri_cornerO/%d_1x-%d_0y", center+i, center+j);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
}
		sprintf(inFilename, "tri_cornerO/%d_1y-%d_0x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
if ( i>=0 ) {
		sprintf(inFilename, "tri_cornerO/%d_1y-%d_0y", center+i, center+j);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);
}

		//===============================================

		bcLocal = 1; innerLocal = 1;
		nCurves = i<=0 ? 4 : 0;
if ( nCurves > 0 )
{
		fout.write((char*)&shift, sizeof(char));
		fout.write((char*)&outerX, sizeof(int));///	x_y_-1_1
		fout.write((char*)&outerY, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&innerX, sizeof(int));///	x_y_0_1
		fout.write((char*)&innerY, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		fout.write((char*)&nCurves, sizeof(int));
		sprintf(inFilename, "tri_cornerO/%d_1x-%d_1x", center+i, center+j);//xx
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "tri_cornerO/%d_1x-%d_1y", center+i, center+j);//xy
		ReadAndWrite(fout, inFilename, '1', nData, false);
		sprintf(inFilename, "tri_cornerO/%d_1y-%d_1x", center+i, center+j);//yx
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "tri_cornerO/%d_1y-%d_1y", center+i, center+j);//yy
		ReadAndWrite(fout, inFilename, '4', nData, false);
}
	}

	fout.close();		
}


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
	double dt=0.001, tall=10.;
	int nsides = 5;
//	int ncorners_ = -6;
//	int ncorners = 10;
	string outFilename = "_3D_Ar";

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
	for(int i=0; i<=nsides; ++i)
	for(int j=0; j<=nsides; ++j)
	for(unsigned bcLocal=0; bcLocal<=3; ++bcLocal)
	for(unsigned innerLocal=0; innerLocal<=3; ++innerLocal) {
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&izero, sizeof(int));///	0_-1_-1_0
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&one, sizeof(int));///	1_j_i_0
		fout.write((char*)&j, sizeof(int));///
		fout.write((char*)&i, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		int nCurves = 9;
		fout.write((char*)&nCurves, sizeof(int));

		sprintf(inFilename, "Ar_face/144_%dx-%d_%dx", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "Ar_face/144_%dx-%d_%dy", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '1', nData, true);
		sprintf(inFilename, "Ar_face/144_%dx-%d_%dz", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '2', nData, true);
		sprintf(inFilename, "Ar_face/144_%dy-%d_%dx", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "Ar_face/144_%dy-%d_%dy", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '4', nData, false);
		sprintf(inFilename, "Ar_face/144_%dy-%d_%dz", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '5', nData, false);
		sprintf(inFilename, "Ar_face/144_%dz-%d_%dx", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '6', nData, false);
		sprintf(inFilename, "Ar_face/144_%dz-%d_%dy", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '7', nData, false);
		sprintf(inFilename, "Ar_face/144_%dz-%d_%dz", bcLocal, 144+i*17+j, innerLocal);
		ReadAndWrite(fout, inFilename, '8', nData, false);
	}
}


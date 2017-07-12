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
	double dt=0.001, tall=200.;
	int nsides = 1,i=0,j=0;
//	int ncorners_ = -6;
//	int ncorners = 10;
	string outFilename = "_2D_Squ";

	ofstream fout(outFilename.c_str(), ios::out | ios::binary);
	fout.write((char*)&dt, sizeof(double));
	fout.write((char*)&tall, sizeof(double));

	char inFilename[100];
	int nData = int(tall/dt+0.001);
	char regular = '0', shift = '1';
	int one = 1;
	int izero = 0;
	int pbc = -1;
	bool positive = false, negative = true;
	
	for(unsigned bcLocal=0; bcLocal<=0; ++bcLocal)
	for(unsigned innerLocal=0; innerLocal<=0; ++innerLocal) {
		fout.write((char*)&regular, sizeof(char));
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&izero, sizeof(int));///	0_-1_0_-1
		fout.write((char*)&pbc, sizeof(int));///
		fout.write((char*)&bcLocal, sizeof(unsigned));
		fout.write((char*)&izero, sizeof(int));///	0_1_i_0
		fout.write((char*)&one, sizeof(int));///
		fout.write((char*)&i, sizeof(int));///
		fout.write((char*)&innerLocal, sizeof(unsigned));

		int nCurves = 9;
		fout.write((char*)&nCurves, sizeof(int));

		sprintf(inFilename, "k1/1_0z-1_0z", bcLocal, 1, innerLocal);
		ReadAndWrite(fout, inFilename, '0', nData, false);
		sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '1', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '2', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '3', nData, false);
		sprintf(inFilename, "k1/1_0z-1_0z", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '4', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '5', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '6', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '7', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '8', nData, false);
	}
	one = 12;
        for(unsigned bcLocal=0; bcLocal<=0; ++bcLocal)
        for(unsigned innerLocal=0; innerLocal<=0; ++innerLocal) {
                fout.write((char*)&regular, sizeof(char));
                fout.write((char*)&pbc, sizeof(int));///
                fout.write((char*)&izero, sizeof(int));///      0_-1_0_-1
                fout.write((char*)&pbc, sizeof(int));///
                fout.write((char*)&bcLocal, sizeof(unsigned));
                fout.write((char*)&izero, sizeof(int));///      0_1_i_0
                fout.write((char*)&one, sizeof(int));///
                fout.write((char*)&i, sizeof(int));///
                fout.write((char*)&innerLocal, sizeof(unsigned));

                int nCurves = 9;
                fout.write((char*)&nCurves, sizeof(int));

                sprintf(inFilename, "kN/1_0z-1_0z", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '0', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '1', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '2', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '3', nData, false);
                sprintf(inFilename, "kN/1_0z-1_0z", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '4', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '5', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '6', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '7', nData, false);
                sprintf(inFilename, "noKernel", bcLocal, 1, innerLocal);
                ReadAndWrite(fout, inFilename, '8', nData, false);
        }        

}


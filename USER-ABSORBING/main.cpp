#include "InterfaceCellIdentifier.h"
#include "Lattice.h"
#include "main.h"
#include "data3d.h"
#include "KernelRelation.h"

#include "KernelFileReader.h"
#include <vector>
#include <iostream>
#include <exception>
using namespace std;

PBC pbc;
PBCDifference nbc;

void CellizeTest_Tri ()
{
	width = 10;
	heigh = 10;
	box[0] = width*lattice;
	box[1] = heigh*sqrt(3.)*lattice;
	vector<data3d> outerPos;
	vector<data3d> innerPos;

//	Create_Side(width, heigh, outerPos, innerPos);
//	Create_SideV(width, heigh, outerPos, innerPos);
//	InterfaceCellIdentifier<Lattice_Triangular> identifier( &pbc );

	Create_Corner1(width, heigh, outerPos, innerPos);
//	Create_Corner1_(width, heigh, outerPos, innerPos);
//	Create_Corner2(width, heigh, outerPos, innerPos);
//	Create_Corner2_(width, heigh, outerPos, innerPos);
	InterfaceCellIdentifier<Lattice_Triangular> identifier( &nbc );

cout<<"outer:\n";
for(const auto& pos : outerPos)
cout<<pos[0]<<" "<<pos[1]<<"\n";
cout<<"\n";
cout<<"inner:\n";
for(const auto& pos : innerPos)
cout<<pos[0]<<" "<<pos[1]<<"\n";
cout<<endl;
	identifier.Identify( innerPos, outerPos );
}

void CellizeTest_Si ()
{
	width=5;
	heigh=5;
	deep=5;
	box[0] = width*lattice;
	box[1] = heigh*lattice;
	box[2] = width*lattice;
	vector<data3d> outerPos;
	vector<data3d> innerPos;

//	Create_SideX_Si( width, heigh, deep, outerPos, innerPos );
//	Create_SideY_Si( width, heigh, deep, outerPos, innerPos );
//	Create_SideZ_Si( width, heigh, deep, outerPos, innerPos );

//	box[0] *= 10.; box[1] *= 10.;
//	Create_EdgeXY_Si( width, heigh, deep, outerPos, innerPos );

//	box[0] *= 10.; box[2] *= 10.;
//	Create_EdgeXZ_Si( width, heigh, deep, outerPos, innerPos );

//	box[1] *= 10.; box[2] *= 10.;
//	Create_EdgeYZ_Si( width, heigh, deep, outerPos, innerPos );

	box[0] *= 10.; box[1] *= 10.; box[2] *= 10.;
	Create_Corner_Si( width, heigh, deep, outerPos, innerPos );

	InterfaceCellIdentifier<Lattice_Si> identifier( &pbc );

/*cout<<"outer:\n";
for(unsigned i=0; i<outerPos.size(); ++i)
cout<<i<<" "<<outerPos[i][0]<<" "<<outerPos[i][1]<<" "<<outerPos[i][2]<<"\n";
cout<<"\n";*/
/*cout<<"inner:\n";
for(unsigned i=0; i<innerPos.size(); ++i)
cout<<i<<" "<<innerPos[i][0]<<" "<<innerPos[i][1]<<" "<<innerPos[i][2]<<"\n";
cout<<endl;*/
	identifier.Identify( innerPos, outerPos );
}

void CellizeTest_FCC ()
{
	width = 5;
	heigh = 5;
	deep = 5;
	box[0] = width*lattice;
	box[1] = heigh*lattice;
	box[2] = deep*lattice;
	vector<data3d> outerPos;
	vector<data3d> innerPos;

	Create_SideX_FCC(width, heigh, deep, outerPos, innerPos);

	InterfaceCellIdentifier<Lattice_FCC> identifier( &pbc );

cout<<"inner:\n";
for(unsigned i=0; i<innerPos.size(); ++i)
cout<<i<<" "<<innerPos[i][0]<<" "<<innerPos[i][1]<<" "<<innerPos[i][2]<<"\n";
cout<<endl;
/*cout<<"outer:\n";
for(unsigned i=0; i<outerPos.size(); ++i)
cout<<i<<" "<<outerPos[i][0]<<" "<<outerPos[i][1]<<" "<<outerPos[i][2]<<"\n";
cout<<"\n";*/
	identifier.Identify( innerPos, outerPos );
}

void PrintConfig ( const InterfacePairConfiguration& c )
{
	cout << "["<<c.orientation<<"] "
	<<c.outerCellID[0]<<"_"<<c.outerCellID[1]<<"_"<<c.outerCellID[2]<<"_"<<c.outerLocal<<" "
	<<c.innerCellID[0]<<"_"<<c.innerCellID[1]<<"_"<<c.innerCellID[2]<<"_"<<c.innerLocal<<" "
	<<"["<<c.sign[0]<<","<<c.sign[1]<<","<<c.sign[2]<<"]\t"
	<<"{"<<c.rotateMatrix[0][0]<<" "<<c.rotateMatrix[1][1]<<" "<<c.rotateMatrix[2][2]<<"}"<<endl;
}

void Outer2InnerPairTest()
{
	vector<data3d> outerPos;
	vector<data3d> innerPos;
	KernelRelation<Lattice_Triangular> relation(10., 10., 10., 0.01);
	NBC npbc;

//	Create_Side(width, heigh, outerPos, innerPos);
//	Create_SideV(width, heigh, outerPos, innerPos);
//	relation.SetInitialState( innerPos, outerPos, &pbc );

	Create_Corner1(width, heigh, outerPos, innerPos);
//	Create_Corner1_(width, heigh, outerPos, innerPos);
//	Create_Corner2(width, heigh, outerPos, innerPos);
//	Create_Corner2_(width, heigh, outerPos, innerPos);
	relation.SetInitialState( innerPos, outerPos, &npbc );
for ( unsigned i=0; i<outerPos.size(); ++i)
{
	const auto& kernelPair = relation.KernelFunctionsRelationOf(i);
	for ( const auto& e : kernelPair )
	{
		const auto& c = e.second.config;
		cout << i <<" " << e.first<<": ";
		PrintConfig( c );
	}
}
}


void ReadKernelFunctionTest ()
{
	KernelFileReader<Lattice_Triangular> reader;
	KernelFunctionsMap kernelFunction;
	reader.ReadKernelFile ( 30, 0.01, kernelFunction );
/*	InterfacePairConfiguration config;
	config.orientation = REGULAR;
	config.outerCellID = {-1, 0, -1};
	config.outerLocal = 1;
	config.innerCellID = {-1, 1, 0};
	config.innerLocal = 0;
	auto kernelPtr = kernelFunction.at(config);
	const auto* function = kernelPtr->yy.get();
	for ( int i=0; i<3000; ++i )
		cout << (i+1)*0.01 <<"\t"<< function[i] << endl;*/
}

template < typename LATTICE >
void ConfigAliasTest2D ()
{
	InterfacePairConfiguration config11a;
	config11a.orientation = REGULAR;
	config11a.outerCellID = {0, -1, -1};
	config11a.outerLocal = 0;
	config11a.innerCellID = {-1, 0, 0};
	config11a.innerLocal = 0;
	PrintConfig( config11a );
	InterfacePairConfigurationTransformer()( config11a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config11a );
	cout<<endl;
	InterfacePairConfiguration config11b;
	config11b.orientation = REGULAR;
	config11b.outerCellID = {0, -1, -1};
	config11b.outerLocal = 0;
	config11b.innerCellID = {-1, 1, 0};
	config11b.innerLocal = 1;
	PrintConfig( config11b );
	InterfacePairConfigurationTransformer()( config11b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config11b );
	cout<<endl;
	InterfacePairConfiguration config11c;
	config11c.orientation = REGULAR;
	config11c.outerCellID = {0, -1, -1};
	config11c.outerLocal = 1;
	config11c.innerCellID = {-1, 0, 0};
	config11c.innerLocal = 0;
	PrintConfig( config11c );
	InterfacePairConfigurationTransformer()( config11c,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config11c );
	cout<<endl;
	InterfacePairConfiguration config11d;
	config11d.orientation = REGULAR;
	config11d.outerCellID = {0, -1, -1};
	config11d.outerLocal = 1;
	config11d.innerCellID = {-1, 2, 0};
	config11d.innerLocal = 1;
	PrintConfig( config11d );
	InterfacePairConfigurationTransformer()( config11d,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config11d );
	cout<<endl;

	InterfacePairConfiguration config12a;
	config12a.orientation = REGULAR;
	config12a.outerCellID = {0, -1, -1};
	config12a.outerLocal = 0;
	config12a.innerCellID = {1, -1, 0};
	config12a.innerLocal = 1;
	PrintConfig( config12a );
	InterfacePairConfigurationTransformer()( config12a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config12a );
	cout << endl;
	InterfacePairConfiguration config12b;
	config12b.orientation = REGULAR;
	config12b.outerCellID = {0, -1, -1};
	config12b.outerLocal = 1;
	config12b.innerCellID = {1, -2, 0};
	config12b.innerLocal = 1;
	PrintConfig( config12b );
	InterfacePairConfigurationTransformer()( config12b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config12b );
	cout << endl;

	InterfacePairConfiguration config13a;
	config13a.orientation = REGULAR;
	config13a.outerCellID = {0, -1, -1};
	config13a.outerLocal = 1;
	config13a.innerCellID = {-1, -1, 0};
	config13a.innerLocal = 0;
	PrintConfig( config13a );
	InterfacePairConfigurationTransformer()( config13a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config13a );
	cout<<endl;
	InterfacePairConfiguration config13b;
	config13b.orientation = REGULAR;
	config13b.outerCellID = {0, -1, -1};
	config13b.outerLocal = 0;
	config13b.innerCellID = {-1, -2, 0};
	config13b.innerLocal = 0;
	PrintConfig( config13b );
	InterfacePairConfigurationTransformer()( config13b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config13b );
	cout<<endl;

	InterfacePairConfiguration config21a;
	config21a.orientation = REGULAR;
	config21a.outerCellID = {-1, 0, -1};
	config21a.outerLocal = 0;
	config21a.innerCellID = {-1, 1, 0};
	config21a.innerLocal = 1;
	PrintConfig( config21a );
	InterfacePairConfigurationTransformer()( config21a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config21a );
	cout << endl;
	InterfacePairConfiguration config21b;
	config21b.orientation = REGULAR;
	config21b.outerCellID = {-1, 0, -1};
	config21b.outerLocal = 1;
	config21b.innerCellID = {-2, 1, 0};
	config21b.innerLocal = 0;
	PrintConfig( config21b );
	InterfacePairConfigurationTransformer()( config21b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config21b );
	cout << endl;

	InterfacePairConfiguration config22a;
	config22a.orientation = REGULAR;
	config22a.outerCellID = {-1, 0, -1};
	config22a.outerLocal = 0;
	config22a.innerCellID = {0, -1, 0};
	config22a.innerLocal = 1;
	PrintConfig( config22a );
	InterfacePairConfigurationTransformer()( config22a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config22a );
	cout<<endl;
	InterfacePairConfiguration config22b;
	config22b.orientation = REGULAR;
	config22b.outerCellID = {-1, 0, -1};
	config22b.outerLocal = 1;
	config22b.innerCellID = {0, -1, 0};
	config22b.innerLocal = 0;
	PrintConfig( config22b );
	InterfacePairConfigurationTransformer()( config22b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config22b );
	cout<<endl;
	InterfacePairConfiguration config22c;
	config22c.orientation = REGULAR;
	config22c.outerCellID = {-1, 0, -1};
	config22c.outerLocal = 1;
	config22c.innerCellID = {1, -1, 0};
	config22c.innerLocal = 1;
	PrintConfig( config22c );
	InterfacePairConfigurationTransformer()( config22c,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config22c );
	cout<<endl;
	InterfacePairConfiguration config22d;
	config22d.orientation = REGULAR;
	config22d.outerCellID = {-1, 0, -1};
	config22d.outerLocal = 0;
	config22d.innerCellID = {2, -1, 0};
	config22d.innerLocal = 1;
	PrintConfig( config22d );
	InterfacePairConfigurationTransformer()( config22d,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config22d );
	cout<<endl;

	InterfacePairConfiguration config23a;
	config23a.orientation = REGULAR;
	config23a.outerCellID = {-1, 0, -1};
	config23a.outerLocal = 0;
	config23a.innerCellID = {-1, -1, 0};
	config23a.innerLocal = 0;
	PrintConfig( config23a );
	InterfacePairConfigurationTransformer()( config23a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config23a );
	cout << endl;
	InterfacePairConfiguration config23b;
	config23b.orientation = REGULAR;
	config23b.outerCellID = {-1, 0, -1};
	config23b.outerLocal = 0;
	config23b.innerCellID = {-2, -1, 0};
	config23b.innerLocal = 0;
	PrintConfig( config23b );
	InterfacePairConfigurationTransformer()( config23b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config23b );
	cout << endl;

	InterfacePairConfiguration config31a;
	config31a.orientation = REGULAR;
	config31a.outerCellID = {1, 0, -1};
	config31a.outerLocal = 0;
	config31a.innerCellID = {0, 2, 0};
	config31a.innerLocal = 1;
	config31a.sign = {1, 0, 0};
	PrintConfig( config31a );
	InterfacePairConfigurationTransformer()( config31a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config31a );
	cout<<endl;
	InterfacePairConfiguration config31b;
	config31b.orientation = REGULAR;
	config31b.outerCellID = {0, 0, -1};
	config31b.outerLocal = 1;
	config31b.innerCellID = {-2, 1, 0};
	config31b.innerLocal = 1;
	config31b.sign = {0, 0, 0};
	PrintConfig( config31b );
	InterfacePairConfigurationTransformer()( config31b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config31b );
	cout<<endl;
	InterfacePairConfiguration config31c;
	config31c.orientation = SHIFT;
	config31c.outerCellID = {2, 0, -1};
	config31c.outerLocal = 1;
	config31c.innerCellID = {-1, 2, 0};
	config31c.innerLocal = 0;
	config31c.sign = {1, 0, 0};
	PrintConfig( config31c );
	InterfacePairConfigurationTransformer()( config31c,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config31c );
	cout<<endl;
	InterfacePairConfiguration config31d;
	config31d.orientation = SHIFT;
	config31d.outerCellID = {0, 2, -1};
	config31d.outerLocal = 1;
	config31d.innerCellID = {-1, 0, 0};
	config31d.innerLocal = 0;
	config31d.sign = {0, 0, 0};
	PrintConfig( config31d );
	InterfacePairConfigurationTransformer()( config31d,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config31d );
	cout<<endl;

	InterfacePairConfiguration config32a;
	config32a.orientation = SHIFT;
	config32a.outerCellID = {0, 2, -1};
	config32a.outerLocal = 1;
	config32a.innerCellID = {1, -1, 0};
	config32a.innerLocal = 0;
	config32a.sign = {0, 1, 0};
	PrintConfig( config32a );
	InterfacePairConfigurationTransformer()( config32a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config32a );
	cout<<endl;
	InterfacePairConfiguration config32b;
	config32b.orientation = REGULAR;
	config32b.outerCellID = {0, 0, -1};
	config32b.outerLocal = 1;
	config32b.innerCellID = {1, -1, 0};
	config32b.innerLocal = 0;
	config32b.sign = {0, 0, 0};
	PrintConfig( config32b );
	InterfacePairConfigurationTransformer()( config32b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config32b );
	cout<<endl;
	InterfacePairConfiguration config32c;
	config32c.orientation = SHIFT;
	config32c.outerCellID = {2, 0, -1};
	config32c.outerLocal = 1;
	config32c.innerCellID = {-1, -2, 0};
	config32c.innerLocal = 0;
	config32c.sign = {0, 0, 0};
	PrintConfig( config32c );
	InterfacePairConfigurationTransformer()( config32c,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config32c );
	cout<<endl;
	InterfacePairConfiguration config32d;
	config32d.orientation = REGULAR;
	config32d.outerCellID = {0, 1, -1};
	config32d.outerLocal = 1;
	config32d.innerCellID = {2, 0, 0};
	config32d.innerLocal = 0;
	config32d.sign = {0, 1, 0};
	PrintConfig( config32d );
	InterfacePairConfigurationTransformer()( config32d,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config32d );
	cout<<endl;


	InterfacePairConfiguration config33a;
	config33a.orientation = SHIFT;
	config33a.outerCellID = {2, 0, -1};
	config33a.outerLocal = 0;
	config33a.innerCellID = {-1, -2, 0};
	config33a.innerLocal = 1;
	config33a.sign = {1,0,0};
	PrintConfig( config33a );
	InterfacePairConfigurationTransformer()( config33a,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config33a );
	cout<<endl;
	InterfacePairConfiguration config33b;
	config33b.orientation = REGULAR;
	config33b.outerCellID = {0, 2, -1};
	config33b.outerLocal = 0;
	config33b.innerCellID = {-2, -1, 0};
	config33b.innerLocal = 1;
	config33b.sign = {0,1,0};
	PrintConfig( config33b );
	InterfacePairConfigurationTransformer()( config33b,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config33b );
	cout<<endl;
	InterfacePairConfiguration config33c;
	config33c.orientation = SHIFT;
	config33c.outerCellID = {0, 0, -1};
	config33c.outerLocal = 0;
	config33c.innerCellID = {-2, -1, 0};
	config33c.innerLocal = 1;
	config33c.sign = {0,0,0};
	PrintConfig( config33c );
	InterfacePairConfigurationTransformer()( config33c,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config33c );
	cout<<endl;
	InterfacePairConfiguration config33d;
	config33d.orientation = REGULAR;
	config33d.outerCellID = {1, 0, -1};
	config33d.outerLocal = 0;
	config33d.innerCellID = {0, -2, 0};
	config33d.innerLocal = 1;
	config33d.sign = {1,0,0};
	PrintConfig( config33d );
	InterfacePairConfigurationTransformer()( config33d,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config33d );
	cout<<endl;
}

template < typename LATTICE >
void ConfigAliasTest3D ()
{
	InterfacePairConfiguration config11;
	config11.orientation = REGULAR;
	config11.outerCellID = {1, -1, -1};
	config11.outerLocal = 0;
	config11.innerCellID = {-1, 1, 0};
	config11.innerLocal = 0;
	config11.sign = {0, 0, 0};
	PrintConfig( config11 );
	InterfacePairConfigurationTransformer()( config11,
		KernelRelation<LATTICE>::ShiftFunc,
		LATTICE::Rotate
	);
	PrintConfig( config11 );
}


int main()
{
try{
//	CellizeTest_Tri();
//	CellizeTest_Si ();
	CellizeTest_FCC ();
//	Outer2InnerPairTest();

	//ReadKernelFunctionTest();

//	ConfigAliasTest2D<Lattice_Triangular> ();
//	ConfigAliasTest2D<Lattice_Si> ();
//	ConfigAliasTest3D<Lattice_Si> ();
//	ConfigAliasTest3D<Lattice_FCC> ();
}
catch ( exception e )
{
	cout<< e.what() << endl;
}
}



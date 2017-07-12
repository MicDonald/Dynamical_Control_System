#include "Potential_HarmonicDihedral.h"
#include "Potential_CosineDihedral.h"
#include "Potential_HarmonicCosineDihedral.h"
#include <iomanip>
#include <iostream>

using namespace std;


//Potential_HarmonicDihedral potential;
//Potential_CosineDihedral potential;
Potential_HarmonicCosineDihedral potential;

array3d rij{0.4, -0.4, 0.9}; // rj-ri
array3d rjk{0.4, 0.38, -0.4}; // rk-rj
array3d rkl{0.1, 0.7, 0.6}; // rl-rk

constexpr double FiniteDist = 1.e-10;
constexpr double InverseFiniteDist = 1.e10;

void test_1stDerivative ()
{
	auto energy = potential.Energy(rij, rjk, rkl);
	cout << energy<<endl;

	auto force = potential.Force(rij, rjk, rkl);
	cout << "potentialPack:" << endl;
	cout << "fi: " << force[0][0] << " " << force[0][1] << " " << force[0][2] << endl;
	cout << "fj: " << -force[0][0]-force[1][0] << " " << -force[0][1]-force[1][1] << " " << -force[0][2]-force[1][2] << endl;
	cout << "fk: " << force[1][0]-force[2][0] << " " << force[1][1]-force[2][1] << " " << force[1][2]-force[2][2] << endl;
	cout << "fl: " << force[2][0] << " " << force[2][1] << " " << force[2][2] << endl;
	cout << endl;

	force = potential.Force(rij, rjk, rkl, FiniteDifference_t());
	cout << "potentialPack (finite difference):" << endl;
	cout << "fi: " << force[0][0] << " " << force[0][1] << " " << force[0][2] << endl;
	cout << "fj: " << -force[0][0]-force[1][0] << " " << -force[0][1]-force[1][1] << " " << -force[0][2]-force[1][2] << endl;
	cout << "fk: " << force[1][0]-force[2][0] << " " << force[1][1]-force[2][1] << " " << force[1][2]-force[2][2] << endl;
	cout << "fl: " << force[2][0] << " " << force[2][1] << " " << force[2][2] << endl;
	cout << endl;

	cout << "finite difference:" << endl;
	cout << "fi: ";
	rij[0] -= FiniteDist; // move ri
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[0] += FiniteDist;
	rij[1] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[1] += FiniteDist;
	rij[2] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rij[2] += FiniteDist;

	cout << "fj: ";
	rij[0] += FiniteDist; // move rj
	rjk[0] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[0] -= FiniteDist;
	rjk[0] += FiniteDist;
	rij[1] += FiniteDist;
	rjk[1] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rij[1] -= FiniteDist;
	rjk[1] += FiniteDist;
	rij[2] += FiniteDist;
	rjk[2] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rij[2] -= FiniteDist;
	rjk[2] += FiniteDist;

	cout << "fk: ";
	rjk[0] += FiniteDist; // move rk
	rkl[0] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rjk[0] -= FiniteDist;
	rkl[0] += FiniteDist;
	rjk[1] += FiniteDist;
	rkl[1] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rjk[1] -= FiniteDist;
	rkl[1] += FiniteDist;
	rjk[2] += FiniteDist;
	rkl[2] -= FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rjk[2] -= FiniteDist;
	rkl[2] += FiniteDist;

	cout << "fl: ";
	rkl[0] += FiniteDist; // move rl
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rkl[0] -= FiniteDist;
	rkl[1] += FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << " ";
	rkl[1] -= FiniteDist;
	rkl[2] += FiniteDist;
	cout << -(potential.Energy(rij, rjk, rkl)-energy) * InverseFiniteDist << endl;
	rkl[2] -= FiniteDist;
}


void test_2ndDerivative ()
{
	cout << "potentialPack:"<<endl;
	auto hessian = potential.Hessian(rij, rjk, rkl);
	cout << "Hij\tHik\tHil:" << endl;
	cout
	<< hessian[0][0][0] << " " << hessian[0][0][1] << " " << hessian[0][0][2] << "\t"
	<< hessian[1][0][0] << " " << hessian[1][0][1] << " " << hessian[1][0][2] << "\t"
	<< hessian[2][0][0] << " " << hessian[2][0][1] << " " << hessian[2][0][2] << "\n"
	<< hessian[0][1][0] << " " << hessian[0][1][1] << " " << hessian[0][1][2] << "\t"
	<< hessian[1][1][0] << " " << hessian[1][1][1] << " " << hessian[1][1][2] << "\t"
	<< hessian[2][1][0] << " " << hessian[2][1][1] << " " << hessian[2][1][2] << "\n"
	<< hessian[0][2][0] << " " << hessian[0][2][1] << " " << hessian[0][2][2] << "\t"
	<< hessian[1][2][0] << " " << hessian[1][2][1] << " " << hessian[1][2][2] << "\t"
	<< hessian[2][2][0] << " " << hessian[2][2][1] << " " << hessian[2][2][2] << "\n";
	cout << "Hjk\tHjl:" << endl;
	cout
	<< hessian[3][0][0] << " " << hessian[3][0][1] << " " << hessian[3][0][2] << "\t"
	<< hessian[4][0][0] << " " << hessian[4][0][1] << " " << hessian[4][0][2] << "\n"
	<< hessian[3][1][0] << " " << hessian[3][1][1] << " " << hessian[3][1][2] << "\t"
	<< hessian[4][1][0] << " " << hessian[4][1][1] << " " << hessian[4][1][2] << "\n"
	<< hessian[3][2][0] << " " << hessian[3][2][1] << " " << hessian[3][2][2] << "\t"
	<< hessian[4][2][0] << " " << hessian[4][2][1] << " " << hessian[4][2][2] << "\n";
	cout << "Hkl:" << endl;
	cout
	<< hessian[5][0][0] << " " << hessian[5][0][1] << " " << hessian[5][0][2] << "\n"
	<< hessian[5][1][0] << " " << hessian[5][1][1] << " " << hessian[5][1][2] << "\n"
	<< hessian[5][2][0] << " " << hessian[5][2][1] << " " << hessian[5][2][2] << "\n" << endl;

	hessian = potential.Hessian(rij, rjk, rkl, FiniteDifference_t());
	cout << "Hij\tHik\tHil:" << endl;
	cout
	<< hessian[0][0][0] << " " << hessian[0][0][1] << " " << hessian[0][0][2] << "\t"
	<< hessian[1][0][0] << " " << hessian[1][0][1] << " " << hessian[1][0][2] << "\t"
	<< hessian[2][0][0] << " " << hessian[2][0][1] << " " << hessian[2][0][2] << "\n"
	<< hessian[0][1][0] << " " << hessian[0][1][1] << " " << hessian[0][1][2] << "\t"
	<< hessian[1][1][0] << " " << hessian[1][1][1] << " " << hessian[1][1][2] << "\t"
	<< hessian[2][1][0] << " " << hessian[2][1][1] << " " << hessian[2][1][2] << "\n"
	<< hessian[0][2][0] << " " << hessian[0][2][1] << " " << hessian[0][2][2] << "\t"
	<< hessian[1][2][0] << " " << hessian[1][2][1] << " " << hessian[1][2][2] << "\t"
	<< hessian[2][2][0] << " " << hessian[2][2][1] << " " << hessian[2][2][2] << "\n";
	cout << "Hjk\tHjl:" << endl;
	cout
	<< hessian[3][0][0] << " " << hessian[3][0][1] << " " << hessian[3][0][2] << "\t"
	<< hessian[4][0][0] << " " << hessian[4][0][1] << " " << hessian[4][0][2] << "\n"
	<< hessian[3][1][0] << " " << hessian[3][1][1] << " " << hessian[3][1][2] << "\t"
	<< hessian[4][1][0] << " " << hessian[4][1][1] << " " << hessian[4][1][2] << "\n"
	<< hessian[3][2][0] << " " << hessian[3][2][1] << " " << hessian[3][2][2] << "\t"
	<< hessian[4][2][0] << " " << hessian[4][2][1] << " " << hessian[4][2][2] << "\n";
	cout << "Hkl:" << endl;
	cout
	<< hessian[5][0][0] << " " << hessian[5][0][1] << " " << hessian[5][0][2] << "\n"
	<< hessian[5][1][0] << " " << hessian[5][1][1] << " " << hessian[5][1][2] << "\n"
	<< hessian[5][2][0] << " " << hessian[5][2][1] << " " << hessian[5][2][2] << "\n" << endl;

	cout << "finite difference:" << endl;
	auto force = potential.Force(rij, rjk, rkl);
	auto forceJ = Scale( Addition( force[0], force[1] ), -1.);
	auto forceK = Addition( force[1], Scale(force[2], -1.) );
	rij[0] -= FiniteDist; // move ri
	auto dForce = potential.Force(rij, rjk, rkl);
	auto dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	auto dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout << "Hii\tHij\tHik\tHil\n"
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rij[0] += FiniteDist;
	rij[1] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rij[1] += FiniteDist;
	rij[2] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rij[2] += FiniteDist;

	rij[0] += FiniteDist; // move rj
	rjk[0] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout << "Hji\tHjj\tHjk\tHjl\n"
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rij[0] -= FiniteDist;
	rjk[0] += FiniteDist;
	rij[1] += FiniteDist;
	rjk[1] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rij[1] -= FiniteDist;
	rjk[1] += FiniteDist;
	rij[2] += FiniteDist;
	rjk[2] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rij[2] -= FiniteDist;
	rjk[2] += FiniteDist;

	rjk[0] += FiniteDist; // move rk
	rkl[0] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout << "Hki\tHkj\tHkk\tHkl:\n"
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rjk[0] -= FiniteDist;
	rkl[0] += FiniteDist;
	rjk[1] += FiniteDist;
	rkl[1] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rjk[1] -= FiniteDist;
	rkl[1] += FiniteDist;
	rjk[2] += FiniteDist;
	rkl[2] -= FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rjk[2] -= FiniteDist;
	rkl[2] += FiniteDist;

	rkl[0] += FiniteDist; // move rl
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout << "Hli\tHlj\tHlk\tHll:\n"
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rkl[0] -= FiniteDist;
	rkl[1] += FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rkl[1] -= FiniteDist;
	rkl[2] += FiniteDist;
	dForce = potential.Force(rij, rjk, rkl);
	dForceJ = Scale( Addition( dForce[0], dForce[1] ), -1.);
	dForceK = Addition( dForce[1], Scale(dForce[2], -1.) );
	cout
	<< -(dForce[0][0]-force[0][0])/FiniteDist << " "
	<< -(dForce[0][1]-force[0][1])/FiniteDist << " "
	<< -(dForce[0][2]-force[0][2])/FiniteDist << "\t"
	<< -(dForceJ[0]-forceJ[0])/FiniteDist << " "
	<< -(dForceJ[1]-forceJ[1])/FiniteDist << " "
	<< -(dForceJ[2]-forceJ[2])/FiniteDist << "\t"
	<< -(dForceK[0]-forceK[0])/FiniteDist << " "
	<< -(dForceK[1]-forceK[1])/FiniteDist << " "
	<< -(dForceK[2]-forceK[2])/FiniteDist << "\t"
	<< -(dForce[2][0]-force[2][0])/FiniteDist << " "
	<< -(dForce[2][1]-force[2][1])/FiniteDist << " "
	<< -(dForce[2][2]-force[2][2])/FiniteDist << "\n";
	rkl[2] -= FiniteDist;
}


int main ()
{
	test_1stDerivative();
cout << endl;
	test_2ndDerivative();
}


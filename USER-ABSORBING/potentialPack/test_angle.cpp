#include "Potential_HarmonicAngle.h"
#include "Potential_HarmonicCosineAngle.h"
#include "Potential_CosineAngle.h"
#include "Potential_QuarticAngle.h"
#include "Potential_StillingerWeber.h"
#include "Potential_ScreenedHarmonic.h"
#include "Potential_Compass.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//Potential_HarmonicAngle potential(0.5);
//Potential_HarmonicCosineAngle potential;
//Potential_CosineAngle potential;
//Potential_QuarticAngle potential;
//Potential_StillingerWeber potential;
//Potential_ScreenedHarmonicAngle potential;
Potential_Compass potential;

array3d rij{-sqrt(3.)*0.5, -0.5, 0.}; // rj-ri
array3d rik{sqrt(3.)*0.5, -0.5, 0.}; // rk-ri
constexpr double finiteDist = 1.e-10;


void test_1stDerivative ()
{
	auto force = potential.Force(rij, rik);
	cout << "potentialPack:" << endl;
	cout << "fj: " << force[0][0] << " " <<force[0][1] << " " << force[0][2] << "\n"
	<< "fk: " << force[1][0] << " " <<force[1][1] << " " << force[1][2] << endl;

	force = potential.Force(rij, rik, FiniteDifference_t());
	cout << "potentialPack (finite difference):" << endl;
	cout << "fj: " << force[0][0] << " " <<force[0][1] << " " << force[0][2] << "\n"
	<< "fk: " << force[1][0] << " " <<force[1][1] << " " << force[1][2] << endl;

	cout << "finite difference:" << endl;
	double energy = potential.Energy(rij, rik);
	double e;
	cout << "fj: ";
	rij[0] += finiteDist; // move rj
	e = potential.Energy(rij, rik);
	cout << -(e-energy)/finiteDist << " ";
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << -(e-energy)/finiteDist << " ";
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << -(e-energy)/finiteDist << "\n";
	rij[2] -= finiteDist;
	cout << "fk: ";
	rik[0] += finiteDist; // move rk
	e = potential.Energy(rij, rik);
	cout << -(e-energy)/finiteDist << " ";
	rik[0] -= finiteDist;
	rik[1] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << -(e-energy)/finiteDist << " ";
	rik[1] -= finiteDist;
	rik[2] += finiteDist;
	e = potential.Energy(rij, rik);
	cout << -(e-energy)/finiteDist << endl;
	rik[2] -= finiteDist;
}

void test_2ndDerivative ()
{
	auto hessian = potential.Hessian(rij, rik);
	cout << "potentialPack:" << endl;
	cout << "Hij:\tHik\tHjk" << endl
	<< hessian[0][0][0] << " " << hessian[0][0][1] << " " << hessian[0][0][2] << "\t"
	<< hessian[1][0][0] << " " << hessian[1][0][1] << " " << hessian[1][0][2] << "\t"
	<< hessian[2][0][0] << " " << hessian[2][0][1] << " " << hessian[2][0][2] << "\n"
	<< hessian[0][1][0] << " " << hessian[0][1][1] << " " << hessian[0][1][2] << "\t"
	<< hessian[1][1][0] << " " << hessian[1][1][1] << " " << hessian[1][1][2] << "\t"
	<< hessian[2][1][0] << " " << hessian[2][1][1] << " " << hessian[2][1][2] << "\n"
	<< hessian[0][2][0] << " " << hessian[0][2][1] << " " << hessian[0][2][2] << "\t"
	<< hessian[1][2][0] << " " << hessian[1][2][1] << " " << hessian[1][2][2] << "\t"
	<< hessian[2][2][0] << " " << hessian[2][2][1] << " " << hessian[2][2][2] << "\n"<< endl;

	hessian = potential.Hessian(rij, rik, FiniteDifference_t());
	cout << "potentialPack (finite difference):" << endl;
	cout << "Hij:\tHik\tHjk" << endl
	<< hessian[0][0][0] << " " << hessian[0][0][1] << " " << hessian[0][0][2] << "\t"
	<< hessian[1][0][0] << " " << hessian[1][0][1] << " " << hessian[1][0][2] << "\t"
	<< hessian[2][0][0] << " " << hessian[2][0][1] << " " << hessian[2][0][2] << "\n"
	<< hessian[0][1][0] << " " << hessian[0][1][1] << " " << hessian[0][1][2] << "\t"
	<< hessian[1][1][0] << " " << hessian[1][1][1] << " " << hessian[1][1][2] << "\t"
	<< hessian[2][1][0] << " " << hessian[2][1][1] << " " << hessian[2][1][2] << "\n"
	<< hessian[0][2][0] << " " << hessian[0][2][1] << " " << hessian[0][2][2] << "\t"
	<< hessian[1][2][0] << " " << hessian[1][2][1] << " " << hessian[1][2][2] << "\t"
	<< hessian[2][2][0] << " " << hessian[2][2][1] << " " << hessian[2][2][2] << "\n"<< endl;

	cout << "finite difference:" << endl;
	auto force = potential.Force(rij, rik);
	auto forceI = Scale( Addition(force[0], force[1]), -1. );
	rij[0] -= finiteDist; // move ri
	rik[0] -= finiteDist;
	auto dForce = potential.Force(rij, rik);
	auto dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout << "Hii:\tHij\tHik\n"
	<<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[0] += finiteDist;
	rik[0] += finiteDist;
	rij[1] -= finiteDist;
	rik[1] -= finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout
	<<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[1] += finiteDist;
	rik[1] += finiteDist;
	rij[2] -= finiteDist;
	rik[2] -= finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout
	<<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[2] += finiteDist;
	rik[2] += finiteDist;

	rij[0] += finiteDist; // move rj
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout << "Hji\tHjj\tHjk:\n"
	<<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rij[2] -= finiteDist;

	rik[0] += finiteDist; // move rk
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	"Hki\tHkj\tHkk\n"
	<<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rik[0] -= finiteDist;
	rik[1] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << "\n";
	rik[1] -= finiteDist;
	rik[2] += finiteDist;
	dForce = potential.Force(rij, rik);
	dForceI = Scale( Addition(dForce[0], dForce[1]), -1.);
	cout <<	-(dForceI[0]-forceI[0])/finiteDist << " "
	<<	-(dForceI[1]-forceI[1])/finiteDist << " "
	<<	-(dForceI[2]-forceI[2])/finiteDist << "\t"
	<<	-(dForce[0][0]-force[0][0])/finiteDist << " "
	<<	-(dForce[0][1]-force[0][1])/finiteDist << " "
	<<	-(dForce[0][2]-force[0][2])/finiteDist << "\t"
	<<	-(dForce[1][0]-force[1][0])/finiteDist << " "
	<<	-(dForce[1][1]-force[1][1])/finiteDist << " "
	<<	-(dForce[1][2]-force[1][2])/finiteDist << endl;
	rik[2] -= finiteDist;
}

int main ()
{
cout<<setprecision(10);
	potential._1stDebug(rij, rik);
	test_1stDerivative();
cout << endl;
	potential._2ndDebug(rij, rik);
	test_2ndDerivative();
}


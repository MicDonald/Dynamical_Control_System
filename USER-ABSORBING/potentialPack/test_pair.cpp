#include "Potential_HarmonicBond.h"
#include "Potential_LennardJones.h"
#include "Potential_Morse.h"
#include "Potential_Buckingham.h"
#include "Potential_StillingerWeber.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

Potential_HarmonicBond potential;
//Potential_LennardJones potential;
//Potential_Morse potential;
//Potential_Buckingham potential;
//Potential_StillingerWeber potential;

array3d rij{-1.4, 1.3, -0.1}; // rj-ri
constexpr double finiteDist = 1.e-10;


void test_1stDerivative ()
{
	auto force = potential.Force(rij);
	cout << "potentialPack:" << endl;
	cout << "fi: " << force[0] << " " << force[1] << " " << force[2] << endl;
	cout << "fj: " << -force[0] << " " << -force[1] << " " << -force[2] << endl;

	force = potential.Force(rij, FiniteDifference_t());
	cout << "potentialPack (finite difference):" << endl;
	cout << "fi: " << force[0] << " " << force[1] << " " << force[2] << endl;
	cout << "fj: " << -force[0] << " " << -force[1] << " " << -force[2] << endl;

	double energy = potential.Energy(rij);
	double e;
	cout << "finite difference:" << endl;
	cout << "fi: ";
	rij[0] -= finiteDist; // move ri
	e = potential.Energy(rij);
	cout << -(e-energy)/finiteDist << " ";
	rij[0] += finiteDist;
	rij[1] -= finiteDist;
	e = potential.Energy(rij);
	cout << -(e-energy)/finiteDist << " ";
	rij[1] += finiteDist;
	rij[2] -= finiteDist;
	e = potential.Energy(rij);
	cout << -(e-energy)/finiteDist << endl;
	rij[2] += finiteDist;
}

void test_2ndDerivative ()
{
	auto hessian = potential.Hessian(rij);
	cout << "potentialPack:" << endl;
	cout << "Hij: " << endl
	<< hessian[0][0] << " " << hessian[0][1] << " " << hessian[0][2] << "\n"
	<< hessian[1][0] << " " << hessian[1][1] << " " << hessian[1][2] << "\n"
	<< hessian[2][0] << " " << hessian[2][1] << " " << hessian[2][2] << endl<<endl;

	hessian = potential.Hessian(rij, FiniteDifference_t());
	cout << "potentialPack (finite difference):" << endl;
	cout << "Hij: " << endl
	<< hessian[0][0] << " " << hessian[0][1] << " " << hessian[0][2] << "\n"
	<< hessian[1][0] << " " << hessian[1][1] << " " << hessian[1][2] << "\n"
	<< hessian[2][0] << " " << hessian[2][1] << " " << hessian[2][2] << endl<<endl;

	// compare to drjdri
	cout << "finite difference:\nHij:"<<endl;
	auto forceI = potential.Force(rij);
	rij[0] += finiteDist; // move rj
	auto dForceI = potential.Force(rij);
	cout << -(dForceI[0]-forceI[0])/finiteDist << " "
		<< -(dForceI[1]-forceI[1])/finiteDist << " "
		<< -(dForceI[2]-forceI[2])/finiteDist << endl;
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	dForceI = potential.Force(rij);
	cout << -(dForceI[0]-forceI[0])/finiteDist << " "
		<< -(dForceI[1]-forceI[1])/finiteDist << " "
		<< -(dForceI[2]-forceI[2])/finiteDist << endl;
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	dForceI = potential.Force(rij);
	cout << -(dForceI[0]-forceI[0])/finiteDist << " "
		<< -(dForceI[1]-forceI[1])/finiteDist << " "
		<< -(dForceI[2]-forceI[2])/finiteDist << endl;
	rij[2] -= finiteDist;
}

int main ()
{
cout << setprecision(10);
	test_1stDerivative();
cout << endl;
	test_2ndDerivative();
}


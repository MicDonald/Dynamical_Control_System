#include "Potential.h"
#include "Potential_Tersoff.h"
#include <iomanip>
#include <iostream>
#include <cmath>

using namespace std;

Potential_Tersoff potential;
array3d rij{-1.5, 1.3, -0.1};
std::vector<array3d> rik{
	{1.4, 1.4, 0.1},
	{0.2, -1.3, -1.5},
	{-0.2, 1.5, 1.3}
};
constexpr double finiteDist = 1.e-10;

void
test_Tersoff_1stDerivative ()
{
	double energy = potential.Energy(rij, rik);

	auto forces = potential.Force(rij, rik);
	for ( auto& f : forces )
		cout << f[0] << " " << f[1] << " " << f[2] << endl;
	cout << endl;

	forces = potential.Force(rij, rik, FiniteDifference_t());
	for ( auto& f : forces )
		cout << f[0] << " " << f[1] << " " << f[2] << endl;
	cout << endl;

	rij[0] += finiteDist;
	cout << (potential.Energy(rij, rik) - energy)/finiteDist << " ";
	rij[0] -= finiteDist;
	rij[1] += finiteDist;
	cout << (potential.Energy(rij, rik) - energy)/finiteDist << " ";
	rij[1] -= finiteDist;
	rij[2] += finiteDist;
	cout << (potential.Energy(rij, rik) - energy)/finiteDist << endl;
	rij[2] -= finiteDist;
	for ( unsigned i=0; i<rik.size(); ++i )
	{
		rik[i][0] += finiteDist;
		cout << (potential.Energy(rij, rik) - energy)/finiteDist << " ";
		rik[i][0] -= finiteDist;
		rik[i][1] += finiteDist;
		cout << (potential.Energy(rij, rik) - energy)/finiteDist << " ";
		rik[i][1] -= finiteDist;
		rik[i][2] += finiteDist;
		cout << (potential.Energy(rij, rik) - energy)/finiteDist << endl;
		rik[i][2] -= finiteDist;
	}
}


void
test_Tersoff_2ndDerivative ()
{
	;
}


int
main ()
{
cout << setprecision(10);
	test_Tersoff_1stDerivative();
cout << endl;
	test_Tersoff_2ndDerivative();
}


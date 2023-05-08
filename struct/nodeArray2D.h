# include <vector>
using namespace std;

# ifndef nodeArray_H
# define nodeArray_H
struct nodeArray
{
	double* coordinate;
	vector<int> neighbourNode {};
	vector<int> neighbourElement {};
	double** conservedVariable; /* O : Ex ; 1 : Ey ; 2 : Hz */
	double* fluxResidual;
	double* dissipation;
	double** firstDerivative;
	double nodeArea;
  	char boundary;
	complexNumber Ex;
	complexNumber Ey;
	complexNumber Ez;
	complexNumber Hx;
	complexNumber Hy;
	complexNumber Hz;
};
# endif

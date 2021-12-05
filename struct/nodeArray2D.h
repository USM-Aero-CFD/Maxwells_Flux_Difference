# ifndef nodeArray_H
# define nodeArray_H
struct nodeArray
{
	double* coordinate;
	int nodeSharedNumber;
	int nodeLinkNumber;
	double** conservedVariable; /* O : Ex ; 1 : Ey ; 2 : Hz */
	double* fluxResidual;
	double nodeArea;
  	char boundary;
	complexNumber Hr;
	complexNumber Hphi;
	complexNumber Ez;
};
# endif

# ifndef nodeArray_H
# define nodeArray_H
struct nodeArray
{
	double* coordinate;
	double* fluxResidual;
	double** conservedVariable; /* O : Ex ; 1 : Ey ; 2 : Hz */
	char boundary;
	double volume;

	complexNumber ERIncident;
	complexNumber ERScattered;

	complexNumber EThetaIncident;
	complexNumber EThetaScattered;

	complexNumber EPhiIncident;
	complexNumber EPhiScattered;

	complexNumber HRIncident;
	complexNumber HRScattered;

	complexNumber HThetaIncident;
	complexNumber HThetaScattered;

	complexNumber HPhiIncident;
	complexNumber HPhiScattered;
};
# endif

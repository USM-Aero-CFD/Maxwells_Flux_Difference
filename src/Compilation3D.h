# include "MeshConstruction3D.h"
# include "Initialization3D.h"
# include "Computation3D.h"
using namespace std;

# ifndef Compilation3D_H
# define Compilation3D_H
class Compilation3D: public Computation3D
{
	public:
		double L1Calculation() const;
		double L2Calculation() const;
		double LInfinityCalculation() const;
		double mean(double* const &, const int &) const;
		double sumOfSquare(double* const &, double* const &, const int &) const;
		void leastSquareCalculation(double* const &, double* const &, const int &);
		void gridIteration(const bool &, const int &, const char &,
			 				const char &, const char &,
							const string &, const string &, const string &, const string &);
		Compilation3D();
		~Compilation3D();
	private:
		string gMshFile;
		string tetrahedronFile;
		string nodeFile;
		string crossSectionFile;
		double logL1;
		double logL2;
		double logLInfinity;
};
# endif

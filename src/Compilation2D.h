# include "MeshConstruction2D.h"
# include "Initialization2D.h"
# include "Computation2D.h"
using namespace std;

# ifndef Compilation2D_H
# define Compilation2D_H
class Compilation2D: public Computation2D
{
	public:
		double L1Calculation() const;
		double L2Calculation() const;
		double LInfinityCalculation() const;
		double mean(double* const &, const int &) const;
		double sumOfSquare(double* const &, double* const &, const int &) const;
		void leastSquareCalculation(double* const &, double* const &, const int &);
		void gridIteration(const bool &, const int &, const char &,
			 				const char &, const char &, const double &, const double &,
							const string &, const string &, const string &);
		Compilation2D();
		~Compilation2D();
	private:
		string gMshFile;
		string elementFile;
		string nodeFile;
		bool randomization;
		double logL1;
		double logL2;
		double logLInfinity;
};
# endif

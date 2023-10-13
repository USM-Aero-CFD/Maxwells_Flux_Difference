# define _USE_MATH_DEFINES
# include <string>
# include <cmath>
# include <iomanip>
# include <fstream>
# include <iostream>
# include "Initialization2D.h"
using namespace std;

# ifndef Computation2D_H
# define Computation2D_H
class Computation2D: public Initialization2D
{
	public:
		complexNumber complexAddition(const complexNumber &, const complexNumber &) const;
		complexNumber complexComplexMultiplication(const complexNumber &, const complexNumber &) const;
		complexNumber scalarComplexMultiplication(const float &, const complexNumber &) const;
		complexNumber imaginaryCoefficient(const int &) const;
		void TMmodeSolution();
		void TEmodeSolution();
		void spatialSolution();
		void timeDependentSolution(const int &, const int &, const double &);
		void fieldInitialization();
		double larger(const double &, const double &) const;
		double smaller(const double &, const double &) const;
		int KroneckerDelta(const int & value1, const int & value2) const;
		double dotProduct(double* const &, double* const &) const;
		void constructNodeNeighbour();
		void constructElementNeighbour();
		void calculateEdgeTangent();
		void calculateMedianDualNormal();
		void constructInwardNormal();
		void medianCellArea();
		void calculateCentroid();
		double determinant(double** const &) const;
		void inverseMatrix(double** &) const;
		double globalTimeStep();
		void GaussElimination(float** const &, float** &, const int &) const;
		void matrixVectorMultiplication(float** const &, float* const &, float* &, const int &) const;
		double LagrangeInterpolation(const double &, const double &, const int &, const int &) const;
		void constructDistributionMatrix();
		void calculateFiniteVolume(const int &);
		void TMFluxDifference(const int &, const int &);
		void TEFluxDifference(const int &, const int &);
		void calculateGradient(const int &);
		void calculateFluxDifference(const int &);
		void calculateFiniteVolumeBoundaryFlux(const int &);
        void calculateFluxDifferenceBoundaryFlux(const int &);
		double toleranceCalculation(double** const &) const;
		void boundaryCondition(const int &);
		void finiteVolumeNodalUpdate();
		void fluxDifferenceNodalUpdate(const int &);
		void errorsCalculation();
		void printResults() const;
		void intervalResults(bool &, bool &, bool &, bool &) const;
		void initializeArray();
		void timeIterations(const char &, const char &, const double &, const double &);
		~Computation2D();
	protected:
		char TMTEMode;
		char method;
		double time;
		double timeDelta;
		double timeLast;
		int timeNumber;
		double permittivity = 1.0;
		double permeability = 1.0;
		double propagationCoefficient = M_PI;
		double angularFrequency = propagationCoefficient / sqrt(permittivity * permeability);
		double speed = 1.0 / sqrt(permittivity * permeability);
		double impedance = sqrt(permeability / permittivity);
		float** amplificationFactor;
};
# endif

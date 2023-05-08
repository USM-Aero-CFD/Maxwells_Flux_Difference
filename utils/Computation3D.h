# define _USE_MATH_DEFINES
# include <string>
# include <cmath>
# include <ctime>
# include <iomanip>
# include <fstream>
# include <iostream>
# include "Initialization3D.h"
using namespace std;

# ifndef Computation3D_H
# define Computation3D_H
class Computation3D: public Initialization3D
{
	public:
		complexNumber complexAddition(const complexNumber &, const complexNumber &) const;
		complexNumber complexMultiplication(const complexNumber &, const complexNumber &) const;
		complexNumber complexConjugate(const complexNumber &) const;
		complexNumber scalarMultiplication(const float &, const complexNumber &) const;
		float radius(const float &, const float &, const float &) const;
		float theta(const float &, const float &, const float &) const;
		float phi(const float &, const float &, const float &) const;
		void TMmodeSolution();
		void TEmodeSolution();
		void spatialSolution();
		void timeDependentSolution(const int &, const int &, const double &) const;
		void fieldInitialization() const;
		double larger(const double &, const double &) const;
		double smaller(const double &, const double &) const;
		int KroneckerDelta(const int &, const int &) const;
		double dotProduct(double* const &, double* const &) const;
		double crossProduct(double* const &, double* const &, const int &) const;
		void constructNodeNeighbour();
		void constructTetrahedronNeighbour();
		void constructMedianDualNormal();
		void calculateInwardNormal();
		void calculateInwardNormalFromMedianDualNormal();
		void calculateSubVolume();
		void medianCellVolume();
		void calculateCentroid();
		void findCrossSection();
		double globalTimeStep();
		void GaussElimination(float** const &, float** &, const int &) const;
		void matrixVectorMultiplication(float** const &, float* const &, float* &, const int &) const;
		void findTetrahedronOrientation();
		double LagrangeInterpolationFunction(double* const &, const int &, const int &) const;
		void constructDistributionMatrix();
		void calculateFiniteVolume(const int &);
		void fluxDifference(const int &, const int &);
		void calculateGradient(const int &);
		void calculateFluxDifference(const int &);
		void boundaryFluxResidual(double* &, double* const &, const double &, const double &, const double &) const;
		void calculateFiniteVolumeBoundaryFlux(const int &);
		void calculateFluxDifferenceBoundaryFlux(const int &);
		void boundaryCondition(const int &);
		void finiteVolumeNodalUpdate();
		void fluxDifferenceNodalUpdate();
		void interpolateCrossSection();
		void errorsCalculation();
		void printResults() const;
		void intervalResults(bool &, bool &, bool &, bool &) const;
		void printCrossSection() const;
		void initializeArray();
		void timeComputations(const char &, const char &, const double &, const double &);
		~Computation3D();
	protected:
		char TMTEMode;
		char method;
		double permittivity = 1.0;
		double permeability = 1.0;
		int timeNumber;
		double time;
		double timeDelta;
		double timeLast;
		double E0 = 1.0;
		double H0 = 1.0;
		bool randomization;
		double speed = 1.0 / sqrt(permittivity * permeability);
		double propagationCoefficient = M_PI;
		double angularFrequency = propagationCoefficient / sqrt(permeability * permittivity);
		double impedance = sqrt(permeability / permittivity);
};
# endif

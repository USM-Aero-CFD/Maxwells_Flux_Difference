# include "Compilation3D.h"

Compilation3D::Compilation3D()
{

}

Compilation3D::~Compilation3D()
{
	cout << "Destruct Compilation3D";
}

// double Compilation3D::L1Calculation() const
// {
// 	double sum = 0.0;
// 	for (int i = 0; i < nodeNumber; i++)
// 		sum += abs(node[i].conservedVariable[2][3]);
// 	return sum / nodeNumber;
// }
//
// double Compilation3D::L2Calculation() const
// {
// 	double sum = 0.0;
// 	for (int i = 0; i < nodeNumber; i++)
// 		sum += pow(node[i].conservedVariable[2][3], 2.0);
// 	return sqrt(sum / nodeNumber);
// }
//
// double Compilation3D::LInfinityCalculation() const
// {
// 	double maximum = - 1.0;
// 	for (int i = 0; i < nodeNumber; i++)
// 		if (abs(node[i].conservedVariable[2][3]) > maximum)
// 			maximum = abs(node[i].conservedVariable[2][3]);
// 	return maximum;
// }
//
// double Compilation3D::mean(double* const & array, const int & size) const
// {
// 	double arraySum = 0.0;
// 	for (int count = 0; count < size; count++)
// 		arraySum += array[count];
//
// 	return static_cast <double> (arraySum / size);
// }
//
// double Compilation3D::sumOfSquare(double* const & xArray, double* const & yArray, const int & size) const
// {
// 	double sum = 0.0;
// 	double xMean = mean(xArray, size);
// 	double yMean = mean(yArray, size);
// 	for (int count = 0; count < size; count++)
// 		sum += (xArray[count] - xMean) * (yArray[count] - yMean);
//
// 	return static_cast <double> (sum);
// }
//
// void Compilation3D::leastSquareCalculation(double* const & xArray, double* const & yArray, const int & size)
// {
// 	double xMean = mean(xArray, size);
// 	double yMean = mean(yArray, size);
//
// 	double fittedGradient = sumOfSquare(xArray, yArray, size) / sumOfSquare(xArray, xArray, size);
// 	double fittedIntercept = yMean - fittedGradient * xMean;
// }

void Compilation3D::gridIteration(const bool & Construct_Grid, const int & Grid_Size, const char & Randomization,
	 								const char & TM_TE_Mode, const char & Method, const double & Time_Last, const double & Time_Delta,
									const string & Gmsh_File, const string & Tetrahedron_File, const string & Node_File, const string & Cross_Section_File)
{
	cout << "Compilation3D" << endl;
	int xNumber, yNumber, zNumber;

	char TMTEMode = TM_TE_Mode;
	char method = Method;
	double timeLast = Time_Last;
	double timeDelta = Time_Delta;

	tetrahedronFile = Tetrahedron_File;
	nodeFile = Node_File;
	crossSectionFile = Cross_Section_File;

	if (Construct_Grid) {
		randomization = Randomization;
		xNumber = Grid_Size;
		yNumber = Grid_Size;
		zNumber = Grid_Size;

		MeshConstruction3D structuredMesh(xNumber, yNumber, zNumber, randomization);
		structuredMesh.generateGrid();
	};

	gMshFile = "Gmsh_3D.msh";
	tetrahedronFile = "Tetrahedron " + to_string(xNumber) + ".txt";
	nodeFile = "Node " + to_string(xNumber) + ".txt";
	crossSectionFile = "CrossSection " + to_string(xNumber) + ".txt";
	gMshFile = Gmsh_File;
	tetrahedronFile = Tetrahedron_File;
	nodeFile = Node_File;
	crossSectionFile = Cross_Section_File;
	Computation3D computation;
	computation.readMesh(gMshFile, tetrahedronFile, nodeFile, crossSectionFile);
	computation.timeComputations(TMTEMode, method, timeLast, timeDelta);

	// // logL1 = log10(L1Calculation());
	// // logL2 = log10(L2Calculation());
	// // logLInfinity = log10(LInfinityCalculation());
	// // cout << logL1 << setw(20) << logL2 << setw(20) << logLInfinity << setw(20) << endl;
}

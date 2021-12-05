# ifndef tetrahedronArray_H
# define tetrahedronArray_H
struct tetrahedronArray
{
	int* globalNode;
	double** medianDualNormal;
	char orientation;
	double** LagrangeInterpolation;
	double** inwardNormal;
	double*** distributionMatrix;
	double volume;
	double* subVolume;
	int boundaryVertex;
	char boundaryType;
};
# endif

# ifndef elementArray_H
# define elementArray_H
struct elementArray
{
	int* globalNode;
	double** edgeTangent;
	double** medianDualNormal;
	double** interpolateVariable;
	double** inwardNormal;
	double*** distributionMatrix;
	double cellArea;
	int boundaryVertex;
	char boundaryType;
};
# endif

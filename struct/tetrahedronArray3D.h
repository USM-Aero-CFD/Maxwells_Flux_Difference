# include <vector>
using namespace std;

# ifndef tetrahedronArray_H
# define tetrahedronArray_H
struct tetrahedronArray
{
	int* globalNode;
	vector<int> neighbourTetrahedron {};
	double** medianDualNormal;
	char orientation;
	double** inwardNormal;
	float* centroid;
	double*** distributionMatrix;
	double** gradient;
	double volume;
	double* subVolume;
	int boundaryVertex;
	char boundaryType;
};
# endif

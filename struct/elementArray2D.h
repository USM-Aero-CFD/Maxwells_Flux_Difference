# include <vector>
using namespace std;

# ifndef elementArray_H
# define elementArray_H
struct elementArray
{
	int* globalNode;
	vector<int> neighbourElement {};
	double** edgeTangent;
	double** medianDualNormal;
	double** inwardNormal;
	float* centroid;
	double*** distributionMatrix;
	double** gradient;
	double cellArea;
	int boundaryVertex;
	char boundaryType;
};
# endif

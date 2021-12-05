# define _USE_MATH_DEFINES
# include <cmath>
# include <fstream>
# include <string>
# include <iostream>
# include <iomanip>
# include "complexNumber.h"
# include "tetrahedronArray3D.h"
# include "nodeArray3D.h"
# include "crossSectionArray3D.h"
using namespace std;

# ifndef Initialization3D_H
# define Initialization3D_H
class Initialization3D
{
	public:
		void initializeGmsh();
		void initializeTetrahedron();
		void initializeNode();
    	void initializeCrossSection();
    	void readMesh(const string &, const string &, const string &, const string &);
		~Initialization3D();
	protected:
		ifstream inputGmsh;
		ifstream inputTetrahedron;
		ifstream inputNode;
		ifstream inputCrossSection;
		int tetrahedronNumber;
		int nodeNumber;
		int crossSectionNumber;
		tetrahedronArray* tetrahedron;
		nodeArray* node;
		crossSectionArray* crossSection;
};
# endif

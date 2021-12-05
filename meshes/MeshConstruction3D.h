# define _USE_MATH_DEFINES
# include <fstream>
# include <iomanip>
# include <vector>
# include <cmath>
# include <iostream>
# include <string>
using std::vector;
using namespace std;

# ifndef MeshConstruction3D_H
# define MeshConstruction3D_H
class MeshConstruction3D
{
	public:
		void constructNode();
		void constructTetrahedron();
		bool verticesOnBoundary(const int &, const int &, const int &, const int &, char &) const;
		void determineBoundaryVertex();
		void mappingCubeToSphere();
		void constructCrossSection();
		void randomizeGrid();
		void calculateTetrahedronVolume();

		void printTetrahedron() const;
		void printNode() const;
		void printCrossSection() const;
		void printGmsh() const;
		void generateGrid();
		MeshConstruction3D(const int &, const int &, const int &, const bool &);
		~MeshConstruction3D();
	private:
		struct globalTetrahedronArray
		{
			int* nodeNumbering;
			double volume;
			int boundaryVertex;
			char boundaryType;
		};
        struct domainNodeArray
		{
			double* coordinate;
			int globalNumbering;
		};
		struct globalNodeArray
		{
			double* coordinate;
			char boundary;
		};
		struct crossSectionArray
        {
            double* coordinate;
            int node;
        };
		double radius;
		double xInitial;
		double xLast;
        double xDelta;
		double yInitial;
		double yLast;
		double yDelta;
        double zInitial;
		double zLast;
		double zDelta;
        vector <domainNodeArray> domainNode;
        vector <globalNodeArray> globalNode;
        vector <globalTetrahedronArray> globalTetrahedron;
        vector <crossSectionArray> crossSection;
        int xNumber;
        int yNumber;
        int zNumber;
		int globalTetrahedronNumber;
        int domainNodeNumber;
		int globalNodeNumber;
		int crossSectionNumber;
		bool randomization;
};
# endif

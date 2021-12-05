# define _USE_MATH_DEFINES
# include <fstream>
# include <iomanip>
# include <vector>
# include <cmath>
# include <iostream>
# include <string>
using std::vector;
using namespace std;


# ifndef MeshConstruction2D_H
# define MeshConstruction2D_H
class MeshConstruction2D
{
	public:
		void constructNode();
		void constructElement();
		void randomizeGrid();
		void calculateCellArea();

		void printElement() const;
		void printNode() const;
		void printGmsh() const;
		void generateGrid();
		MeshConstruction2D(const int &, const int &, const bool &);
		~MeshConstruction2D();
	private:
		struct elementArray
		{
			int* globalNode;
			double cellArea;
			int boundaryVertex;
            char boundaryType;
		};
		struct nodeArray
		{
			double* coordinate;
			char boundary;
		};
		int xNumber;
		int yNumber;
		bool randomization;
		double xInitial;
		double xLast;
        double yInitial;
		double yLast;
        double xDelta;
		double yDelta;
		double wedgeAngle;
		elementArray* element;
		nodeArray* node;
		int** nodeLink;
		int elementNumber;
		int nodeNumber;
};
# endif

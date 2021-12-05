# define _USE_MATH_DEFINES
# include <cmath>
# include <fstream>
# include <string>
# include <iostream>
# include <iomanip>
# include "complexNumber.h"
# include "elementArray2D.h"
# include "nodeArray2D.h"
using namespace std;

# ifndef Initialization2D_H
# define Initialization2D_H
class Initialization2D
{
	public:
		void initializeGmsh();
		void initializeElement();
		void initializeNode();
        void readMesh(const string &, const string &, const string &);
		~Initialization2D();
	protected:
		ifstream inputGmsh;
		ifstream inputElement;
		ifstream inputNode;
		int elementNumber;
		int nodeNumber;
		elementArray* element;
		nodeArray* node;
};
# endif

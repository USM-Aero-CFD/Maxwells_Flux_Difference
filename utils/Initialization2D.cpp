# include "Initialization2D.h"

Initialization2D::~Initialization2D()
{
	cout << "Destruct Initialization2D" << endl;

	for (int e = 0; e < elementNumber; e++)
	{
		delete [] element[e].globalNode; element[e].globalNode = NULL;
	}
	delete [] element; element = NULL;

	for (int i = 0; i < nodeNumber; i++)
	{
		delete [] node[i].coordinate; node[i].coordinate = NULL;
	};
	delete [] node; node = NULL;
}

void Initialization2D::initializeGmsh()
{
	int block, minIndex, maxIndex;
	int dimension, tag, parameter;
	int index;

	string line;
	if ( inputGmsh.is_open() )
	{
		cout << "The Gmsh file exists!!" << endl;

		while( getline(inputGmsh, line) )
		{
			if ( line.find("$Nodes") == static_cast<string::size_type>(0) ) {
				inputGmsh >> block >> nodeNumber >> minIndex >> maxIndex;
				inputGmsh >> dimension >> tag >> parameter >> nodeNumber;

				for (int i = 0; i < nodeNumber; i++)
					inputGmsh >> index;

				for (int i = 0; i < nodeNumber; i++) {
					node[i].coordinate = new double [2];
					inputGmsh >> node[i].coordinate[0] >> node[i].coordinate[1];
				};
			}
			else if ( line.find("$Elements") == static_cast<string::size_type>(0) ) {
				inputGmsh >> block >> elementNumber >> minIndex >> maxIndex;
				inputGmsh >> dimension >> tag >> parameter >> elementNumber;

				for (int e = 0; e < elementNumber; e++) {
					element[e].globalNode = new int [3];
					inputGmsh >> index >> element[e].globalNode[0] >> element[e].globalNode[1] >> element[e].globalNode[2];
				};
			}
		}
	}
	else
		cout << "The Gmsh file does not exist" << endl;
}

void Initialization2D::initializeElement()
{
	element = new elementArray [elementNumber];

	for (int e = 0; e < elementNumber; e++)
	{
		inputElement >> element[e].cellArea >> element[e].boundaryVertex >> element[e].boundaryType;
	};
	inputElement.close();
}

void Initialization2D::initializeNode()
{
	node = new nodeArray [nodeNumber];

	for (int i = 0; i < nodeNumber; i++)
	{
		/* node[i].UVARIABLELEVEL[ki][0] - initial */
		/* node[i].UVARIABLELEVEL[ki][1] - (n - 1) */
		/* node[i].UVARIABLELEVEL[ki][2] - n */
		/* node[i].UVARIABLELEVEL[ki][3] - error */
		/* node[i].UVARIABLELEVEL[ki][4] - exact */

		inputNode >> node[i].boundary;
	};
	inputNode.close();
}

void Initialization2D::readMesh(const string & gmsh_File, const string & element_File, const string & node_File)
{
	inputGmsh.open(gmsh_File.c_str());
	inputElement.open(element_File.c_str());
	inputNode.open(node_File.c_str());
	inputElement >> elementNumber;
	inputNode >> nodeNumber;

	initializeElement();
	initializeNode();
	initializeGmsh();
}

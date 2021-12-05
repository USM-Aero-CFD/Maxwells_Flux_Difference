# include "Initialization3D.h"

Initialization3D::~Initialization3D()
{
	for (int h = 0; h < tetrahedronNumber; h++)
	{
		delete [] tetrahedron[h].globalNode; tetrahedron[h].globalNode = NULL;
	};
	delete [] tetrahedron; tetrahedron = NULL;

	for (int i = 0; i < nodeNumber; i++)
	{
		delete [] node[i].coordinate; node[i].coordinate = NULL;
	};
	delete [] node; node = NULL;

	for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
	{
		delete [] crossSection[crosssection].coordinate; crossSection[crosssection].coordinate = NULL;
	};
	delete [] crossSection; crossSection = NULL;

}

void Initialization3D::initializeGmsh()
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

				for (int i = 0; i < nodeNumber; i++)
					inputGmsh >> node[i].coordinate[0] >> node[i].coordinate[1] >> node[i].coordinate[2];
			}
			else if ( line.find("$Elements") == static_cast<string::size_type>(0) ) {
				inputGmsh >> block >> tetrahedronNumber >> minIndex >> maxIndex;
				inputGmsh >> dimension >> tag >> parameter >> tetrahedronNumber;

				for (int h = 0; h < tetrahedronNumber; h++)
					inputGmsh >> index;

				for (int h = 0; h < tetrahedronNumber; h++)
					inputGmsh >> tetrahedron[h].globalNode[0] >> tetrahedron[h].globalNode[1]
			                  >> tetrahedron[h].globalNode[2] >> tetrahedron[h].globalNode[3];
			}
		}
	}
	else
		cout << "The Gmsh file does not exist" << endl;
}

void Initialization3D::initializeTetrahedron()
{
	for (int h = 0; h < tetrahedronNumber; h++)
	{
		inputTetrahedron >> tetrahedron[h].volume >> tetrahedron[h].boundaryVertex >> tetrahedron[h].boundaryType;

	};
	inputTetrahedron.close();
}

void Initialization3D::initializeNode()
{
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

void Initialization3D::initializeCrossSection()
{
    for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
        inputCrossSection >> crossSection[crosssection].coordinate[0] >> crossSection[crosssection].coordinate[1]
                            >> crossSection[crosssection].coordinate[2] >> crossSection[crosssection].node;
	inputCrossSection.close();
}

void Initialization3D::readMesh(const string & gmsh_File, const string & element_File, const string & node_File, const string & crossSection_File)
{
	inputGmsh.open(gmsh_File.c_str());
	inputTetrahedron.open(element_File.c_str());
	inputNode.open(node_File.c_str());
	inputCrossSection.open(crossSection_File.c_str());
	inputTetrahedron >> tetrahedronNumber;
	inputNode >> nodeNumber;
	inputCrossSection >> crossSectionNumber;

	initializeTetrahedron();
	initializeNode();
	initializeCrossSection();
	initializeGmsh();
}

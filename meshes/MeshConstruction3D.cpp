#include "MeshConstruction3D.h"

MeshConstruction3D::MeshConstruction3D(const int & x_Number, const int & y_Number, const int & z_Number, const bool & Randomization)
{
	radius = 0.2;
	xInitial = - 1.0;
	xLast = 1.0;
	yInitial = - 1.0;
	yLast = 1.0;
	zInitial = - 1.0;
	zLast = 1.0;
	xNumber = x_Number;
	yNumber = y_Number;
	zNumber = z_Number;
	randomization = Randomization;
	xDelta = (xLast - xInitial) / xNumber;
	yDelta = (yLast - yInitial) / yNumber;
	zDelta = (zLast - zInitial) / zNumber;
}

MeshConstruction3D::~MeshConstruction3D()
{
	cout << xNumber << "     Destruct GRID" << endl;
	for (int i = 0; i < domainNodeNumber; i++)
	{
		delete [] domainNode[i].coordinate; domainNode[i].coordinate = NULL;
	};
	domainNode.erase(domainNode.begin());

    for (int i = 0; i < globalNodeNumber; i++)
	{
		delete [] globalNode[i].coordinate; globalNode[i].coordinate = NULL;
	};
	globalNode.erase(globalNode.begin());

	for (int h = 0; h < globalTetrahedronNumber; h++)
	{
		delete [] globalTetrahedron[h].nodeNumbering; globalTetrahedron[h].nodeNumbering = NULL;
	};
	globalTetrahedron.erase(globalTetrahedron.begin());

    for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
	{
		delete [] crossSection[crosssection].coordinate; crossSection[crosssection].coordinate = NULL;
	};
	crossSection.erase(crossSection.begin());
}

void MeshConstruction3D::constructNode()
{
    // globalNode[i].boundary = 'A' is inlet
    // globalNode[i].boundary = 'B' is PEC
    // globalNode[i].boundary = 'C' is outlet
    // globalNode[i].boundary = 'D' are the side boundaries
    // globalNode[i].boundary = 'N' is not on the boundary

    for (int k = 0; k <= zNumber; k++)
        for (int j = 0; j <= yNumber; j++)
            for (int i = 0; i <= xNumber; i++)
            {
                domainNodeNumber = k * (xNumber + 1) * (yNumber + 1) + j * (xNumber + 1) + i;
                domainNode.resize(domainNodeNumber + 1);
                domainNode[domainNodeNumber].coordinate = new double [3];
                domainNode[domainNodeNumber].coordinate[0] = xInitial + i * xDelta;
                domainNode[domainNodeNumber].coordinate[1] = yInitial + j * yDelta;
                domainNode[domainNodeNumber].coordinate[2] = zInitial + k * zDelta;
                domainNode[domainNodeNumber].globalNumbering = - 1;
            };


    globalNodeNumber = 0;
    for (int k = 0; k <= zNumber; k++)
        for (int j = 0; j <= yNumber; j++)
            for (int i = 0; i <= xNumber; i++)
            {
                domainNodeNumber = k * (xNumber + 1) * (yNumber + 1) + j * (xNumber + 1) + i;

                if (! (domainNode[domainNodeNumber].coordinate[0] > (- radius + xDelta / 10.0) && domainNode[domainNodeNumber].coordinate[0] < (radius - xDelta / 10.0)
                       && domainNode[domainNodeNumber].coordinate[1] > (- radius + yDelta / 10.0) && domainNode[domainNodeNumber].coordinate[1] < (radius - yDelta / 10.0)
                       && domainNode[domainNodeNumber].coordinate[2] > (- radius + zDelta / 10.0) && domainNode[domainNodeNumber].coordinate[2] < (radius - zDelta / 10.0)))
                {
                    globalNode.resize(globalNodeNumber + 1);
                    globalNode[globalNodeNumber].coordinate = new double [3];
                    globalNode[globalNodeNumber].coordinate[0] = domainNode[domainNodeNumber].coordinate[0];
                    globalNode[globalNodeNumber].coordinate[1] = domainNode[domainNodeNumber].coordinate[1];
                    globalNode[globalNodeNumber].coordinate[2] = domainNode[domainNodeNumber].coordinate[2];
                    domainNode[domainNodeNumber].globalNumbering = globalNodeNumber;

                    if (k == 0)
                        globalNode[globalNodeNumber].boundary = 'A';
                    else if (i == 0 || i == xNumber || j == 0 || j == yNumber)
                        globalNode[globalNodeNumber].boundary = 'D';
                    else if (k == zNumber)
                        globalNode[globalNodeNumber].boundary = 'C';
                    else if (domainNode[domainNodeNumber].coordinate[0] >= (- radius - xDelta / 10.0) && domainNode[domainNodeNumber].coordinate[0] <= (- radius + xDelta / 10.0)
                            || domainNode[domainNodeNumber].coordinate[0] >= (radius - xDelta / 10.0) && domainNode[domainNodeNumber].coordinate[0] <= (radius + xDelta / 10.0)
                            || domainNode[domainNodeNumber].coordinate[1] >= (- radius - yDelta / 10.0) && domainNode[domainNodeNumber].coordinate[1] <= (- radius + yDelta / 10.0)
                            || domainNode[domainNodeNumber].coordinate[1] >= (radius - yDelta / 10.0) && domainNode[domainNodeNumber].coordinate[1] <= (radius + yDelta / 10.0)
                            || domainNode[domainNodeNumber].coordinate[2] >= (- radius - zDelta / 10.0) && domainNode[domainNodeNumber].coordinate[2] <= (- radius + zDelta / 10.0)
                            || domainNode[domainNodeNumber].coordinate[2] >= (radius - zDelta / 10.0) && domainNode[domainNodeNumber].coordinate[2] <= (radius + zDelta / 10.0))
                        globalNode[globalNodeNumber].boundary = 'B';
                    else
                        globalNode[globalNodeNumber].boundary = 'N';

                    globalNodeNumber++;
                };
            };

    domainNodeNumber = (xNumber + 1) * (yNumber + 1) * (zNumber + 1);
}

void MeshConstruction3D::constructTetrahedron()
{
    globalTetrahedronNumber = 0;
    for (int k = 0; k < zNumber; k++)
        for (int j = 0; j < yNumber; j++)
            for (int i = 0; i < xNumber; i++)
            {
                domainNodeNumber = k * (xNumber + 1) * (yNumber + 1) + j * (xNumber + 1) + i;
                if (!(((xInitial + i * xDelta) + xDelta / 10.0) >= - radius && ((xInitial + i * xDelta) + xDelta / 10.0) <= radius
                    && ((yInitial + j * yDelta) + yDelta / 10.0) >= - radius && ((yInitial + j * yDelta) + yDelta / 10.0) <= radius
                    && ((zInitial + k * zDelta) + zDelta / 10.0) >= - radius && ((zInitial + k * zDelta) + zDelta / 10.0) <= radius))
                {
                    // Local tetrahedron 0
                    globalTetrahedron.resize(globalTetrahedronNumber + 1);
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering = new int [4];
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[0] = domainNode[domainNodeNumber].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[1] = domainNode[domainNodeNumber + (xNumber + 1)].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[2] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[3] = domainNode[domainNodeNumber + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedronNumber++;

                    // Local tetrahedron 1
                    globalTetrahedron.resize(globalTetrahedronNumber + 1);
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering = new int [4];
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[0] = domainNode[domainNodeNumber].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[1] = domainNode[domainNodeNumber + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[2] = domainNode[domainNodeNumber + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[3] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedronNumber++;

                    // Local tetrahedron 2
                    globalTetrahedron.resize(globalTetrahedronNumber + 1);
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering = new int [4];
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[0] = domainNode[domainNodeNumber].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[1] = domainNode[domainNodeNumber + (xNumber + 1)].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[2] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1)].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[3] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedronNumber++;

                    // Local tetrahedron 3
                    globalTetrahedron.resize(globalTetrahedronNumber + 1);
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering = new int [4];
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[0] = domainNode[domainNodeNumber].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[1] = domainNode[domainNodeNumber + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[2] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[3] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + 1].globalNumbering;
                    globalTetrahedronNumber++;

                    // Local tetrahedron 4
                    globalTetrahedron.resize(globalTetrahedronNumber + 1);
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering = new int [4];
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[0] = domainNode[domainNodeNumber].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[1] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1)].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[2] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[3] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1)].globalNumbering;
                    globalTetrahedronNumber++;

                    // Local tetrahedron 5
                    globalTetrahedron.resize(globalTetrahedronNumber + 1);
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering = new int [4];
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[0] = domainNode[domainNodeNumber].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[1] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1)].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[2] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + 1].globalNumbering;
                    globalTetrahedron[globalTetrahedronNumber].nodeNumbering[3] = domainNode[domainNodeNumber + (xNumber + 1) * (yNumber + 1) + (xNumber + 1) + 1].globalNumbering;
                    globalTetrahedronNumber++;
                }
            };

    domainNodeNumber = (xNumber + 1) * (yNumber + 1) * (zNumber + 1);
}

bool MeshConstruction3D::verticesOnBoundary(const int & h, const int & j0, const int & j1, const int & j2, char & boundary_Type) const
{
    bool* criteria;
    criteria = new bool [12];

    criteria[0] = ((abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - (- radius)) <= xDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - yDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - (- radius)) <= xDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - yDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - (- radius)) <= xDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - yDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - zDelta / 10.0) <= radius));

    criteria[1] = ((abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - radius) <= xDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - yDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - radius) <= xDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - yDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - radius) <= xDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - yDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - zDelta / 10.0) <= radius));

    criteria[2] = ((abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - (- radius)) <= yDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - (- radius)) <= yDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - (- radius)) <= yDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - zDelta / 10.0) <= radius));

    criteria[3] = ((abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - radius) <= yDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - radius) <= yDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - zDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - radius) <= yDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] + zDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - zDelta / 10.0) <= radius));

    criteria[4] = ((abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - (- radius)) <= zDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - yDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - (- radius)) <= zDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - yDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - (- radius)) <= zDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - yDelta / 10.0) <= radius));

    criteria[5] = ((abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - radius) <= zDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - yDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - radius) <= zDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - yDelta / 10.0) <= radius)
                && (abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - radius) <= zDelta / 10.0
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] + xDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - xDelta / 10.0) <= radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] + yDelta / 10.0) >= - radius
                        && (globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - yDelta / 10.0) <= radius));

    criteria[6] = (abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - xInitial) <= xDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - xInitial) <= xDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - xInitial) <= xDelta / 10.0);

    criteria[7] = (abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[0] - xLast) <= xDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[0] - xLast) <= xDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[0] - xLast) <= xDelta / 10.0);

    criteria[8] = (abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - yInitial) <= yDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - yInitial) <= yDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - yInitial) <= yDelta / 10.0);

    criteria[9] = (abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[1] - yLast) <= yDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[1] - yLast) <= yDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[1] - yLast) <= yDelta / 10.0);

    criteria[10] = (abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - zInitial) <= zDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - zInitial) <= zDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - zInitial) <= zDelta / 10.0);

    criteria[11] = (abs(globalNode[globalTetrahedron[h].nodeNumbering[j0]].coordinate[2] - zLast) <= zDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j1]].coordinate[2] - zLast) <= zDelta / 10.0
                    && abs(globalNode[globalTetrahedron[h].nodeNumbering[j2]].coordinate[2] - zLast) <= zDelta / 10.0);


    if (criteria[10])
        boundary_Type = 'A';
    else if (criteria[0] || criteria[1] || criteria[2] || criteria[3] || criteria[4] || criteria[5])
        boundary_Type = 'B';
    else if (criteria[11])
        boundary_Type = 'C';
    else if (criteria[6] || criteria[7] || criteria[8] || criteria[9])
        boundary_Type = 'D';
    else
        boundary_Type = 'N';

    return criteria[0] || criteria[1] || criteria[2] || criteria[3] || criteria[4] || criteria[5] ||
            criteria[6] || criteria[7] || criteria[8] || criteria[9] || criteria[10] || criteria[11];

    delete [] criteria; criteria = NULL;
}

void MeshConstruction3D::determineBoundaryVertex()
{
    // globalTetrahedron[h].boundaryType = 'A' is inlet
    // globalTetrahedron[h].boundaryType = 'B' is PEC
    // globalTetrahedron[h].boundaryType = 'C' is outlet
    // globalTetrahedron[h].boundaryType = 'D' are the side boundaries
    // globalTetrahedron[h].boundaryType = 'N' is not on the boundary

    for (int h = 0; h < globalTetrahedronNumber; h++)
    {
        char boundaryType = 'N';
        globalTetrahedron[h].boundaryVertex = - 1;

        if (verticesOnBoundary(h, 1, 2, 3, boundaryType))
        {
            globalTetrahedron[h].boundaryVertex = 0;
        }
        else if (verticesOnBoundary(h, 0, 2, 3, boundaryType))
        {
            globalTetrahedron[h].boundaryVertex = 1;
        }
        else if (verticesOnBoundary(h, 0, 1, 3, boundaryType))
        {
            globalTetrahedron[h].boundaryVertex = 2;
        }
        else if (verticesOnBoundary(h, 0, 1, 2, boundaryType))
        {
            globalTetrahedron[h].boundaryVertex = 3;
        };

        globalTetrahedron[h].boundaryType = boundaryType;
    };

    /*/
    ofstream Test_Tetrahedron;
    Test_Tetrahedron.open("./Output/Test_Tetrahedron.txt");
    Test_Tetrahedron << globalTetrahedronNumber << endl;
    for (int h = 0; h < globalTetrahedronNumber; h++)
        Test_Tetrahedron << setw(10) << globalTetrahedron[h].nodeNumbering[0] << setw(10) << globalTetrahedron[h].nodeNumbering[1]
                            << setw(10) << globalTetrahedron[h].nodeNumbering[2] << setw(10) << globalTetrahedron[h].nodeNumbering[3]
                            << setw(10) << globalTetrahedron[h].boundaryVertex << endl;
    Test_Tetrahedron.close();
    /**/
}

void MeshConstruction3D::mappingCubeToSphere()
{
    double theta;
    double phi;
    double arbitraryRadius;

    //ofstream Testing_Node;
    //Testing_Node.open("./Output/Testing_Node.txt");
    //Testing_Node << globalNodeNumber << endl;

    for (int i = 0; i < globalNodeNumber; i++)
    {
        arbitraryRadius = abs(globalNode[i].coordinate[0]);
        if (abs(globalNode[i].coordinate[1]) > arbitraryRadius)
            arbitraryRadius = abs(globalNode[i].coordinate[1]);

        if (abs(globalNode[i].coordinate[2]) > arbitraryRadius)
            arbitraryRadius = abs(globalNode[i].coordinate[2]);


        if (globalNode[i].coordinate[0] >= 0 && globalNode[i].coordinate[1] >= 0)
            phi = atan2(abs(globalNode[i].coordinate[1]), abs(globalNode[i].coordinate[0]));
        else if (globalNode[i].coordinate[0] < 0 && globalNode[i].coordinate[1] >= 0)
            phi = M_PI - atan2(abs(globalNode[i].coordinate[1]), abs(globalNode[i].coordinate[0]));
        else if (globalNode[i].coordinate[0] < 0 && globalNode[i].coordinate[1] < 0)
            phi = M_PI + atan2(abs(globalNode[i].coordinate[1]), abs(globalNode[i].coordinate[0]));
        else if (globalNode[i].coordinate[0] >= 0 && globalNode[i].coordinate[1] < 0)
            phi = 2.0 * M_PI - atan2(abs(globalNode[i].coordinate[1]), abs(globalNode[i].coordinate[0]));

        if (globalNode[i].coordinate[2] >= 0)
            theta = atan2(sqrt(pow(globalNode[i].coordinate[0], 2) + pow(globalNode[i].coordinate[1], 2)), abs(globalNode[i].coordinate[2]));
        else if (globalNode[i].coordinate[2] < 0)
            theta = M_PI - atan2(sqrt(pow(globalNode[i].coordinate[0], 2) + pow(globalNode[i].coordinate[1], 2)), abs(globalNode[i].coordinate[2]));

        globalNode[i].coordinate[0] = arbitraryRadius * sin(theta) * cos(phi);
        globalNode[i].coordinate[1] = arbitraryRadius * sin(theta) * sin(phi);
        globalNode[i].coordinate[2] = arbitraryRadius * cos(theta);

        //Testing_Node << setw(20) << phi << setw(20) << theta
        //                << setw(20) << globalNode[i].coordinate[0] << setw(20) << globalNode[i].coordinate[1] << setw(20) << globalNode[i].coordinate[2] << endl;
    };
    //Testing_Node.close();
}

void MeshConstruction3D::constructCrossSection()
{
	crossSectionNumber = 0;
	for (int i = 0; i < globalNodeNumber; i++)
        if (abs(globalNode[i].coordinate[1] - 0.0) <= yDelta / 10.0)
        {
            crossSection.resize(crossSectionNumber + 1);
            crossSection[crossSectionNumber].coordinate = new double [3];
            crossSection[crossSectionNumber].coordinate[0] = globalNode[i].coordinate[0];
            crossSection[crossSectionNumber].coordinate[1] = globalNode[i].coordinate[1];
            crossSection[crossSectionNumber].coordinate[2] = globalNode[i].coordinate[2];
            crossSection[crossSectionNumber].node = i;
            crossSectionNumber++;
        };
}

void MeshConstruction3D::randomizeGrid()
{
    for (int i = 0; i < globalNodeNumber; i++)
    {
        if (globalNode[i].boundary == 'N')
        {
            double displacingRadius = 0.35 * (static_cast<double> (rand()) / RAND_MAX) * xDelta;
            double displacingTheta = (static_cast<double> (rand()) / RAND_MAX) * M_PI;
            double displacingPhi = (static_cast<double> (rand()) / RAND_MAX) * 2.0 * M_PI;

            globalNode[i].coordinate[0] = globalNode[i].coordinate[0] + displacingRadius * sin(displacingTheta) * cos(displacingPhi);
            globalNode[i].coordinate[1] = globalNode[i].coordinate[1] + displacingRadius * sin(displacingTheta) * sin(displacingPhi);
            globalNode[i].coordinate[2] = globalNode[i].coordinate[2] + displacingRadius * cos(displacingTheta);
        };
    }
}

void MeshConstruction3D::calculateTetrahedronVolume()
{
	for (int h = 0; h < globalTetrahedronNumber; h++)
	{
		int j = 3;
		globalTetrahedron[h].volume = (1.0 / 6.0) * abs((globalNode[globalTetrahedron[h].nodeNumbering[j]].coordinate[0] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[0])
                                                    * ((globalNode[globalTetrahedron[h].nodeNumbering[(j + 2) % 4]].coordinate[1] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[1])
                                                            * (globalNode[globalTetrahedron[h].nodeNumbering[(j - 1) % 4]].coordinate[2] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[2])
                                                        - (globalNode[globalTetrahedron[h].nodeNumbering[(j - 1) % 4]].coordinate[1] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[1])
                                                            * (globalNode[globalTetrahedron[h].nodeNumbering[(j + 2) % 4]].coordinate[2] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[2]))
                                                + (globalNode[globalTetrahedron[h].nodeNumbering[j]].coordinate[1] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[1])
                                                    * ((globalNode[globalTetrahedron[h].nodeNumbering[(j - 1) % 4]].coordinate[0] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[0])
                                                            * (globalNode[globalTetrahedron[h].nodeNumbering[(j + 2) % 4]].coordinate[2] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[2])
                                                        - (globalNode[globalTetrahedron[h].nodeNumbering[(j + 2) % 4]].coordinate[0] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[0])
                                                            * (globalNode[globalTetrahedron[h].nodeNumbering[(j - 1) % 4]].coordinate[2] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[2]))
                                                + (globalNode[globalTetrahedron[h].nodeNumbering[j]].coordinate[2] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[2])
                                                    * ((globalNode[globalTetrahedron[h].nodeNumbering[(j + 2) % 4]].coordinate[0] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[0])
                                                            * (globalNode[globalTetrahedron[h].nodeNumbering[(j - 1) % 4]].coordinate[1] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[1])
                                                        - (globalNode[globalTetrahedron[h].nodeNumbering[(j - 1) % 4]].coordinate[0] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[0])
                                                            * (globalNode[globalTetrahedron[h].nodeNumbering[(j + 2) % 4]].coordinate[1] - globalNode[globalTetrahedron[h].nodeNumbering[(j + 1) % 4]].coordinate[1])));
	}
}

void MeshConstruction3D::printTetrahedron() const
{
	ofstream outputTetrahedron;
	string emptyspace = " ";
	outputTetrahedron.open("Tetrahedron " + to_string(globalTetrahedronNumber) + ".txt");
	outputTetrahedron << globalTetrahedronNumber << endl;
	for (int h = 0; h < globalTetrahedronNumber; h++)
		outputTetrahedron << showpoint << setprecision(10) << setw(20) << globalTetrahedron[h].volume
						<< setw(10) << globalTetrahedron[h].boundaryVertex << setw(10) << globalTetrahedron[h].boundaryType << endl;
	outputTetrahedron.close();
}

void MeshConstruction3D::printNode() const
{
	ofstream outputNode;
	string emptyspace = " ";
	outputNode.open("Node " + to_string(globalNodeNumber) + ".txt");
	outputNode << globalNodeNumber << endl;
    for (int i = 0; i < globalNodeNumber; i++)
		outputNode << showpoint << setprecision(10) << setw(6) << globalNode[i].boundary << endl;
	outputNode.close();
}

void MeshConstruction3D::printCrossSection() const
{
    ofstream outputCrossSection;
	string emptyspace = " ";
	outputCrossSection.open("CrossSection " + to_string(xNumber) + ".txt");
	outputCrossSection << crossSectionNumber << endl;
    for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
        outputCrossSection << showpoint << setprecision(10)
                        << setw(20) << crossSection[crosssection].coordinate[0]
                        << setw(20) << crossSection[crosssection].coordinate[1]
                        << setw(20) << crossSection[crosssection].coordinate[2]
                        << setw(10) << crossSection[crosssection].node << endl;
	outputCrossSection.close();
}

void MeshConstruction3D::printGmsh() const
{
	ofstream outputGmsh;
	outputGmsh.open("Gmsh_3D.msh");

	outputGmsh << "$MeshFormat" << endl;
	outputGmsh << "4.1" << setw(12) << "0" << setw(12) << "8" << endl;
	outputGmsh << "$EndMeshFormat" << endl;

	outputGmsh << "$Nodes" << endl;
	outputGmsh << "1" << setw(12) << globalNodeNumber << setw(12) << "1" << setw(12) << globalNodeNumber << endl;
	outputGmsh << "3" << setw(12) << "1" << setw(12) << "0" << setw(12) << globalNodeNumber << endl;
	for (int i = 0; i < globalNodeNumber; i++)
		outputGmsh << (i + 1) << endl;
	for (int i = 0; i < globalNodeNumber; i++)
		outputGmsh << showpoint << setprecision(10)
					<< globalNode[i].coordinate[0] << setw(20) << globalNode[i].coordinate[1] << setw(20) << globalNode[i].coordinate[2] << endl;
	outputGmsh << "$EndNodes" << endl;

	outputGmsh << "$Elements" << endl;
	outputGmsh << "1" << setw(12) << globalTetrahedronNumber << setw(12) << "1" << setw(12) << globalTetrahedronNumber << endl;
	outputGmsh << "3" << setw(12) << "1" << setw(12) << "4" << setw(12) << globalTetrahedronNumber << endl;
	for (int h = 0; h < globalTetrahedronNumber; h++)
		outputGmsh << (h + 1)
						<< setw(12) << globalTetrahedron[h].nodeNumbering[0] << setw(12) << globalTetrahedron[h].nodeNumbering[1]
                        << setw(12) << globalTetrahedron[h].nodeNumbering[2] << setw(12) << globalTetrahedron[h].nodeNumbering[3] << endl;
	outputGmsh << "$EndElements" << endl;
}

void MeshConstruction3D::generateGrid()
{
	constructNode();
	constructTetrahedron();
	determineBoundaryVertex();
	mappingCubeToSphere();
	constructCrossSection();
	if (randomization)
        randomizeGrid();
	calculateTetrahedronVolume();

	printTetrahedron();
	printNode();
    printCrossSection();
	printGmsh();
}

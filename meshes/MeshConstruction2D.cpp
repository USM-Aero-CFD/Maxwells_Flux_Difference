#include "MeshConstruction2D.h"

MeshConstruction2D::MeshConstruction2D(const int & x_Number, const int & y_Number, const bool & Randomization)
{
	xNumber = x_Number;
	yNumber = y_Number;
	randomization = Randomization;
	xInitial = - 1.0;
	xLast = 1.0;
    yInitial = - sqrt(3.0) / 2.0;
	yLast = sqrt(3.0) / 2.0;
	xDelta = static_cast <double> ((xLast - xInitial) / xNumber);
	yDelta = static_cast <double> ((yLast - yInitial) / yNumber);

    wedgeAngle = M_PI / 20.0;
	//elementNumber = 2 * xNumber * yNumber;
	//nodeNumber = (xNumber + 1) * (yNumber + 1) + (xNumber / 2);
	elementNumber = 2 * xNumber * yNumber + yNumber;
	nodeNumber = (xNumber + 1) * (yNumber + 1) + xNumber / 2 + yNumber / 2;

	element = new elementArray [elementNumber];
	node = new nodeArray [nodeNumber];
}

MeshConstruction2D::~MeshConstruction2D()
{
	cout << xNumber << "     Destruct GRID" << endl;
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

void MeshConstruction2D::constructNode()
{
		// node[i].boundary = 'A' is inlet
		// node[i].boundary = 'B' is PEC
		// node[i].boundary = 'C' is outlet
		// node[i].boundary = 'D' is side boundary
		// node[i].boundary = 'E' is the singular point
		// node[i].boundary = 'N' is not on the boundary

	int i = 0;
    int iy = 0;
    for (int ix = 0; ix <= xNumber; ix++)
    {
        node[i].coordinate = new double [2];
        node[i].coordinate[0] = xInitial + ix * xDelta;
        node[i].coordinate[1] = yInitial + iy * yDelta;
        node[i].boundary = 'D';
        i++;
    };

    for (int iy = 1; iy <= yNumber / 2; iy++)
    {
        if (iy % 2 == 0)
        {
            for (int ix = 0; ix <= xNumber / 2; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta;
                node[i].coordinate[1] = yInitial + iy * yDelta;
                if (ix == 0)
                    node[i].boundary = 'A';
				else if (ix == xNumber / 2 && iy == yNumber / 2)
					node[i].boundary = 'E';
                else
                    node[i].boundary = 'N';
                i++;
            };
            for (int ix = xNumber / 2 + 1; ix <= xNumber; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta;
                double tangentAngle = ((yInitial + iy * yDelta) - (yInitial + iy *  ((0.0 - yInitial) - (0.0 - xInitial) * tan(wedgeAngle)) / (yNumber / 2))) / (xLast - 0.0);
                node[i].coordinate[1] = (yInitial + iy * yDelta) - ((ix - xNumber / 2) * xDelta) * tangentAngle;
				if (ix == xNumber)
                    node[i].boundary = 'C';
                else if ((ix > xNumber / 2 && ix < xNumber) && iy == yNumber / 2)
                    node[i].boundary = 'B';
                else
                    node[i].boundary = 'N';
                i++;
            };
        }
        else if (iy % 2 == 1)
        {
            for (int ix = 0; ix <= xNumber / 2; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta - xDelta / 2.0;
                node[i].coordinate[1] = yInitial + iy * yDelta;
                if (ix == 0)
                    node[i].boundary = 'A';
                else
                    node[i].boundary = 'N';
                i++;
            };
            for (int ix = xNumber / 2 + 1; ix <= xNumber + 1; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta - xDelta / 2.0;
                double tangentAngle = ((yInitial + iy * yDelta) - (yInitial + iy *  ((0.0 - yInitial) - (0.0 - xInitial) * tan(wedgeAngle)) / (yNumber / 2))) / (xLast - 0.0);
                node[i].coordinate[1] = (yInitial + iy * yDelta) - ((ix - xNumber / 2) * xDelta - xDelta / 2.0) * tangentAngle;
				if (ix == xNumber + 1)
                    node[i].boundary = 'C';
                else
                    node[i].boundary = 'N';
                i++;
            };
        };
    };

    iy = yNumber / 2;
    for (int ix = xNumber / 2 + 1; ix <= xNumber; ix++)
    {
        node[i].coordinate = new double [2];
        node[i].coordinate[0] = xInitial + ix * xDelta;
        double tangentAngle = ((yLast - (yNumber - iy) * (((yLast - 0.0) - (xLast - 0.0) * tan(wedgeAngle)) / (yNumber / 2))) - (yInitial + iy * yDelta)) / (xLast - 0.0);
        node[i].coordinate[1] = (yInitial + iy * yDelta) + ((ix - xNumber / 2) * xDelta) * tangentAngle;
		if (ix < xNumber)
            node[i].boundary = 'B';
        else if (ix == xNumber)
            node[i].boundary = 'C';
        i++;
    };

    for (int iy = yNumber / 2 + 1; iy < yNumber; iy++)
    {
        if (iy % 2 == 0)
        {
            for (int ix = 0; ix <= xNumber / 2; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta;
                node[i].coordinate[1] = yInitial + iy * yDelta;
                if (ix == 0)
                    node[i].boundary = 'A';
                else
                    node[i].boundary = 'N';
                i++;
            };
            for (int ix = xNumber / 2 + 1; ix <= xNumber; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta;
                double tangentAngle = ((yLast - (yNumber - iy) * (((yLast - 0.0) - (xLast - 0.0) * tan(wedgeAngle)) / (yNumber / 2))) - (yInitial + iy * yDelta)) / (xLast - 0.0);
                node[i].coordinate[1] = (yInitial + iy * yDelta) + ((ix - xNumber / 2) * xDelta) * tangentAngle;
				if (ix == xNumber)
                    node[i].boundary = 'C';
                else
                    node[i].boundary = 'N';
                i++;
            };
        }
        else if (iy % 2 == 1)
        {
            for (int ix = 0; ix <= xNumber / 2; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta - xDelta / 2.0;
                node[i].coordinate[1] = yInitial + iy * yDelta;
                if (ix == 0)
                    node[i].boundary = 'A';
                else
                    node[i].boundary = 'N';
                i++;
            };
            for (int ix = xNumber / 2 + 1; ix <= xNumber + 1; ix++)
            {
                node[i].coordinate = new double [2];
                node[i].coordinate[0] = xInitial + ix * xDelta - xDelta / 2.0;
                double tangentAngle = ((yLast - (yNumber - iy) * (((yLast - 0.0) - (xLast - 0.0) * tan(wedgeAngle)) / (yNumber / 2))) - (yInitial + iy * yDelta)) / (xLast - 0.0);
                node[i].coordinate[1] = (yInitial + iy * yDelta) + ((ix - xNumber / 2) * xDelta - xDelta / 2.0) * tangentAngle;
				if (ix == xNumber + 1)
                    node[i].boundary = 'C';
                else
                    node[i].boundary = 'N';
                i++;
            };
        };
    };

    iy = yNumber;
    for (int ix = 0; ix <= xNumber; ix++)
    {
        node[i].coordinate = new double [2];
        node[i].coordinate[0] = xInitial + ix * xDelta;
        node[i].coordinate[1] = yInitial + iy * yDelta;
        node[i].boundary = 'D';
        i++;
    };
}

void MeshConstruction2D::constructElement()
{
		// element[e].boundaryType = 'A' is inlet
		// element[e].boundaryType = 'B' is PEC
		// element[e].boundaryType = 'C' is outlet
		// element[e].boundaryType = 'D' is side boundary
		// element[e].boundaryType = 'N' is not on the boundary

	int e = 0;
    for (int iy = 0; iy < yNumber / 2; iy++)
    {
        if (iy % 2 == 0)
        {
            for (int ix = 0; ix < xNumber; ix++)
            {
                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + ix;
                element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + ix;
                if (ix == 0)
                {
                    element[e].boundaryVertex = 1;
                    element[e].boundaryType = 'A';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;

                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + ix;
                element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (ix + 1);
                if (iy == 0)
                {
                    element[e].boundaryVertex = 2;
                    element[e].boundaryType = 'D';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;
            };

            element[e].globalNode = new int [3];
            element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + xNumber;
            element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber + 1);
            element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + xNumber;
            element[e].boundaryVertex = 2;
            element[e].boundaryType = 'C';
            e++;
        }
        else if (iy % 2 == 1)
        {
            for (int ix = 0; ix < xNumber; ix++)
            {
                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + ix;
                element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + ix;
                if (ix == 0)
                {
                    element[e].boundaryVertex = 1;
                    element[e].boundaryType = 'A';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;

                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (ix + 1);
                element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + ix;
                if (iy == (yNumber / 2 - 1) && ix >= xNumber / 2)
                {
                    element[e].boundaryVertex = 0;
                    element[e].boundaryType = 'B';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;
            };

            element[e].globalNode = new int [3];
            element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + xNumber;
            element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber + 1);
            element[e].globalNode[2] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + xNumber;
            element[e].boundaryVertex = 0;
            element[e].boundaryType = 'C';
            e++;
        };
    };


    int iy = yNumber / 2;
    for (int ix = 0; ix < xNumber / 2; ix++)
    {
        element[e].globalNode = new int [3];
        element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + ix;
        element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
        element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
        if (ix == 0)
        {
            element[e].boundaryVertex = 1;
            element[e].boundaryType = 'A';
        }
        else
        {
            element[e].boundaryVertex = - 1;
            element[e].boundaryType = 'N';
        };
        e++;

        element[e].globalNode = new int [3];
        element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + ix;
        element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (ix + 1);
        element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
        element[e].boundaryVertex = - 1;
        element[e].boundaryType = 'N';
        e++;
    };

    int ix = xNumber / 2;
	element[e].globalNode = new int [3];
	element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix - xNumber / 2);
	element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
	element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
	element[e].boundaryVertex = - 1;
	element[e].boundaryType = 'N';
	e++;

	element[e].globalNode = new int [3];
	element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix - xNumber / 2);
	element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
	element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
	element[e].boundaryVertex = 2;
	element[e].boundaryType = 'B';
	e++;

    for (int ix = xNumber / 2 + 1; ix < xNumber; ix++)
    {
		element[e].globalNode = new int [3];
		element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
		element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
		element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
		element[e].boundaryVertex = - 1;
		element[e].boundaryType = 'N';
		e++;

		element[e].globalNode = new int [3];
		element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
		element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
		element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
		element[e].boundaryVertex = 2;
		element[e].boundaryType = 'B';
		e++;
    };

    element[e].globalNode = new int [3];
    element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + xNumber;
    element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (xNumber + 1);
    element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + xNumber;
    element[e].boundaryVertex = 2;
    element[e].boundaryType = 'C';
    e++;



    for (int iy = yNumber / 2 + 1; iy < yNumber; iy++)
    {
        if (iy % 2 == 0)
        {
            for (int ix = 0; ix < xNumber; ix++)
            {
                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
                element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
                if (ix == 0)
                {
                    element[e].boundaryVertex = 1;
                    element[e].boundaryType = 'A';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;

                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
                element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
                element[e].boundaryVertex = - 1;
                element[e].boundaryType = 'N';
                e++;
            };

            element[e].globalNode = new int [3];
            element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + xNumber;
            element[e].globalNode[1] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (xNumber + 1);
            element[e].globalNode[2] = ((iy + 1) / 2 + 1) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + xNumber;
            element[e].boundaryVertex = 2;
            element[e].boundaryType = 'C';
            e++;
        }
        else if (iy % 2 == 1)
        {
            for (int ix = 0; ix < xNumber; ix++)
            {
                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + ix;
                element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + (xNumber / 2) + ix;
                if (ix == 0)
                {
                    element[e].boundaryVertex = 1;
                    element[e].boundaryType = 'A';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;

                element[e].globalNode = new int [3];
                element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
                element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + (xNumber / 2) + (ix + 1);
                element[e].globalNode[2] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + (xNumber / 2) + ix;
                if (iy == yNumber - 1)
                {
                    element[e].boundaryVertex = 0;
                    element[e].boundaryType = 'D';
                }
                else
                {
                    element[e].boundaryVertex = - 1;
                    element[e].boundaryType = 'N';
                };
                e++;
            };

            element[e].globalNode = new int [3];
            element[e].globalNode[0] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + xNumber;
            element[e].globalNode[1] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2) * (xNumber + 2) + (xNumber / 2) + (xNumber + 1);
            element[e].globalNode[2] = ((iy + 1) / 2) * (xNumber + 1) + (iy / 2 + 1) * (xNumber + 2) + (xNumber / 2) + xNumber;
            element[e].boundaryVertex = 0;
            element[e].boundaryType = 'C';
            e++;
        };
    };
}

void MeshConstruction2D::randomizeGrid()
{
	for (int i = 0; i < nodeNumber; i++)
	{
		if (node[i].boundary == 'N')
		{
			double displacingRadius = 0.3 * (static_cast<double> (rand()) / RAND_MAX) * xDelta;
			double displacingAngle = (static_cast<double> (rand()) / RAND_MAX) * 2.0 * M_PI;

			node[i].coordinate[0] = node[i].coordinate[0] + displacingRadius * cos(displacingAngle);
			node[i].coordinate[1] = node[i].coordinate[1] + displacingRadius * sin(displacingAngle);
		};
	}
}

void MeshConstruction2D::calculateCellArea()
{
	for (int e = 0; e < elementNumber; e++)
	{
		element[e].cellArea = 0.5 * abs(node[element[e].globalNode[0]].coordinate[0] * node[element[e].globalNode[1]].coordinate[1]
										+ node[element[e].globalNode[1]].coordinate[0] * node[element[e].globalNode[2]].coordinate[1]
										+ node[element[e].globalNode[2]].coordinate[0] * node[element[e].globalNode[0]].coordinate[1]
										- node[element[e].globalNode[0]].coordinate[1] * node[element[e].globalNode[1]].coordinate[0]
										- node[element[e].globalNode[1]].coordinate[1] * node[element[e].globalNode[2]].coordinate[0]
										- node[element[e].globalNode[2]].coordinate[1] * node[element[e].globalNode[0]].coordinate[0]);
	}
}

void MeshConstruction2D::printElement() const
{
	ofstream outputElement;
	string emptyspace = " ";
	outputElement.open("Element " + to_string(elementNumber) + ".txt");
	outputElement << elementNumber << endl;
	for (int e = 0; e < elementNumber; e++)
		outputElement << showpoint << setprecision(10) << setw(20) << element[e].cellArea
						<< setw(10) << element[e].boundaryVertex << setw(10) << element[e].boundaryType << endl;
	outputElement.close();
}

void MeshConstruction2D::printNode() const
{
	ofstream outputNode;
	string emptyspace = " ";
	outputNode.open("Node " + to_string(xNumber) + ".txt");
	outputNode << nodeNumber << endl;
	for (int i = 0; i < nodeNumber; i++)
		outputNode << showpoint << setprecision(10) << setw(6) << node[i].boundary << endl;
	outputNode.close();
}

void MeshConstruction2D::printGmsh() const
{
	ofstream outputGmsh;
	outputGmsh.open("Gmsh_2D.msh");

	outputGmsh << "$MeshFormat" << endl;
	outputGmsh << "4.1" << setw(12) << "0" << setw(12) << "8" << endl;
	outputGmsh << "$EndMeshFormat" << endl;

	outputGmsh << "$Nodes" << endl;
	outputGmsh << "1" << setw(12) << nodeNumber << setw(12) << "1" << setw(12) << nodeNumber << endl;
	outputGmsh << "2" << setw(12) << "1" << setw(12) << "0" << setw(12) << nodeNumber << endl;
	for (int i = 0; i < nodeNumber; i++)
		outputGmsh << (i + 1) << endl;
	for (int i = 0; i < nodeNumber; i++)
		outputGmsh << showpoint << setprecision(10)
					<< node[i].coordinate[0] << setw(20) << node[i].coordinate[1] << endl;
	outputGmsh << "$EndNodes" << endl;

	outputGmsh << "$Elements" << endl;
	outputGmsh << "1" << setw(12) << elementNumber << setw(12) << "1" << setw(12) << elementNumber << endl;
	outputGmsh << "2" << setw(12) << "1" << setw(12) << "2" << setw(12) << elementNumber << endl;
	for (int e = 0; e < elementNumber; e++)
		outputGmsh << (e + 1)
						<< setw(12) << element[e].globalNode[0] << setw(12) << element[e].globalNode[1] << setw(12) << element[e].globalNode[2] << endl;
	outputGmsh << "$EndElements" << endl;
}

void MeshConstruction2D::generateGrid()
{
	constructNode();
	constructElement();
	if (randomization)
		randomizeGrid();
    calculateCellArea();

    printElement();
	printNode();
	printGmsh();
}

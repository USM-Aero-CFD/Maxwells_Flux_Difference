# include "Computation2D.h"

Computation2D::~Computation2D()
{
    cout << "Destruct Computation2D" << endl;

    for (int e = 0; e < elementNumber; e++)
    {
        if (method == 'A')
        {
            for (int j = 0; j < 3; j++)
            {
                delete [] element[e].medianDualNormal[j]; element[e].medianDualNormal[j] = NULL;
                delete [] element[e].edgeTangent[j]; element[e].edgeTangent[j] = NULL;
            };
            delete [] element[e].medianDualNormal; element[e].medianDualNormal = NULL;
            delete [] element[e].edgeTangent; element[e].edgeTangent = NULL;
			delete [] element[e].centroid; element[e].centroid = NULL;
        }
        else if (method == 'B')
        {
            for (int j = 0; j < 3; j++)
            {
                delete [] element[e].inwardNormal[j]; element[e].inwardNormal[j] = NULL;
            };
            delete [] element[e].inwardNormal; element[e].inwardNormal = NULL;
			
        };

        if (method == 'B')
        {
			delete [] element[e].centroid; element[e].centroid = NULL;

			if (element[e].boundaryType == 'B' || element[e].boundaryType == 'C')
			{
				for (int ki = 0; ki < 3; ki++)
            	{
                	delete [] element[e].gradient[ki]; element[e].gradient[ki] = NULL;
            	};
            	delete [] element[e].gradient; element[e].gradient = NULL;
			}
        };
    };

    for (int i = 0; i < nodeNumber; i++)
    {
        delete [] node[i].fluxResidual; node[i].fluxResidual = NULL;
		delete [] node[i].dissipation; node[i].dissipation = NULL;

        for (int ki = 0; ki < 3; ki++)
        {
            delete [] node[i].conservedVariable[ki]; node[i].conservedVariable[ki] = NULL;
        };
        delete [] node[i].conservedVariable; node[i].conservedVariable = NULL;

		if (method == 'A')
        {
            node[i].neighbourNode.erase(node[i].neighbourNode.begin(), node[i].neighbourNode.end());
            node[i].neighbourElement.erase(node[i].neighbourElement.begin(), node[i].neighbourElement.end());

            for (int ki = 0; ki < 3; ki++)
            {
                delete [] node[i].firstDerivative[ki]; node[i].firstDerivative[ki] = NULL;
            };
            delete [] node[i].firstDerivative; node[i].firstDerivative = NULL;
        };

        if (method == 'B')
        {
            node[i].neighbourNode.erase(node[i].neighbourNode.begin(), node[i].neighbourNode.end());
            node[i].neighbourElement.erase(node[i].neighbourElement.begin(), node[i].neighbourElement.end());
        };
    };
}

complexNumber Computation2D::complexAddition(const complexNumber & complexNumber_1,
                                    const complexNumber & complexNumber_2) const
{
    complexNumber newComplexNumber;

    newComplexNumber.real = complexNumber_1.real + complexNumber_2.real;
    newComplexNumber.imaginary = complexNumber_1.imaginary + complexNumber_2.imaginary;

    return newComplexNumber;
}

complexNumber Computation2D::complexComplexMultiplication(const complexNumber & complexNumber_1,
                                    const complexNumber & complexNumber_2) const
{
    complexNumber newComplexNumber;

    newComplexNumber.real = complexNumber_1.real * complexNumber_2.real
                              - complexNumber_1.imaginary * complexNumber_2.imaginary;
    newComplexNumber.imaginary = complexNumber_1.real * complexNumber_2.imaginary
                              + complexNumber_1.imaginary * complexNumber_2.real;

    return newComplexNumber;
}

complexNumber Computation2D::scalarComplexMultiplication(const float & scalar,
                                    const complexNumber & complexNumber_1) const
{
    complexNumber newComplexNumber;

    newComplexNumber.real = scalar * complexNumber_1.real;
    newComplexNumber.imaginary = scalar * complexNumber_1.imaginary;

    return newComplexNumber;
}

void Computation2D::TMmodeSolution()
{
	for (int i = 0; i < nodeNumber; i++)
	{
		complexNumber Hx, Hy, Ez;

		Hx.real = 0.0;
		Hx.imaginary = 0.0;
		Hy.real = (- propagationCoefficient / (angularFrequency * permeability)) * cos(- propagationCoefficient * node[i].coordinate[0]);
		Hy.imaginary = (- propagationCoefficient / (angularFrequency * permeability)) * sin(- propagationCoefficient * node[i].coordinate[0]);
		Ez.real = cos(- propagationCoefficient * node[i].coordinate[0]);
		Ez.imaginary = sin(- propagationCoefficient * node[i].coordinate[0]);

		node[i].Hx = Hx;
		node[i].Hy = Hy;
		node[i].Ez = Ez;
	}
}

void Computation2D::TEmodeSolution()
{
	for (int i = 0; i < nodeNumber; i++)
	{
		complexNumber Ex, Ey, Hz;

		Ex.real = 0.0;
		Ex.imaginary = 0.0;
		Ey.real = (propagationCoefficient / (angularFrequency * permittivity)) * cos(- propagationCoefficient * node[i].coordinate[0]);
		Ey.imaginary = (propagationCoefficient / (angularFrequency * permittivity)) * sin(- propagationCoefficient * node[i].coordinate[0]);
		Hz.real = cos(- propagationCoefficient * node[i].coordinate[0]);
		Hz.imaginary = sin(- propagationCoefficient * node[i].coordinate[0]);

		node[i].Ex = Ex;
		node[i].Ey = Ey;
		node[i].Hz = Hz;
	}
}

void Computation2D::spatialSolution()
{
    switch (TMTEMode)
    {
        case 'A':
            TMmodeSolution();
            break;
        case 'B':
            TEmodeSolution();
            break;
    };
}

void Computation2D::timeDependentSolution(const int & I, const int & UVARIABLE_LEVEL, const double & Time)
{
	complexNumber timeHarmonic;
	timeHarmonic.real = cos(angularFrequency * Time);
	timeHarmonic.imaginary = sin(angularFrequency * Time);

    switch (TMTEMode)
    {
        case 'A':
            complexNumber Hx, Hy, Ez;
            Hx = complexComplexMultiplication(timeHarmonic, node[I].Hx);
            Hy = complexComplexMultiplication(timeHarmonic, node[I].Hy);
            Ez = complexComplexMultiplication(timeHarmonic, node[I].Ez);
            node[I].conservedVariable[0][UVARIABLE_LEVEL] = Hx.real;
            node[I].conservedVariable[1][UVARIABLE_LEVEL] = Hy.real;
            node[I].conservedVariable[2][UVARIABLE_LEVEL] = Ez.real;
            break;
        case 'B':
            complexNumber Ex, Ey, Hz;
            Ex = complexComplexMultiplication(timeHarmonic, node[I].Ex);
            Ey = complexComplexMultiplication(timeHarmonic, node[I].Ey);
            Hz = complexComplexMultiplication(timeHarmonic, node[I].Hz);
            node[I].conservedVariable[0][UVARIABLE_LEVEL] = Ex.real;
            node[I].conservedVariable[1][UVARIABLE_LEVEL] = Ey.real;
            node[I].conservedVariable[2][UVARIABLE_LEVEL] = Hz.real;
            break;
    };


}

void Computation2D::fieldInitialization()
{
	for (int i = 0; i < nodeNumber; i++)
	{
		// node[i].UVARIABLELEVEL[ki][0] - initial
		// node[i].UVARIABLELEVEL[ki][1] - (n - 1)
		// node[i].UVARIABLELEVEL[ki][2] - n
		// node[i].UVARIABLELEVEL[ki][3] - error
		// node[i].UVARIABLELEVEL[ki][4] - exact at 0.25 * timelast
		// node[i].UVARIABLELEVEL[ki][5] - exact at 0.5 * timelast
		// node[i].UVARIABLELEVEL[ki][6] - exact at 0.75 * timelast
		// node[i].UVARIABLELEVEL[ki][7] - exact at timelast

		// INITIAL
		timeDependentSolution(i, 0, 0.0);

		// first time step
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][0];

		// ANALYTICAL
		timeDependentSolution(i, 4, timeLast / 4.0);
		timeDependentSolution(i, 5, timeLast / 2.0);
		timeDependentSolution(i, 6, 3.0 * timeLast / 4.0);
		timeDependentSolution(i, 7, timeLast);
	}
}

double Computation2D::larger(const double & value1, const double & value2) const
{
	if (value1 >= value2)
		return value1;
	else
		return value2;
}

double Computation2D::smaller(const double & value1, const double & value2) const
{
	if (value1 <= value2)
		return value1;
	else
		return value2;
}

int Computation2D::KroneckerDelta(const int & value1, const int & value2) const
{
	if (value1 == value2)
		return 1;
	else
		return 0;
}

double Computation2D::dotProduct(double* const & vector1, double* const & vector2) const
{
	double product = 0;
	for (int ix = 0; ix < 2; ix++)
	{
		product += vector1[ix] * vector2[ix];
	};
	return product;
}

void Computation2D::constructNodeNeighbour()
{
    for (int e = 0; e < elementNumber; e++)
    {
        for (int j = 0; j < 3; j++)
        {
            int i = element[e].globalNode[j];
            int size = node[i].neighbourElement.size();

            if (size == 0)
            {
                node[i].neighbourElement.resize(1);
                node[i].neighbourElement[0] = e;
            }
            else
            {
                int neighbourCell [] = { e };
                int n = 0;
                while(n < size)
                {
                    if (neighbourCell[0] == node[i].neighbourElement[n])
                        break;
                    if (neighbourCell[0] < node[i].neighbourElement[n])
                    {
                        node[i].neighbourElement.insert(node[i].neighbourElement.begin() + n, neighbourCell, neighbourCell + 1);
                        break;
                    };
                    n++;
                };
                if (n == size)
                    node[i].neighbourElement.push_back(neighbourCell[0]);
            };
        }
    }

    for (int e = 0; e < elementNumber; e++)
    {
        for (int j = 0; j < 3; j++)
        {
            int i = element[e].globalNode[j];

            for (int local = 0; local < 3; local++)
            {
                int global = element[e].globalNode[local];
                int size = node[i].neighbourNode.size();

                if (size == 0)
                {
                    node[i].neighbourNode.resize(1);
                    node[i].neighbourNode[0] = global;
                }
                else
                {
                    int neighbourNode [] = { global };
                    int n = 0;
                    while(n < size)
                    {
                        if (neighbourNode[0] == node[i].neighbourNode[n])
                            break;
                        if (neighbourNode[0] < node[i].neighbourNode[n])
                        {
                            node[i].neighbourNode.insert(node[i].neighbourNode.begin() + n, neighbourNode, neighbourNode + 1);
                            break;
                        };
                        n++;
                    };
                    if (n == size)
                        node[i].neighbourNode.push_back(neighbourNode[0]);
                }
            }
        }
    }
}

void Computation2D::constructElementNeighbour()
{
	vector<int>* neighbourCell;
    neighbourCell = new vector<int> [nodeNumber];

	for (int e = 0; e < elementNumber; e++)
	{
		for (int j = 0; j < 3; j++)
		{
			int i = element[e].globalNode[j];
            int size = neighbourCell[i].size();

		    if (size == 0)
            {
                neighbourCell[i].resize(1);
                neighbourCell[i][0] = e;
            }
            else
            {
                int n = 0;
                while(n < size)
                {
					if (e == neighbourCell[i][n])
                        break;
                    if (e < neighbourCell[i][n])
                    {
                        neighbourCell[i].insert(neighbourCell[i].begin() + n, e);
                        break;
                    }
                    n++;
                };
                if (n == size)
                    neighbourCell[i].push_back(e);
            }
		}	
	};

	for (int e = 0; e < elementNumber; e++)
	{
		vector<int> elementNeighbour;

		for (int j = 0; j < 3; j++)
		{
			int i = element[e].globalNode[j];

			int sizeNode = neighbourCell[i].size();
			for (int m = 0; m < sizeNode; m++)
			{
               	int n = 0;
				int sizeElement = elementNeighbour.size();	// the sizeElement will increase by 1 for every m-iteration

                while(n < sizeElement)
                {
					if (neighbourCell[i][m] == elementNeighbour[n])
                    {
                       	elementNeighbour.insert(elementNeighbour.begin() + n, neighbourCell[i][m]);
						if (e != neighbourCell[i][m])
							element[e].neighbourElement.push_back(neighbourCell[i][m]);
                       	break;
                    }
					else if (neighbourCell[i][m] < elementNeighbour[n])
                    {
                       	elementNeighbour.insert(elementNeighbour.begin() + n, neighbourCell[i][m]);
                       	break;
                    };
					n++;
				};
				if (n == sizeElement)
                   	elementNeighbour.push_back(neighbourCell[i][m]);
            };
		};

		elementNeighbour.erase(elementNeighbour.begin(), elementNeighbour.end());	
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		neighbourCell[i].erase(neighbourCell[i].begin(), neighbourCell[i].end());
	};
	delete [] neighbourCell; neighbourCell = NULL;
}

void Computation2D::calculateEdgeTangent()
{
	for (int e = 0; e < elementNumber; e++)
	{
		element[e].edgeTangent[0][0] = node[element[e].globalNode[2]].coordinate[0] - node[element[e].globalNode[1]].coordinate[0];
		element[e].edgeTangent[0][1] = node[element[e].globalNode[2]].coordinate[1] - node[element[e].globalNode[1]].coordinate[1];
		element[e].edgeTangent[1][0] = node[element[e].globalNode[0]].coordinate[0] - node[element[e].globalNode[2]].coordinate[0];
		element[e].edgeTangent[1][1] = node[element[e].globalNode[0]].coordinate[1] - node[element[e].globalNode[2]].coordinate[1];
		element[e].edgeTangent[2][0] = node[element[e].globalNode[1]].coordinate[0] - node[element[e].globalNode[0]].coordinate[0];
		element[e].edgeTangent[2][1] = node[element[e].globalNode[1]].coordinate[1] - node[element[e].globalNode[0]].coordinate[1];
	}
}

void Computation2D::calculateMedianDualNormal()
{
	for (int e = 0; e < elementNumber; e++)
	{
		double* midpoint;
		midpoint = new double [2];
		midpoint[0] = (node[element[e].globalNode[0]].coordinate[0] + node[element[e].globalNode[1]].coordinate[0] + node[element[e].globalNode[2]].coordinate[0]) / 3.0;
		midpoint[1] = (node[element[e].globalNode[0]].coordinate[1] + node[element[e].globalNode[1]].coordinate[1] + node[element[e].globalNode[2]].coordinate[1]) / 3.0;

		element[e].medianDualNormal[0][0] = midpoint[1] - (node[element[e].globalNode[2]].coordinate[1] + node[element[e].globalNode[1]].coordinate[1]) / 2.0;
		element[e].medianDualNormal[0][1] = - (midpoint[0] - (node[element[e].globalNode[2]].coordinate[0] + node[element[e].globalNode[1]].coordinate[0]) / 2.0);
		element[e].medianDualNormal[1][0] = midpoint[1] - (node[element[e].globalNode[0]].coordinate[1] + node[element[e].globalNode[2]].coordinate[1]) / 2.0;
		element[e].medianDualNormal[1][1] = - (midpoint[0] - (node[element[e].globalNode[0]].coordinate[0] + node[element[e].globalNode[2]].coordinate[0]) / 2.0);
		element[e].medianDualNormal[2][0] = midpoint[1] - (node[element[e].globalNode[1]].coordinate[1] + node[element[e].globalNode[0]].coordinate[1]) / 2.0;
		element[e].medianDualNormal[2][1] = - (midpoint[0] - (node[element[e].globalNode[1]].coordinate[0] + node[element[e].globalNode[0]].coordinate[0]) / 2.0);
		delete [] midpoint; midpoint = NULL;
	}
}

void Computation2D::constructInwardNormal()
{
	for (int e = 0; e < elementNumber; e++)
	{
		element[e].inwardNormal[0][0] = node[element[e].globalNode[1]].coordinate[1] - node[element[e].globalNode[2]].coordinate[1];
		element[e].inwardNormal[0][1] = node[element[e].globalNode[2]].coordinate[0] - node[element[e].globalNode[1]].coordinate[0];
		element[e].inwardNormal[1][0] = node[element[e].globalNode[2]].coordinate[1] - node[element[e].globalNode[0]].coordinate[1];
		element[e].inwardNormal[1][1] = node[element[e].globalNode[0]].coordinate[0] - node[element[e].globalNode[2]].coordinate[0];
		element[e].inwardNormal[2][0] = node[element[e].globalNode[0]].coordinate[1] - node[element[e].globalNode[1]].coordinate[1];
		element[e].inwardNormal[2][1] = node[element[e].globalNode[1]].coordinate[0] - node[element[e].globalNode[0]].coordinate[0];
	}
}

void Computation2D::medianCellArea()
{
	for (int i = 0; i < nodeNumber; i++)
        node[i].nodeArea = 0.0;

    for (int e = 0; e < elementNumber; e++)
		for (int j = 0; j < 3; j++)
        {
            int i = element[e].globalNode[j];
            node[i].nodeArea = node[i].nodeArea + (1.0 / 3.0) * element[e].cellArea;
        };
}

void Computation2D::calculateCentroid()
{
	for (int e = 0; e < elementNumber; e++)
    {
        for (int coor = 0; coor < 2; coor++)
        {
            element[e].centroid[coor] = 0.0;

            for (int j = 0; j < 3; j++)
                element[e].centroid[coor] = element[e].centroid[coor] + node[element[e].globalNode[j]].coordinate[coor];
            element[e].centroid[coor] = element[e].centroid[coor] / 3.0;
        }
    }
}

double Computation2D::determinant(double** const & matrix) const
{
	return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2])
			- matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[2][0] * matrix[1][2])
			+ matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]);
}

void Computation2D::inverseMatrix(double** & matrix) const
{
	double determinant = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2])
												- matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[2][0] * matrix[1][2])
												+ matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]);

	double m00, m01, m02;
	double m10, m11, m12;
	double m20, m21, m22;

	if (determinant == 0.0)
	{
		m00 = 0.0;
		m01 = 0.0;
		m02 = 0.0;

		m10 = 0.0;
		m11 = 0.0;
		m12 = 0.0;

		m20 = 0.0;
		m21 = 0.0;
		m22 = 0.0;
	}
	else
	{
		m00 = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant;
		m01 = - (matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]) / determinant;
		m02 = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant;

		m10 = - (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) / determinant;
		m11 = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant;
		m12 = - (matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0]) / determinant;

		m20 = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant;
		m21 = - (matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]) / determinant;
		m22 = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant;
	};

	matrix[0][0] = m00;
	matrix[0][1] = m01;
	matrix[0][2] = m02;

	matrix[1][0] = m10;
	matrix[1][1] = m11;
	matrix[1][2] = m12;

	matrix[2][0] = m20;
	matrix[2][1] = m21;
	matrix[2][2] = m22;
}

double Computation2D::globalTimeStep()
{
	double* elementTimeStep;
	elementTimeStep = new double [elementNumber];

	for (int e = 0; e < elementNumber; e++)
		elementTimeStep[e] = (2.0 / 3.0) * element[e].cellArea / speed;


	double minimum = 10000.00;
	double CFL = 0.8;
	for (int e = 0; e < elementNumber; e++)
		if (elementTimeStep[e] < minimum)
			minimum = elementTimeStep[e];

	minimum = CFL * minimum;

	delete [] elementTimeStep; elementTimeStep = NULL;

	timeNumber = static_cast <int> (timeLast / minimum) + 1;

	return static_cast <double> (timeLast / timeNumber);
}

void Computation2D::GaussElimination(float** const & matrix, float** & inverse, const int & size) const
{
    // initialize inverse matrix
    for (int mi = 0; mi < size; mi++)
        for (int mj = 0; mj < size; mj++)
            if (mi == mj)
                inverse[mi][mj] = 1.0;
            else
                inverse[mi][mj] = 0.0;

    // eliminate to obtain upper triangular (U decomposition) matrix
    for (int mi = 0; mi < size; mi++)
    {
        float multiplier;

        for (int mi_ = mi + 1; mi_ < size; mi_++)
        {
            multiplier = - matrix[mi_][mi] / matrix[mi][mi];
            for (int mj = 0; mj < size; mj++)
            {
                matrix[mi_][mj] = matrix[mi_][mj] + multiplier * matrix[mi][mj];
                inverse[mi_][mj] = inverse[mi_][mj] + multiplier * inverse[mi][mj];
            }
        }
    }

    // eliminate to obtain lower triangular (L decomposition) matrix
    for (int mi = size - 1; mi >= 0; mi--)
    {
        float multiplier;

        for (int mi_ = mi - 1; mi_ >= 0; mi_--)
        {
            multiplier = - matrix[mi_][mi] / matrix[mi][mi];
            for (int mj = 0; mj < size; mj++)
            {
                matrix[mi_][mj] = matrix[mi_][mj] + multiplier * matrix[mi][mj];
                inverse[mi_][mj] = inverse[mi_][mj] + multiplier * inverse[mi][mj];
            }
        }
    }

    // diagnolize the matrix
    for (int mi = 0; mi < size; mi++)
    {
        for (int mj = 0; mj < size; mj++)
            inverse[mi][mj] = inverse[mi][mj] / matrix[mi][mi];
        matrix[mi][mi] = 1.0;
    }
}

void Computation2D::matrixVectorMultiplication(float** const & matrix, float* const & vector, float* & result, const int & size) const
{
    for (int mi = 0; mi < size; mi++)
    {
        float sum = 0.0;
        for (int mj = 0; mj < size; mj++)
            sum += (matrix[mi][mj] * vector[mj]);
        result[mi] = sum;
    }
}

double Computation2D::LagrangeInterpolation(const double & xCoordinate, const double & yCoordinate, const int & e, const int & vertex) const
{
	double value;
	double x0 = node[element[e].globalNode[0]].coordinate[0];
	double y0 = node[element[e].globalNode[0]].coordinate[1];
	double x1 = node[element[e].globalNode[1]].coordinate[0];
	double y1 = node[element[e].globalNode[1]].coordinate[1];
	double x2 = node[element[e].globalNode[2]].coordinate[0];
	double y2 = node[element[e].globalNode[2]].coordinate[1];

	switch (vertex)
	{
		case 0:
			value = (0.5 / element[e].cellArea) * ((x2 - x1) * (yCoordinate - y1) - (y2 - y1) * (xCoordinate - x1));
			break;
		case 1:
			value = (0.5 / element[e].cellArea) * ((x0 - x2) * (yCoordinate - y2) - (y0 - y2) * (xCoordinate - x2));
			break;
		case 2:
			value = (0.5 / element[e].cellArea) * ((x1 - x0) * (yCoordinate - y0) - (y1 - y0) * (xCoordinate - x0));
			break;
	};

	return value;
}

void Computation2D::constructDistributionMatrix()
{
	/* case A : TM mode */
	/* case B : TE mode */
	switch (TMTEMode)
	{
		case 'A':
			for (int e = 0; e < elementNumber; e++)
			{
				double*** inflowMatrix;
				inflowMatrix = new double** [3];
				for (int j = 0; j < 3; j++)
				{
					inflowMatrix[j] = new double* [3];
					for (int ki = 0; ki < 3; ki++)
						inflowMatrix[j][ki] = new double [3];
				};

				for (int j = 0; j < 3; j++)
				{
					// Eigenvalues & Eigenvectors
					double xNormal = element[e].inwardNormal[j][0];
					double yNormal = element[e].inwardNormal[j][1];

					// INFLOW MATRIX
					inflowMatrix[j][0][0] = 0.0;
					inflowMatrix[j][0][1] = 0.0;
					inflowMatrix[j][0][2] = 0.5 * yNormal / permeability;

					inflowMatrix[j][1][0] = 0.0;
					inflowMatrix[j][1][1] = 0.0;
					inflowMatrix[j][1][2] = - 0.5 * xNormal / permeability;

					inflowMatrix[j][2][0] = 0.5 * yNormal / permittivity;
					inflowMatrix[j][2][1] = - 0.5 * xNormal / permittivity;
					inflowMatrix[j][2][2] = 0.0;
				};


				// Distribution Matrix
				for (int j = 0; j < 3; j++)
					for (int ki = 0; ki < 3; ki++)
						for (int kj = 0; kj < 3; kj++)
						{
							element[e].distributionMatrix[j][ki][kj] = (1.0 / 3.0) * KroneckerDelta(ki, kj)
																		+ (timeDelta / (2.0 * element[e].cellArea)) * inflowMatrix[j][ki][kj];
						};


				for (int j = 0; j < 3; j++)
				{
					for (int ki = 0; ki < 3; ki++)
					{
						delete [] inflowMatrix[j][ki]; inflowMatrix[j][ki] = NULL;
					};
					delete [] inflowMatrix[j]; inflowMatrix[j] = NULL;
				};
				delete [] inflowMatrix; inflowMatrix = NULL;
			};
			break;
		case 'B':
			for (int e = 0; e < elementNumber; e++)
			{
				double*** inflowMatrix;
				inflowMatrix = new double** [3];
				for (int j = 0; j < 3; j++)
				{
					inflowMatrix[j] = new double* [3];
					for (int ki = 0; ki < 3; ki++)
						inflowMatrix[j][ki] = new double [3];
				};

				for (int j = 0; j < 3; j++)
				{
					// Eigenvalues & Eigenvectors
					double xNormal = element[e].inwardNormal[j][0];
					double yNormal = element[e].inwardNormal[j][1];

					// INFLOW MATRIX
					inflowMatrix[j][0][0] = 0.0;
					inflowMatrix[j][0][1] = 0.0;
					inflowMatrix[j][0][2] = - 0.5 * yNormal / permittivity;

					inflowMatrix[j][1][0] = 0.0;
					inflowMatrix[j][1][1] = 0.0;
					inflowMatrix[j][1][2] = 0.5 * xNormal / permittivity;

					inflowMatrix[j][2][0] = - 0.5 * yNormal / permeability;
					inflowMatrix[j][2][1] = 0.5 * xNormal / permeability;
					inflowMatrix[j][2][2] = 0.0;
				};


				// Distribution Matrix
				for (int j = 0; j < 3; j++)
					for (int ki = 0; ki < 3; ki++)
						for (int kj = 0; kj < 3; kj++)
						{
							element[e].distributionMatrix[j][ki][kj] = (1.0 / 3.0) * KroneckerDelta(ki, kj)
																		+ (timeDelta / (2.0 * element[e].cellArea)) * inflowMatrix[j][ki][kj];
						};


				for (int j = 0; j < 3; j++)
				{
					for (int ki = 0; ki < 3; ki++)
					{
						delete [] inflowMatrix[j][ki]; inflowMatrix[j][ki] = NULL;
					};
					delete [] inflowMatrix[j]; inflowMatrix[j] = NULL;
				};
				delete [] inflowMatrix; inflowMatrix = NULL;
			};
			break;
	};
}

void Computation2D::calculateFiniteVolume(const int & UVARIABLE_LEVEL)
{
	for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 3; ki++)
			node[i].fluxResidual[ki] = 0.0;

	double* displacement = new double [2];
	double* vector = new double [2];
	float* fluxResidual = new float [3];

	switch (TMTEMode)
	{
        case 'A':
            for (int e = 0; e < elementNumber; e++)
            {
                for (int j = 0; j < 3; j++)
                {
					double xMedianDualNormal;
                    double yMedianDualNormal;
					float xCoordinate, yCoordinate;
                    double xMagneticField;
                    double yMagneticField;
                    double zElectricField;

                    int i = element[e].globalNode[j];
					int i1 = element[e].globalNode[(j + 1) % 3];
					int i2 = element[e].globalNode[(j + 2) % 3];

					// local node j
					vector[0] = element[e].edgeTangent[j][0] / 2.0 - element[e].medianDualNormal[j][1] / 2.0;
					vector[1] = element[e].edgeTangent[j][1] / 2.0 + element[e].medianDualNormal[j][0] / 2.0;

					// positive side of median dual normal
					xMedianDualNormal = element[e].medianDualNormal[j][0];
                    yMedianDualNormal = element[e].medianDualNormal[j][1];
					xCoordinate = (node[i1].coordinate[0] + vector[0]);
					yCoordinate = (node[i1].coordinate[1] + vector[1]);
					for (int coor = 0; coor < 2; coor++)
						displacement[coor] = (node[i1].coordinate[coor] + vector[coor]) - node[i1].coordinate[coor];
                    xMagneticField = node[i1].conservedVariable[0][UVARIABLE_LEVEL] + dotProduct(displacement, node[i1].firstDerivative[0]);
                    yMagneticField = node[i1].conservedVariable[1][UVARIABLE_LEVEL] + dotProduct(displacement, node[i1].firstDerivative[1]);
                    zElectricField = node[i1].conservedVariable[2][UVARIABLE_LEVEL] + dotProduct(displacement, node[i1].firstDerivative[2]);
					// xMagneticField = node[i1].conservedVariable[0][UVARIABLE_LEVEL];
                    // yMagneticField = node[i1].conservedVariable[1][UVARIABLE_LEVEL];
                    // zElectricField = node[i1].conservedVariable[2][UVARIABLE_LEVEL];
					fluxResidual[0] = (- 0.0 * xMedianDualNormal + (yMedianDualNormal / permeability) * zElectricField);
                    fluxResidual[1] = (- (xMedianDualNormal / permeability) * zElectricField + 0.0 * yMedianDualNormal);
                    fluxResidual[2] = (- (xMedianDualNormal / permittivity) * yMagneticField + (yMedianDualNormal / permittivity) * xMagneticField);
					// flux out from local node (j + 1) % 3
					node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + 0.5 * fluxResidual[0];
                    node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + 0.5 * fluxResidual[1];
                    node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + 0.5 * fluxResidual[2];
					// flux into local node (j + 2) % 3
					node[i2].fluxResidual[0] = node[i2].fluxResidual[0] - 0.5 * fluxResidual[0];
                    node[i2].fluxResidual[1] = node[i2].fluxResidual[1] - 0.5 * fluxResidual[1];
                    node[i2].fluxResidual[2] = node[i2].fluxResidual[2] - 0.5 * fluxResidual[2];

					// negative side of median dual normal
					xMedianDualNormal = - element[e].medianDualNormal[j][0];
                    yMedianDualNormal = - element[e].medianDualNormal[j][1];
					for (int coor = 0; coor < 2; coor++)
						displacement[coor] = (node[i1].coordinate[coor] + vector[coor]) - node[i2].coordinate[coor];
                    xMagneticField = node[i2].conservedVariable[0][UVARIABLE_LEVEL] + dotProduct(displacement, node[i2].firstDerivative[0]);
                    yMagneticField = node[i2].conservedVariable[1][UVARIABLE_LEVEL] + dotProduct(displacement, node[i2].firstDerivative[1]);
                    zElectricField = node[i2].conservedVariable[2][UVARIABLE_LEVEL] + dotProduct(displacement, node[i2].firstDerivative[2]);
					// xMagneticField = node[i2].conservedVariable[0][UVARIABLE_LEVEL];
                    // yMagneticField = node[i2].conservedVariable[1][UVARIABLE_LEVEL];
                    // zElectricField = node[i2].conservedVariable[2][UVARIABLE_LEVEL];
					fluxResidual[0] = (- 0.0 * xMedianDualNormal + (yMedianDualNormal / permeability) * zElectricField);
                    fluxResidual[1] = (- (xMedianDualNormal / permeability) * zElectricField + 0.0 * yMedianDualNormal);
                    fluxResidual[2] = (- (xMedianDualNormal / permittivity) * yMagneticField + (yMedianDualNormal / permittivity) * xMagneticField);
					// flux into local node (j + 1) % 3
					node[i1].fluxResidual[0] = node[i1].fluxResidual[0] - 0.5 * fluxResidual[0];
                    node[i1].fluxResidual[1] = node[i1].fluxResidual[1] - 0.5 * fluxResidual[1];
                    node[i1].fluxResidual[2] = node[i1].fluxResidual[2] - 0.5 * fluxResidual[2];
					// flux out from local node (j + 2) % 3
					node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + 0.5 * fluxResidual[0];
                    node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + 0.5 * fluxResidual[1];
                    node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + 0.5 * fluxResidual[2];
                }
            };
            break;
        case 'B':
            for (int e = 0; e < elementNumber; e++)
            {
                for (int j = 0; j < 3; j++)
                {
					double xMedianDualNormal;
                    double yMedianDualNormal;
					float xCoordinate, yCoordinate;
                    double xElectricField;
                    double yElectricField;
                    double zMagneticField;

                    int i = element[e].globalNode[j];
					int i1 = element[e].globalNode[(j + 1) % 3];
					int i2 = element[e].globalNode[(j + 2) % 3];

					// local node j
					vector[0] = element[e].edgeTangent[j][0] / 2.0 - element[e].medianDualNormal[j][1] / 2.0;
					vector[1] = element[e].edgeTangent[j][1] / 2.0 + element[e].medianDualNormal[j][0] / 2.0;

					// positive side of median dual normal
					xMedianDualNormal = element[e].medianDualNormal[j][0];
                    yMedianDualNormal = element[e].medianDualNormal[j][1];
					xCoordinate = (node[i1].coordinate[0] + vector[0]);
					yCoordinate = (node[i1].coordinate[1] + vector[1]);
					for (int coor = 0; coor < 2; coor++)
						displacement[coor] = (node[i1].coordinate[coor] + vector[coor]) - node[i1].coordinate[coor];
                    xElectricField = node[i1].conservedVariable[0][UVARIABLE_LEVEL] + dotProduct(displacement, node[i1].firstDerivative[0]);
                    yElectricField = node[i1].conservedVariable[1][UVARIABLE_LEVEL] + dotProduct(displacement, node[i1].firstDerivative[1]);
                    zMagneticField = node[i1].conservedVariable[2][UVARIABLE_LEVEL] + dotProduct(displacement, node[i1].firstDerivative[2]);
					// xElectricField = node[i1].conservedVariable[0][UVARIABLE_LEVEL];
                    // yElectricField = node[i1].conservedVariable[1][UVARIABLE_LEVEL];
                    // zMagneticField = node[i1].conservedVariable[2][UVARIABLE_LEVEL];
					fluxResidual[0] = (0.0 * xMedianDualNormal - (yMedianDualNormal / permittivity) * zMagneticField);
                    fluxResidual[1] = ((xMedianDualNormal / permittivity) * zMagneticField - 0.0 * yMedianDualNormal);
                    fluxResidual[2] = (- (yMedianDualNormal / permeability) * xElectricField + (xMedianDualNormal / permeability) * yElectricField);
					// flux out from local node (j + 1) % 3
					node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + 0.5 * fluxResidual[0];
                    node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + 0.5 * fluxResidual[1];
                    node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + 0.5 * fluxResidual[2];
					// flux into local node (j + 2) % 3
					node[i2].fluxResidual[0] = node[i2].fluxResidual[0] - 0.5 * fluxResidual[0];
                    node[i2].fluxResidual[1] = node[i2].fluxResidual[1] - 0.5 * fluxResidual[1];
                    node[i2].fluxResidual[2] = node[i2].fluxResidual[2] - 0.5 * fluxResidual[2];

					// negative side of median dual normal
					xMedianDualNormal = - element[e].medianDualNormal[j][0];
                    yMedianDualNormal = - element[e].medianDualNormal[j][1];
					for (int coor = 0; coor < 2; coor++)
						displacement[coor] = (node[i1].coordinate[coor] + vector[coor]) - node[i2].coordinate[coor];
                    xElectricField = node[i2].conservedVariable[0][UVARIABLE_LEVEL] + dotProduct(displacement, node[i2].firstDerivative[0]);
                    yElectricField = node[i2].conservedVariable[1][UVARIABLE_LEVEL] + dotProduct(displacement, node[i2].firstDerivative[1]);
                    zMagneticField = node[i2].conservedVariable[2][UVARIABLE_LEVEL] + dotProduct(displacement, node[i2].firstDerivative[2]);
					// xElectricField = node[i2].conservedVariable[0][UVARIABLE_LEVEL];
                    // yElectricField = node[i2].conservedVariable[1][UVARIABLE_LEVEL];
                    // zMagneticField = node[i2].conservedVariable[2][UVARIABLE_LEVEL];
					fluxResidual[0] = (0.0 * xMedianDualNormal - (yMedianDualNormal / permittivity) * zMagneticField);
                    fluxResidual[1] = ((xMedianDualNormal / permittivity) * zMagneticField - 0.0 * yMedianDualNormal);
                    fluxResidual[2] = (- (yMedianDualNormal / permeability) * xElectricField + (xMedianDualNormal / permeability) * yElectricField);
					// flux into local node (j + 1) % 3
					node[i1].fluxResidual[0] = node[i1].fluxResidual[0] - 0.5 * fluxResidual[0];
                    node[i1].fluxResidual[1] = node[i1].fluxResidual[1] - 0.5 * fluxResidual[1];
                    node[i1].fluxResidual[2] = node[i1].fluxResidual[2] - 0.5 * fluxResidual[2];
					// flux out from local node (j + 2) % 3
					node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + 0.5 * fluxResidual[0];
                    node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + 0.5 * fluxResidual[1];
                    node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + 0.5 * fluxResidual[2];
                }
            };
            break;
	}

	delete [] displacement; displacement = NULL;
	delete [] vector; vector = NULL;
	delete [] fluxResidual; fluxResidual = NULL;
}

void Computation2D::calculateGradient(const int & UVARIABLE_LEVEL)
{
	if (method == 'A')
	{
		for (int i = 0; i < nodeNumber; i++)
    	{
        	float Sxx = 0.0;
        	float Syy = 0.0;
        	float Sxy = 0.0;
        	float Sxu = 0.0;
        	float Syu = 0.0;
        	float determinant;

        	float* xk = new float [node[i].neighbourNode.size()];
        	float* yk = new float [node[i].neighbourNode.size()];

        	for (int n = 0; n < node[i].neighbourNode.size(); n++)
        	{
        	    xk[n] = node[node[i].neighbourNode[n]].coordinate[0] - node[i].coordinate[0];
        	    yk[n] = node[node[i].neighbourNode[n]].coordinate[1] - node[i].coordinate[1];
        	};

        	// calculations for Sxx, Sxy, Syy
        	for (int n = 0; n < node[i].neighbourNode.size(); n++)
        	{
        	    if (node[i].neighbourNode[n] != i)
        	    {
        	        Sxx += (xk[n] * xk[n]);
        	        Syy += (yk[n] * yk[n]);
        	        Sxy += (xk[n] * yk[n]);
        	    };
        	};
        	determinant = Sxx * Syy - Sxy * Sxy;

        	// calculations for first derivative
        	for (int ki = 0; ki < 3; ki++)
        	{
        	    for (int n = 0; n < node[i].neighbourNode.size(); n++)
        	    {
        	        float uk = node[node[i].neighbourNode[n]].conservedVariable[ki][UVARIABLE_LEVEL] - node[i].conservedVariable[ki][UVARIABLE_LEVEL];

        	        if (node[i].neighbourNode[n] != i)
        	        {
        	            Sxu += (xk[n] * uk);
        	            Syu += (yk[n] * uk);
        	        };
        	    };
        	    // inverse estimator matrix mulply with Sxu, Syu
				// node[i].firstDerivative[ki][0] = (Syy * Sxu + (- Sxy) * Syu) / determinant;
        	    // node[i].firstDerivative[ki][1] = ((- Sxy) * Sxu + Sxx * Syu) / determinant;

				// gradient limiter
				if (ki == 1)
				{
					node[i].firstDerivative[ki][0] = (Syy * Sxu + (- Sxy) * Syu) / determinant;
				}
				else if (ki == 2)
				{
					node[i].firstDerivative[ki][0] = - node[i].firstDerivative[ki - 1][0];
				}
				else
				{
					node[i].firstDerivative[ki][0] = 0.0;
					node[i].firstDerivative[ki][1] = 0.0;
				}
        	};

        	delete [] xk; xk = NULL;
        	delete [] yk; yk = NULL;

        	// ///////////////////////////////////////////////////////////////
        	// // analytical first derivatives
        	// complexNumber Hxx, Hxy, Hyx, Hyy, Ezx, Ezy;
        	// Hxx.real = 0.0;
        	// Hxx.imaginary = 0.0;
        	// Hxy.real = 0.0;
        	// Hxy.imaginary = 0.0;

        	// Hyx.real = (- pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * sin(- propagationCoefficient * node[i].coordinate[0]);
        	// Hyx.imaginary = (pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * cos(- propagationCoefficient * node[i].coordinate[0]);
        	// Hyy.real = 0.0;
        	// Hyy.imaginary = 0.0;

        	// Ezx.real = propagationCoefficient * sin(- propagationCoefficient * node[i].coordinate[0]);
        	// Ezx.imaginary = - propagationCoefficient * cos(- propagationCoefficient * node[i].coordinate[0]);
        	// Ezy.real = 0.0;
        	// Ezy.imaginary = 0.0;

        	// complexNumber timeHarmonic;
        	// timeHarmonic.real = cos(angularFrequency * time);
        	// timeHarmonic.imaginary = sin(angularFrequency * time);

        	// node[i].firstDerivative[0][0] = complexComplexMultiplication(timeHarmonic, Hxx).real;
        	// node[i].firstDerivative[0][1] = complexComplexMultiplication(timeHarmonic, Hxy).real;
        	// node[i].firstDerivative[1][0] = complexComplexMultiplication(timeHarmonic, Hyx).real;
        	// node[i].firstDerivative[1][1] = complexComplexMultiplication(timeHarmonic, Hyy).real;
        	// node[i].firstDerivative[2][0] = complexComplexMultiplication(timeHarmonic, Ezx).real;
        	// node[i].firstDerivative[2][1] = complexComplexMultiplication(timeHarmonic, Ezy).real;
    	}
	}
	else if (method == 'B')
	{
    	for (int e = 0; e < elementNumber; e++)
    	{
			if (element[e].boundaryType == 'B' || element[e].boundaryType == 'C')
			{
        		for (int ki = 0; ki < 3; ki++)
        		{
            		element[e].gradient[ki][0] = 0.0;
            		element[e].gradient[ki][1] = 0.0;
            		for (int j = 0; j < 3; j++)
            		{
                		int i = element[e].globalNode[j];
                		element[e].gradient[ki][0] = element[e].gradient[ki][0] + 0.5 * (element[e].inwardNormal[j][0] / element[e].cellArea) * node[i].conservedVariable[ki][UVARIABLE_LEVEL];
                		element[e].gradient[ki][1] = element[e].gradient[ki][1] + 0.5 * (element[e].inwardNormal[j][1] / element[e].cellArea) * node[i].conservedVariable[ki][UVARIABLE_LEVEL];
            		}
       	 		}
			}
    	};

		// for (int i = 0; i < nodeNumber; i++)
    	// {
        // 	float Sxx = 0.0;
        // 	float Syy = 0.0;
        // 	float Sxy = 0.0;
        // 	float Sxu = 0.0;
        // 	float Syu = 0.0;
        // 	float determinant;

        // 	float* xk = new float [node[i].neighbourNode.size()];
        // 	float* yk = new float [node[i].neighbourNode.size()];

        // 	for (int n = 0; n < node[i].neighbourNode.size(); n++)
        // 	{
        // 	    xk[n] = node[node[i].neighbourNode[n]].coordinate[0] - node[i].coordinate[0];
        // 	    yk[n] = node[node[i].neighbourNode[n]].coordinate[1] - node[i].coordinate[1];
        // 	};

        // 	// calculations for Sxx, Sxy, Syy
        // 	for (int n = 0; n < node[i].neighbourNode.size(); n++)
        // 	{
        // 	    if (node[i].neighbourNode[n] != i)
        // 	    {
        // 	        Sxx += (xk[n] * xk[n]);
        // 	        Syy += (yk[n] * yk[n]);
        // 	        Sxy += (xk[n] * yk[n]);
        // 	    };
        // 	};
        // 	determinant = Sxx * Syy - Sxy * Sxy;

        // 	// calculations for first derivative
        // 	for (int ki = 0; ki < 3; ki++)
        // 	{
        // 	    for (int n = 0; n < node[i].neighbourNode.size(); n++)
        // 	    {
        // 	        float uk = node[node[i].neighbourNode[n]].conservedVariable[ki][UVARIABLE_LEVEL] - node[i].conservedVariable[ki][UVARIABLE_LEVEL];

        // 	        if (node[i].neighbourNode[n] != i)
        // 	        {
        // 	            Sxu += (xk[n] * uk);
        // 	            Syu += (yk[n] * uk);
        // 	        };
        // 	    };
        // 	    // inverse estimator matrix mulply with Sxu, Syu
        // 	    node[i].firstDerivative[ki][0] = (Syy * Sxu + (- Sxy) * Syu) / determinant;
        // 	    node[i].firstDerivative[ki][1] = ((- Sxy) * Sxu + Sxx * Syu) / determinant;
        // 	};

        // 	delete [] xk; xk = NULL;
        // 	delete [] yk; yk = NULL;

        // 	///////////////////////////////////////////////////////////////
        // 	// analytical first derivatives
        // 	complexNumber Hxx, Hxy, Hyx, Hyy, Ezx, Ezy;
        // 	Hxx.real = 0.0;
        // 	Hxx.imaginary = 0.0;
        // 	Hxy.real = 0.0;
        // 	Hxy.imaginary = 0.0;

        // 	Hyx.real = (- pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * sin(- propagationCoefficient * node[i].coordinate[0]);
        // 	Hyx.imaginary = (pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * cos(- propagationCoefficient * node[i].coordinate[0]);
        // 	Hyy.real = 0.0;
        // 	Hyy.imaginary = 0.0;

        // 	Ezx.real = propagationCoefficient * sin(- propagationCoefficient * node[i].coordinate[0]);
        // 	Ezx.imaginary = - propagationCoefficient * cos(- propagationCoefficient * node[i].coordinate[0]);
        // 	Ezy.real = 0.0;
        // 	Ezy.imaginary = 0.0;

        // 	complexNumber timeHarmonic;
        // 	timeHarmonic.real = cos(angularFrequency * time);
        // 	timeHarmonic.imaginary = sin(angularFrequency * time);

        // 	node[i].firstDerivative[0][0] = complexComplexMultiplication(timeHarmonic, Hxx).real;
        // 	node[i].firstDerivative[0][1] = complexComplexMultiplication(timeHarmonic, Hxy).real;
        // 	node[i].firstDerivative[1][0] = complexComplexMultiplication(timeHarmonic, Hyx).real;
        // 	node[i].firstDerivative[1][1] = complexComplexMultiplication(timeHarmonic, Hyy).real;
        // 	node[i].firstDerivative[2][0] = complexComplexMultiplication(timeHarmonic, Ezx).real;
        // 	node[i].firstDerivative[2][1] = complexComplexMultiplication(timeHarmonic, Ezy).real;
    	// }
	}
}

void Computation2D::TMFluxDifference(const int & E, const int & UVARIABLE_LEVEL)
{
	///////////////////////////////////////////////////////////
    // calculate flux residual
	double xMagneticField = 0.0;
    double yMagneticField = 0.0;
    double zElectricField = 0.0;

    for (int j = 0; j < 3; j++)
	{		
		int i = element[E].globalNode[j];
		
		xMagneticField = xMagneticField + node[i].conservedVariable[0][UVARIABLE_LEVEL];
		yMagneticField = yMagneticField + node[i].conservedVariable[1][UVARIABLE_LEVEL];
		zElectricField = zElectricField + node[i].conservedVariable[2][UVARIABLE_LEVEL];
    };
    xMagneticField = xMagneticField / 3.0;
    yMagneticField = yMagneticField / 3.0;
    zElectricField = zElectricField / 3.0;

	for (int j = 0; j < 3; j++)
	{
		int i = element[E].globalNode[j];

		double xMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][0];
		double yMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][1];

        node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- 0.0 * xMedianDualNormal + (zElectricField / permeability) * yMedianDualNormal);
		node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- (zElectricField / permeability) * xMedianDualNormal + 0.0 * yMedianDualNormal);
		node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- (yMagneticField / permittivity) * xMedianDualNormal + (xMagneticField / permittivity) * yMedianDualNormal);
	};

	// ///////////////////////////////////////////////////////////
	// // calculate numerical dissipation
    // double* secondDerivative = new double [3];
    // for (int j = 0; j < 3; j++)
	// {
	// 	int i = element[E].globalNode[j];
	// 	int i1 = element[E].globalNode[(j + 1) % 3];
	// 	int i2 = element[E].globalNode[(j + 2) % 3];

	// 	for (int ki = 0; ki < 3; ki++)
    //     {
    //         secondDerivative[ki] = (1.0 / 3.0) * (0.5 * pow(node[i1].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
    //                                 				+ 0.5 * (node[i1].coordinate[0] - node[i].coordinate[0]) * (node[i1].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
    //                                 				+ 0.5 * (node[i1].coordinate[1] - node[i].coordinate[1]) * (node[i1].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
    //                                 				+ 0.5 * pow(node[i1].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]
										
	// 												+ 0.5 * pow(node[i2].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
    //                                 				+ 0.5 * (node[i2].coordinate[0] - node[i].coordinate[0]) * (node[i2].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
    //                                 				+ 0.5 * (node[i2].coordinate[1] - node[i].coordinate[1]) * (node[i2].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
    //                                 				+ 0.5 * pow(node[i2].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]);
    //     }

	// 	double xMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][0];
	// 	double yMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][1];
		
	// 	node[i].dissipation[0] = node[i].dissipation[0] + (- 0.0 * xMedianDualNormal + (secondDerivative[2] / permeability) * yMedianDualNormal);
	// 	node[i].dissipation[1] = node[i].dissipation[1] + (- (secondDerivative[2] / permeability) * xMedianDualNormal + 0.0 * yMedianDualNormal);
	// 	node[i].dissipation[2] = node[i].dissipation[2] + (- (secondDerivative[1] / permittivity) * xMedianDualNormal + (secondDerivative[0] / permittivity) * yMedianDualNormal);

	// };
	// delete [] secondDerivative; secondDerivative = NULL;
}

void Computation2D::TEFluxDifference(const int & E, const int & UVARIABLE_LEVEL)
{
    ///////////////////////////////////////////////////////////
    // calculate flux residual
	double xElectricField = 0.0;
    double yElectricField = 0.0;
    double zMagneticField = 0.0;

    for (int j = 0; j < 3; j++)
	{		
		int i = element[E].globalNode[j];
		
		xElectricField = xElectricField + node[i].conservedVariable[0][UVARIABLE_LEVEL];
		yElectricField = yElectricField + node[i].conservedVariable[1][UVARIABLE_LEVEL];
		zMagneticField = zMagneticField + node[i].conservedVariable[2][UVARIABLE_LEVEL];
    };
    xElectricField = xElectricField / 3.0;
    yElectricField = yElectricField / 3.0;
    zMagneticField = zMagneticField / 3.0;

	for (int j = 0; j < 3; j++)
	{
		int i = element[E].globalNode[j];

		double xMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][0];
		double yMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][1];

		node[i].fluxResidual[0] = node[i].fluxResidual[0] + (0.0 * xMedianDualNormal - (zMagneticField / permittivity) * yMedianDualNormal);
		node[i].fluxResidual[1] = node[i].fluxResidual[1] + ((zMagneticField / permittivity) * xMedianDualNormal - 0.0 * yMedianDualNormal);
		node[i].fluxResidual[2] = node[i].fluxResidual[2] + ((yElectricField / permeability) * xMedianDualNormal - (xElectricField / permeability) * yMedianDualNormal);
	};

	// ///////////////////////////////////////////////////////////
	// // calculate numerical dissipation
    // double* secondDerivative = new double [3];
    // for (int j = 0; j < 3; j++)
	// {
	// 	int i = element[E].globalNode[j];
	// 	int i1 = element[E].globalNode[(j + 1) % 3];
	// 	int i2 = element[E].globalNode[(j + 2) % 3];

	// 	for (int ki = 0; ki < 3; ki++)
    //     {
    //         secondDerivative[ki] = (1.0 / 3.0) * (0.5 * pow(node[i1].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
    //                                 				+ 0.5 * (node[i1].coordinate[0] - node[i].coordinate[0]) * (node[i1].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
    //                                 				+ 0.5 * (node[i1].coordinate[1] - node[i].coordinate[1]) * (node[i1].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
    //                                 				+ 0.5 * pow(node[i1].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]
										
	// 												+ 0.5 * pow(node[i2].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
    //                                 				+ 0.5 * (node[i2].coordinate[0] - node[i].coordinate[0]) * (node[i2].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
    //                                 				+ 0.5 * (node[i2].coordinate[1] - node[i].coordinate[1]) * (node[i2].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
    //                                 				+ 0.5 * pow(node[i2].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]);
    //     }

	// 	double xMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][0];
	// 	double yMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][1];
		
	// 	node[i].dissipation[0] = node[i].dissipation[0] + (0.0 * xMedianDualNormal - (secondDerivative[2] / permittivity) * yMedianDualNormal);
	// 	node[i].dissipation[1] = node[i].dissipation[1] + ((secondDerivative[2] / permittivity) * xMedianDualNormal - 0.0 * yMedianDualNormal);
	// 	node[i].dissipation[2] = node[i].dissipation[2] + ((secondDerivative[1] / permeability) * xMedianDualNormal - (secondDerivative[0] / permeability) * yMedianDualNormal);

	// };
	// delete [] secondDerivative; secondDerivative = NULL;
}

void Computation2D::calculateFluxDifference(const int & UVARIABLE_LEVEL)
{
	for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 3; ki++)
		{
			node[i].fluxResidual[ki] = 0.0;
			node[i].dissipation[ki] = 0.0;
		}            
    
	switch (TMTEMode)
	{
		case 'A':
			for (int e = 0; e < elementNumber; e++)
			{
				TMFluxDifference(e, UVARIABLE_LEVEL);
			};
			break;
		case 'B':
			for (int e = 0; e < elementNumber; e++)
			{
				TEFluxDifference(e, UVARIABLE_LEVEL);
			};
			break;
	};
}

void Computation2D::calculateFiniteVolumeBoundaryFlux(const int & UVARIABLE_LEVEL)
{
	double* interpolateVariable = new double [3];
    double* displacement = new double [2];
	double* firstDerivative = new double [3];

	switch (TMTEMode)
	{
		case 'A':
		    for (int e = 0; e < elementNumber; e++)
            {
				int j = element[e].boundaryVertex;
				int i1 = element[e].globalNode[(j + 1) % 3];
				int i2 = element[e].globalNode[(j + 2) % 3];

				double* interpolateVariable;
				interpolateVariable = new double [3];

				double xEdgeNormal = (node[i2].coordinate[1] - node[i1].coordinate[1]) / 2.0;
				double yEdgeNormal = (node[i1].coordinate[0] - node[i2].coordinate[0]) / 2.0;

				/*********************************************/
				// Boundary on vertex (j + 2) % 3
				double xCoordinate = (3.0 * node[i2].coordinate[0] + node[i1].coordinate[0]) / 4.0;
				double yCoordinate = (3.0 * node[i2].coordinate[1] + node[i1].coordinate[1]) / 4.0;

                displacement[0] = xCoordinate - node[i2].coordinate[0];
                displacement[1] = yCoordinate - node[i2].coordinate[1];

				for (int ki = 0; ki < 3; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[i2].firstDerivative[ki]);

				for (int ki = 0; ki < 3; ki++)
					interpolateVariable[ki] = node[i2].conservedVariable[ki][UVARIABLE_LEVEL] + firstDerivative[ki];

				if (element[e].boundaryType == 'B')
				{
					// node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					// node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
				}
				else if (element[e].boundaryType == 'C')
				{
					node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
				};

				/*********************************************/
				// Boundary on vertex (j + 1) % 3
				xCoordinate = (node[i2].coordinate[0] + 3.0 * node[i1].coordinate[0]) / 4.0;
				yCoordinate = (node[i2].coordinate[1] + 3.0 * node[i1].coordinate[1]) / 4.0;

                displacement[0] = xCoordinate - node[i1].coordinate[0];
                displacement[1] = yCoordinate - node[i1].coordinate[1];

				for (int ki = 0; ki < 3; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[i1].firstDerivative[ki]);

				for (int ki = 0; ki < 3; ki++)
					interpolateVariable[ki] = node[i1].conservedVariable[ki][UVARIABLE_LEVEL] + firstDerivative[ki];

				if (element[e].boundaryType == 'B')
				{
					// node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					// node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
				}
				else if (element[e].boundaryType == 'C')
				{
					node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
				};

				delete [] interpolateVariable; interpolateVariable = NULL;
            };
		    break;
        case 'B':
            for (int e = 0; e < elementNumber; e++)
            {
				int j = element[e].boundaryVertex;
				int i1 = element[e].globalNode[(j + 1) % 3];
				int i2 = element[e].globalNode[(j + 2) % 3];

				double* interpolateVariable;
				interpolateVariable = new double [3];

				double xEdgeNormal = (node[i2].coordinate[1] - node[i1].coordinate[1]) / 2.0;
				double yEdgeNormal = (node[i1].coordinate[0] - node[i2].coordinate[0]) / 2.0;

				/*********************************************/
				// Boundary on vertex (j + 2) % 3
				double xCoordinate = (3.0 * node[i2].coordinate[0] + node[i1].coordinate[0]) / 4.0;
				double yCoordinate = (3.0 * node[i2].coordinate[1] + node[i1].coordinate[1]) / 4.0;

				displacement[0] = xCoordinate - node[i2].coordinate[0];
                displacement[1] = yCoordinate - node[i2].coordinate[1];

				for (int ki = 0; ki < 3; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[i2].firstDerivative[ki]);

				for (int ki = 0; ki < 3; ki++)
					interpolateVariable[ki] = node[i2].conservedVariable[ki][UVARIABLE_LEVEL] + firstDerivative[ki];

				if (element[e].boundaryType == 'B')
				{
					// node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					// node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
				}
				else if (element[e].boundaryType == 'C')
				{
					node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
				};

				/*********************************************/
				// Boundary on vertex (j + 1) % 3
				xCoordinate = (node[i2].coordinate[0] + 3.0 * node[i1].coordinate[0]) / 4.0;
				yCoordinate = (node[i2].coordinate[1] + 3.0 * node[i1].coordinate[1]) / 4.0;

                displacement[0] = xCoordinate - node[i1].coordinate[0];
                displacement[1] = yCoordinate - node[i1].coordinate[1];

				for (int ki = 0; ki < 3; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[i1].firstDerivative[ki]);

				for (int ki = 0; ki < 3; ki++)
					interpolateVariable[ki] = node[i1].conservedVariable[ki][UVARIABLE_LEVEL] + firstDerivative[ki];

				if (element[e].boundaryType == 'B')
				{
					// node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					// node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
				}
				else if (element[e].boundaryType == 'C')
				{
					node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
				};

				delete [] interpolateVariable; interpolateVariable = NULL;
            };
            break;
	};
	delete [] interpolateVariable; interpolateVariable = NULL;
    delete [] displacement; displacement = NULL;
	delete [] firstDerivative; firstDerivative = NULL;
}

void Computation2D::calculateFluxDifferenceBoundaryFlux(const int & UVARIABLE_LEVEL)
{
	double* interpolateVariable = new double [3];
	double* conservedVariableDelta = new double [3];
	
	switch (TMTEMode)
	{
		case 'A':
		    for (int e = 0; e < elementNumber; e++)
            {
				int j = element[e].boundaryVertex;
				int i1 = element[e].globalNode[(j + 1) % 3];
				int i2 = element[e].globalNode[(j + 2) % 3];

			    if (element[e].boundaryType == 'B' || element[e].boundaryType == 'C')
				{
					/*********************************************/
					// Boundary on vertex (j + 2) % 3
					for (int ki = 0; ki < 3; ki++)
						conservedVariableDelta[ki] = (node[i2].conservedVariable[ki][UVARIABLE_LEVEL] - node[i1].conservedVariable[ki][UVARIABLE_LEVEL]);

					for (int ki = 0; ki < 3; ki++)
						interpolateVariable[ki] = (3.0 / 4.0) * conservedVariableDelta[ki] + node[i1].conservedVariable[ki][UVARIABLE_LEVEL];

					double xEdgeNormal = (node[i2].coordinate[1] - node[i1].coordinate[1]) / 2.0;
					double yEdgeNormal = (node[i1].coordinate[0] - node[i2].coordinate[0]) / 2.0;

					if (element[e].boundaryType == 'B')
              		{
						// node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
						// node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
						node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
              		}
              		else if (element[e].boundaryType == 'C')
              		{
						node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
						node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
						node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
              		}

					/*********************************************/
					// Boundary on vertex (j + 1) % 3
                    for (int ki = 0; ki < 3; ki++)
						interpolateVariable[ki] = (1.0 / 4.0) * conservedVariableDelta[ki] + node[i1].conservedVariable[ki][UVARIABLE_LEVEL];

					if (element[e].boundaryType == 'B')
              		{
						// node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
						// node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
						node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
              		}
              		else if (element[e].boundaryType == 'C')
              		{
						node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
						node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
						node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];
              		}
			    }
            };
		    break;
        case 'B':
            for (int e = 0; e < elementNumber; e++)
            {
				int j = element[e].boundaryVertex;
				int i1 = element[e].globalNode[(j + 1) % 3];
				int i2 = element[e].globalNode[(j + 2) % 3];

				if (element[e].boundaryType == 'B' || element[e].boundaryType == 'C')
				{
					/*********************************************/
					// Boundary on vertex (j + 2) % 3
					for (int ki = 0; ki < 3; ki++)
						conservedVariableDelta[ki] = (node[i2].conservedVariable[ki][UVARIABLE_LEVEL] - node[i1].conservedVariable[ki][UVARIABLE_LEVEL]);

					for (int ki = 0; ki < 3; ki++)
						interpolateVariable[ki] = (3.0 / 4.0) * conservedVariableDelta[ki] + node[i1].conservedVariable[ki][UVARIABLE_LEVEL];

					double xEdgeNormal = (node[i2].coordinate[1] - node[i1].coordinate[1]) / 2.0;
					double yEdgeNormal = (node[i1].coordinate[0] - node[i2].coordinate[0]) / 2.0;

					if (element[e].boundaryType == 'B')
              		{
						//node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
						//node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
						node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
              		}
              		else if (element[e].boundaryType == 'C')
              		{
						node[i2].fluxResidual[0] = node[i2].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
						node[i2].fluxResidual[1] = node[i2].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
						node[i2].fluxResidual[2] = node[i2].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
              		}

					/*********************************************/
					// Boundary on vertex (j + 1) % 3
					for (int ki = 0; ki < 3; ki++)
						interpolateVariable[ki] = (1.0 / 4.0) * conservedVariableDelta[ki] + node[i1].conservedVariable[ki][UVARIABLE_LEVEL];

					if (element[e].boundaryType == 'B')
              		{
						//node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
						//node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
						node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
              		}
              		else if (element[e].boundaryType == 'C')
              		{
						node[i1].fluxResidual[0] = node[i1].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
						node[i1].fluxResidual[1] = node[i1].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
						node[i1].fluxResidual[2] = node[i1].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];
              		}
				}
            };
            break;
	};
	delete [] interpolateVariable; interpolateVariable = NULL;
	delete [] conservedVariableDelta; conservedVariableDelta = NULL;
}

double Computation2D::toleranceCalculation(double** const & totalResidual) const
{
	double variance = 0.0;
	for (int i = 0; i < nodeNumber; i++)
		variance += pow(totalResidual[i][3], 2);

	return sqrt(variance / nodeNumber);
}

void Computation2D::boundaryCondition(const int & UVARIABLE_LEVEL)
{
	for (int i = 0; i < nodeNumber; i++)
		if (node[i].boundary == 'A' || node[i].boundary == 'C' || node[i].boundary == 'D')
		    timeDependentSolution(i, 2, time);
}

void Computation2D::finiteVolumeNodalUpdate()
{
	/*/ (1/3) Simpson's Rule /**/
    double** k1 = new double* [nodeNumber];
    double** k2 = new double* [nodeNumber];
    double** k3 = new double* [nodeNumber];
    for (int i = 0; i < nodeNumber; i++) {
        k1[i] = new double [3];
        k2[i] = new double [3];
        k3[i] = new double [3];
    };
    
    /*/ Simpson's Rule Stage-1: (tn, yn) /**/
    calculateGradient(1);
    calculateFiniteVolume(1);
    calculateFiniteVolumeBoundaryFlux(1);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++) {
    		k1[i][ki] = node[i].fluxResidual[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * (k1[i][ki]);
    	};
	time -= (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-2: (t + tdelta / 2, y + tdelta * k1) /**/
    calculateGradient(2);
    calculateFiniteVolume(2);
    calculateFiniteVolumeBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++) {
    		k2[i][ki] = node[i].fluxResidual[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k2[i][ki]);
    	};
	time += (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-3: (t + tdelta, y + tdelta * k2) /**/
    calculateGradient(2);
    calculateFiniteVolume(2);
    calculateFiniteVolumeBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++) {
    		k3[i][ki] = node[i].fluxResidual[ki];
    	};
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++)
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k1[i][ki] + 4.0 * k2[i][ki] + k3[i][ki]) / 6.0;
    boundaryCondition(2);
    
    for (int i = 0; i < nodeNumber; i++) {
        delete [] k1[i]; k1[i] = NULL;
        delete [] k2[i]; k2[i] = NULL;
        delete [] k3[i]; k3[i] = NULL;
    };
    delete [] k1; k1 = NULL;
    delete [] k2; k2 = NULL;
    delete [] k3; k3 = NULL;    
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++)
    		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];


	// /*/ Runge-Kutta Method Stage-1 /**/
	// calculateGradient(1);
	// calculateFiniteVolume(1);
	// calculateFiniteVolumeBoundaryFlux(1);
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 3; ki++)
	// 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * node[i].fluxResidual[ki];
	// time -= timeDelta / 2.0;
	// boundaryCondition(2);

	// /*/ Runge-Kutta Method Stage-2 /**/
	// calculateGradient(2);
	// calculateFiniteVolume(2);
	// calculateFiniteVolumeBoundaryFlux(2);
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 3; ki++)
	// 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * node[i].fluxResidual[ki];
	// time += timeDelta / 2.0;
	// boundaryCondition(2);

	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 3; ki++)
	// 		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation2D::fluxDifferenceNodalUpdate()
{
   	// /*/ Runge-Kutta Stage-1 /**/
	// calculateGradient(1);
	// // calculateHessian(1);
	// calculateFluxDifference(1);
	// calculateFluxDifferenceBoundaryFlux(1);
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 3; ki++)
	// 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * (node[i].fluxResidual[ki] - node[i].dissipation[ki]);
	// time -= timeDelta / 2.0;
	// boundaryCondition(2);
	
	// /*/ Runge-Kutta Stage-2 /**/
	// calculateGradient(2);
	// // calculateHessian(2);
	// calculateFluxDifference(2);
	// calculateFluxDifferenceBoundaryFlux(2);
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 3; ki++)
	// 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (node[i].fluxResidual[ki] - node[i].dissipation[ki]);
	// time += timeDelta / 2.0;
	// boundaryCondition(2);
	
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 3; ki++)
	// 		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];


    // /*/ 4-stage Runge-Kutta /**/
    // double** k1 = new double* [nodeNumber];
    // double** k2 = new double* [nodeNumber];
    // double** k3 = new double* [nodeNumber];
    // double** k4 = new double* [nodeNumber];
    // for (int i = 0; i < nodeNumber; i++) {
    //     k1[i] = new double [3];
    //     k2[i] = new double [3];
    //     k3[i] = new double [3];
    //     k4[i] = new double [3];
    // };
    
    // /*/ Runge-Kutta Stage-1: (tn, yn) /**/
    // calculateGradient(1);
    // calculateHessian(1);
    // calculateFluxDifference(1);
    // calculateFluxDifferenceBoundaryFlux(1);
    // for (int i = 0; i < nodeNumber; i++)
    // 	for (int ki = 0; ki < 3; ki++) {
    // 		k1[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    // 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k1[i][ki] / 2.0);
    // 	};
	// time -= (timeDelta / 2.0);
    // boundaryCondition(2);
    
    // /*/ Runge-Kutta Stage-2: (t + tdelta / 2, y + tdelta * k1 / 2) /**/
    // calculateGradient(2);
    // calculateHessian(2);
    // calculateFluxDifference(2);
    // calculateFluxDifferenceBoundaryFlux(2);
    // for (int i = 0; i < nodeNumber; i++)
    // 	for (int ki = 0; ki < 3; ki++) {
    // 		k2[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    // 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k2[i][ki] / 2.0);
    // 	};
    // boundaryCondition(2);
    
    // /*/ Runge-Kutta Stage-3: (t + tdelta / 2, y + tdelta * k2 / 2) /**/
    // calculateGradient(2);
    // calculateHessian(2);
    // calculateFluxDifference(2);
    // calculateFluxDifferenceBoundaryFlux(2);
    // for (int i = 0; i < nodeNumber; i++)
    // 	for (int ki = 0; ki < 3; ki++) {
    // 		k3[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    // 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * k3[i][ki];
    // 	};
	// time += (timeDelta / 2.0);
    // boundaryCondition(2);
    
    // /*/ Runge-Kutta Stage-4: (t + tdelta, y + tdelta * k3) /**/
    // calculateGradient(2);
    // calculateHessian(2);
    // calculateFluxDifference(2);
    // calculateFluxDifferenceBoundaryFlux(2);
    // for (int i = 0; i < nodeNumber; i++)
    // 	for (int ki = 0; ki < 3; ki++) {
    // 		k4[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    // 	};
    
    // for (int i = 0; i < nodeNumber; i++)
    // 	for (int ki = 0; ki < 3; ki++)
    // 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k1[i][ki] + 2.0 * k2[i][ki] + 2.0 * k3[i][ki] + k4[i][ki]) / 6.0;
    // boundaryCondition(2);
    
    // for (int i = 0; i < nodeNumber; i++) {
    //     delete [] k1[i]; k1[i] = NULL;
    //     delete [] k2[i]; k2[i] = NULL;
    //     delete [] k3[i]; k3[i] = NULL;
    //     delete [] k4[i]; k4[i] = NULL;
    // };
    // delete [] k1; k1 = NULL;
    // delete [] k2; k2 = NULL;
    // delete [] k3; k3 = NULL;
    // delete [] k4; k4 = NULL;
    
    
    // for (int i = 0; i < nodeNumber; i++)
    // 	for (int ki = 0; ki < 3; ki++)
    // 		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];


	/*/ (1/3) Simpson's Rule /**/
    double** k1 = new double* [nodeNumber];
    double** k2 = new double* [nodeNumber];
    double** k3 = new double* [nodeNumber];
    for (int i = 0; i < nodeNumber; i++) {
        k1[i] = new double [3];
        k2[i] = new double [3];
        k3[i] = new double [3];
    };
    
    /*/ Simpson's Rule Stage-1: (tn, yn) /**/
    calculateGradient(1);
    // calculateHessian(1);
    calculateFluxDifference(1);
    calculateFluxDifferenceBoundaryFlux(1);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++) {
    		k1[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * (k1[i][ki]);
    	};
	time -= (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-2: (t + tdelta / 2, y + tdelta * k1) /**/
    calculateGradient(2);
    // calculateHessian(2);
    calculateFluxDifference(2);
    calculateFluxDifferenceBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++) {
    		k2[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k2[i][ki]);
    	};
	time += (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-3: (t + tdelta, y + tdelta * k2) /**/
    calculateGradient(2);
    // calculateHessian(2);
    calculateFluxDifference(2);
    calculateFluxDifferenceBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++) {
    		k3[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    	};
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++)
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * (k1[i][ki] + 4.0 * k2[i][ki] + k3[i][ki]) / 6.0;
    boundaryCondition(2);
    
    for (int i = 0; i < nodeNumber; i++) {
        delete [] k1[i]; k1[i] = NULL;
        delete [] k2[i]; k2[i] = NULL;
        delete [] k3[i]; k3[i] = NULL;
    };
    delete [] k1; k1 = NULL;
    delete [] k2; k2 = NULL;
    delete [] k3; k3 = NULL;    
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 3; ki++)
    		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation2D::errorsCalculation()
{
	for (int i = 0; i < nodeNumber; i++)
	{
		node[i].conservedVariable[0][3] = node[i].conservedVariable[0][2] - node[i].conservedVariable[0][7];
		node[i].conservedVariable[1][3] = node[i].conservedVariable[1][2] - node[i].conservedVariable[1][7];
		node[i].conservedVariable[2][3] = node[i].conservedVariable[2][2] - node[i].conservedVariable[2][7];
	};
}

void Computation2D::printResults() const
{
    string stringTime;
	if (abs(time - 0.5) <= timeDelta / 50.0)
        stringTime = "t05";
	else if (abs(time - 1.0) <= timeDelta / 50.0)
        stringTime = "t1";
    else if (abs(time - 1.5) <= timeDelta / 50.0)
        stringTime = "t15";
    else if (abs(time - 2.0) <= timeDelta / 50.0)
        stringTime = "t2";

	ofstream results;

	switch (method)
	{
        case 'A':
            results.open("Second_Order_Finite_Volume_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
        case 'B':
            results.open("Flux_Difference_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
        case 'C':
            results.open("RD_Galerkin_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
        case 'D':
            results.open("RD_Lax_Wendroff_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
	};
    results << "# vtk DataFile Version 2.0" << endl;
    results << "Wedge Scattering" << endl;
    results << "ASCII" << endl << endl;

    results << "DATASET UNSTRUCTURED_GRID" << endl;
    results << "POINTS " << nodeNumber << " float" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10)
                            << setw(20) << node[i].coordinate[0]
                            << setw(20) << node[i].coordinate[1]
                            << setw(20) << 0.0 << endl;
    results << endl;

    results << "CELLS " << elementNumber << setw(10) << elementNumber * (3 + 1) << endl;
    for (int e = 0; e < elementNumber; e++)
        results << setw(12) << "3"
                            << setw(10) << element[e].globalNode[0]
                            << setw(10) << element[e].globalNode[1]
                            << setw(10) << element[e].globalNode[2] << endl;
    results << endl;

    results << "CELL_TYPES " << elementNumber << endl;
    for (int e = 0; e < elementNumber; e++)
        results << setw(12) << "5" << endl;
    results << endl;

    results << "POINT_DATA " << nodeNumber << endl;
    results << "SCALARS Hx_Numerical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[0][2] << endl;
    results << endl;

    results << "SCALARS Hx_Analytical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[0][4] << endl;
    results << endl;

    results << "SCALARS Hy_Numerical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[1][2] << endl;
    results << endl;

    results << "SCALARS Hy_Analytical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[1][4] << endl;
    results << endl;

    results << "SCALARS Ez_Numerical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][2] << endl;
    results << endl;

    results << "SCALARS Ez_Analytical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nodeNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][4] << endl;
    results << endl;

	results.close();


  // For cross-sectional plots in Asymptote
	///////////////////////////////////////////

    if (method == 'C')
    {
        results.open("Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_Exact.txt");
        results << nodeNumber;
        for (int i = 0; i < nodeNumber; i++)
        {
            results << setw(20) << node[i].coordinate[0]
                    << setw(20) << node[i].coordinate[1]
                    << setw(20) << node[i].conservedVariable[2][5] << endl;
        };
        results.close();
    };


    switch (method)
	{
        case 'A':
            results.open("Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_SecondOrderFV.txt");
            break;
        case 'B':
            results.open("Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_FluxDifference.txt");
            break;
        case 'C':
            results.open("Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_RDGalerkin.txt");
            break;
        case 'D':
            results.open("Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_RDLW.txt");
            break;
	};
    results << nodeNumber << endl;
    for (int i = 0; i < nodeNumber; i++)
    {
        results << setw(20) << node[i].coordinate[0]
                << setw(20) << node[i].coordinate[1]
                << setw(20) << node[i].conservedVariable[2][2] << endl;
    };
    results.close();



    // to calculate L2-errors
    ///////////////////////////////////
	double averageArea = 0.0;
	for (int e = 0; e < elementNumber; e++)
		averageArea += element[e].cellArea;
	averageArea = averageArea / elementNumber;

    if (abs(time - 0.5) <= timeDelta / 50.0 || abs(time - 1.0) <= timeDelta / 50.0
        || abs(time - 1.5) <= timeDelta / 50.0 || abs(time - 2.0) <= timeDelta / 50.0)
    {
        ofstream L2Errors;
        switch (method)
        {
            case 'A':
                L2Errors.open("L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_SecondOrderFV.txt");
                break;
            case 'B':
                L2Errors.open("L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_FluxDifference.txt");
                break;
            case 'C':
                L2Errors.open("L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_RDGalerkin.txt");
                break;
            case 'D':
                L2Errors.open("L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_RDLW.txt");
                break;
        };

        double sum = 0.0;
        if (abs(time - 0.5) <= timeDelta / 50.0)
        {
            for (int i = 0; i < nodeNumber; i++)
                sum += pow(node[i].conservedVariable[2][2] - node[i].conservedVariable[2][4], 2);
        }
        else if (abs(time - 1.0) <= timeDelta / 50.0)
        {
            for (int i = 0; i < nodeNumber; i++)
                sum += pow(node[i].conservedVariable[2][2] - node[i].conservedVariable[2][5], 2);
        }
        else if (abs(time - 1.5) <= timeDelta / 50.0)
        {
            for (int i = 0; i < nodeNumber; i++)
                sum += pow(node[i].conservedVariable[2][2] - node[i].conservedVariable[2][6], 2);
        }
        else if (abs(time - 2.0) <= timeDelta / 50.0)
        {
            for (int i = 0; i < nodeNumber; i++)
                sum += pow(node[i].conservedVariable[2][2] - node[i].conservedVariable[2][7], 2);
        };
        sum = sqrt(sum / nodeNumber);

        L2Errors << setw(20) << "log XDELTA" << setw(20) << "log L2-Errors" << endl;
        L2Errors << setw(20) << log10(averageArea) << setw(20) << log10(sum) << endl;

        L2Errors.close();
    };
}

void Computation2D::intervalResults(bool & time_05, bool & time_1, bool & time_15, bool & time_2) const
{
	if (abs(time - 0.5) <= timeDelta / 50.0 && time_05)
	{
		printResults();
		time_05 = false;
	}
	else if (abs(time - 1.0) <= timeDelta / 50.0 && time_1)
	{
		printResults();
		time_1 = false;
	}
	else if (abs(time - 1.5) <= timeDelta / 50.0 && time_15)
	{
		printResults();
		time_15 = false;
	}
	else if (abs(time - 2.0) <= timeDelta / 50.0 && time_2)
	{
		printResults();
		time_2 = false;
	}
}

void Computation2D::initializeArray()
{
    for (int e = 0; e < elementNumber; e++)
	{
		if (method == 'A')
		{
			element[e].medianDualNormal = new double* [3];
			element[e].edgeTangent = new double* [3];
			for (int j = 0; j < 3; j++)
			{
				element[e].medianDualNormal[j] = new double [2];
				element[e].edgeTangent[j] = new double [2];
			};
			element[e].centroid = new float [2];
		}
		else if (method == 'B')
		{
			element[e].inwardNormal = new double* [3];
			for (int j = 0; j < 3; j++)
				element[e].inwardNormal[j] = new double [2];
		};

        if (method == 'B')
        {
            element[e].centroid = new float [2];
			if (element[e].boundaryType == 'B' || element[e].boundaryType == 'C')
            {
				element[e].gradient = new double* [3];
            	for (int ki = 0; ki < 3; ki++)
                	element[e].gradient[ki] = new double [2];
			}
        };
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		node[i].fluxResidual = new double [3];
		node[i].dissipation = new double [3];

        node[i].conservedVariable = new double* [3];
		for (int ki = 0; ki < 3; ki++)
		{
			node[i].conservedVariable[ki] = new double [8];
		};

		if (method == 'A')
        {
            node[i].firstDerivative = new double* [3];
            for (int ki = 0; ki < 3; ki++)
                node[i].firstDerivative[ki] = new double [2];
        };
	};
}

void Computation2D::timeIterations(const char & TM_TE_Mode, const char & Method, const double & Time_Last, const double & Time_Delta)
{
    TMTEMode = TM_TE_Mode;
    method = Method;
    initializeArray();

	if (method == 'A')
    {
		constructNodeNeighbour();
        calculateEdgeTangent();
		calculateMedianDualNormal();
		calculateCentroid();
    }
    else if (method == 'B')
    {
        constructNodeNeighbour();
        constructInwardNormal();
        calculateCentroid();
    }
	medianCellArea();


	timeLast = Time_Last;
	timeDelta = globalTimeStep();
    if (Time_Delta == 0.0)
	{
		timeDelta = globalTimeStep();
	}
	else
	{
		timeDelta = Time_Delta;
	}
	// timeNumber = static_cast<int> (timeLast / timeDelta) + 1;

	// timeNumber = static_cast <int> ((timeLast / 4.0) / timeDelta + 1) * 4;
    timeNumber = static_cast<int> (timeLast / timeDelta);
	time = 0.0;

	spatialSolution();
	fieldInitialization();

    bool time05 = true;
    bool time1 = true;
    bool time15 = true;
    bool time2 = true;
    ofstream outputTime;
	clock_t STARTTIME;
    switch (method)
	{
        case 'A':
            outputTime.open("Time Record (Second-Order FV) " + to_string(elementNumber) + ".txt");
            break;
        case 'B':
            outputTime.open("Time Record (Flux-Difference) " + to_string(elementNumber) + ".txt");
            break;
	};
	STARTTIME = clock();

	if (method == 'A')
	{
		for (int t = 1; t <= timeNumber; t++)
        {
            time = t * timeDelta;
			cout << "Time is " << time << endl;
            finiteVolumeNodalUpdate();

            intervalResults(time05, time1, time15, time2);
        };
		double sum = 0.0;
        ofstream L2Errors;
        L2Errors.open("L2-Errors-2D-RK2.txt");
        for (int i = 0; i < nodeNumber; i++) {
            sum += pow(node[i].conservedVariable[2][7] - node[i].conservedVariable[2][2], 2);
        };
        L2Errors << log10(sqrt(sum / nodeNumber)) << endl;
	}
    else if (method == 'B')
	{
        for (int t = 1; t <= timeNumber; t++)
        {
            time = t * timeDelta;
            cout << "Time is " << time << endl;
            fluxDifferenceNodalUpdate();

            intervalResults(time05, time1, time15, time2);
        };

		ofstream results;
        // print 2D results
        switch (method)
        {
            case 'A':
                results.open("Second_Order_Finite_Volume_2D.vtk");
                break;
            case 'B':
                results.open("Flux_Difference_2D.vtk");
                break;
        };
        results << "# vtk DataFile Version 2.0" << endl;
        results << "Wedge Scattering" << endl;
        results << "ASCII" << endl << endl;

        results << "DATASET UNSTRUCTURED_GRID" << endl;
        results << "POINTS " << nodeNumber << " float" << endl;
        for (int i = 0; i < nodeNumber; i++)
            results << showpoint << setprecision(10)
                                << setw(20) << node[i].coordinate[0]
                                << setw(20) << node[i].coordinate[1]
                                << setw(20) << 0.0 << endl;
        results << endl;

        results << "CELLS " << elementNumber << setw(10) << elementNumber * (3 + 1) << endl;
        for (int e = 0; e < elementNumber; e++)
            results << setw(12) << "3"
                                << setw(10) << element[e].globalNode[0]
                                << setw(10) << element[e].globalNode[1]
                                << setw(10) << element[e].globalNode[2] << endl;
        results << endl;

        results << "CELL_TYPES " << elementNumber << endl;
        for (int e = 0; e < elementNumber; e++)
            results << setw(12) << "5" << endl;
        results << endl;

        results << "POINT_DATA " << nodeNumber << endl;
        results << "SCALARS Ez_Numerical float 1" << endl;
        results << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nodeNumber; i++)
            results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][2] << endl;
        results << endl;

        results << "SCALARS Ez_Analytical_Q1 float 1" << endl;
        results << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nodeNumber; i++)
            results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][4] << endl;
        results << endl;

        results << "SCALARS Ez_Analytical_Q2 float 1" << endl;
        results << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nodeNumber; i++)
            results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][5] << endl;
        results << endl;

        results << "SCALARS Ez_Analytical_Q3 float 1" << endl;
        results << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nodeNumber; i++)
            results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][6] << endl;
        results << endl;

        results << "SCALARS Ez_Analytical_Q4 float 1" << endl;
        results << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nodeNumber; i++)
            results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][7] << endl;
        results << endl;
        results.close();

        double sum = 0.0;
        ofstream L2Errors;
        L2Errors.open("L2-Errors-2D-RK2.txt");
        for (int i = 0; i < nodeNumber; i++) {
            sum += pow(node[i].conservedVariable[2][7] - node[i].conservedVariable[2][2], 2);
        };
        L2Errors << log10(sqrt(sum / nodeNumber)) << endl;
	};

	outputTime << "The global time step is " << timeDelta << endl;
	outputTime << "The execution time is " << static_cast <double> ((clock() - STARTTIME) / static_cast <double> (CLOCKS_PER_SEC)) << ". " << endl;
	outputTime.close();
	errorsCalculation();
}
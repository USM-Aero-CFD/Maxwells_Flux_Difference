# include "Computation3D.h"

Computation3D::~Computation3D()
{
	cout << "Destruct Computation3D" << endl;

	for (int h = 0; h < tetrahedronNumber; h++)
	{
		for (int mediandual = 0; mediandual < 6; mediandual++)
		{
			delete [] tetrahedron[h].medianDualNormal[mediandual]; tetrahedron[h].medianDualNormal[mediandual] = NULL;
		};
		delete [] tetrahedron[h].medianDualNormal; tetrahedron[h].medianDualNormal = NULL;

		if (method == 'A')
		{
			if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
			{
				for (int j = 0; j < 4; j++)
				{
					delete [] tetrahedron[h].inwardNormal[j]; tetrahedron[h].inwardNormal[j] = NULL;
				};
				delete [] tetrahedron[h].inwardNormal; tetrahedron[h].inwardNormal = NULL;
			};
			delete [] tetrahedron[h].centroid; tetrahedron[h].centroid = NULL;
		}
		else if (method == 'B')
		{
			for (int j = 0; j < 4; j++)
			{
				delete [] tetrahedron[h].inwardNormal[j]; tetrahedron[h].inwardNormal[j] = NULL;
			};
			delete [] tetrahedron[h].inwardNormal; tetrahedron[h].inwardNormal = NULL;
		};

		if (method == 'B')
		{
			delete [] tetrahedron[h].centroid; tetrahedron[h].centroid = NULL;

			if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
			{
				for (int ki = 0; ki < 6; ki++)
				{
					delete [] tetrahedron[h].gradient[ki]; tetrahedron[h].gradient[ki] = NULL;
				};
				delete [] tetrahedron[h].gradient; tetrahedron[h].gradient = NULL;
			}
		};

		delete [] tetrahedron[h].subVolume; tetrahedron[h].subVolume = NULL;
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		delete [] node[i].fluxResidual; node[i].fluxResidual = NULL;
		delete [] node[i].dissipation; node[i].dissipation = NULL;

		for (int ki = 0; ki < 6; ki++)
		{
			delete [] node[i].conservedVariable[ki]; node[i].conservedVariable[ki] = NULL;
		};
		delete [] node[i].conservedVariable; node[i].conservedVariable = NULL;

		if (method == 'A')
		{
			node[i].neighbourNode.erase(node[i].neighbourNode.begin(), node[i].neighbourNode.end());
			node[i].neighbourTetrahedron.erase(node[i].neighbourTetrahedron.begin(), node[i].neighbourTetrahedron.end());

			for (int ki = 0; ki < 6; ki++)
			{
				delete [] node[i].firstDerivative[ki]; node[i].firstDerivative[ki] = NULL;
			};
			delete [] node[i].firstDerivative; node[i].firstDerivative = NULL;
		};

		if (method == 'B')
		{
			node[i].neighbourNode.erase(node[i].neighbourNode.begin(), node[i].neighbourNode.end());
			node[i].neighbourTetrahedron.erase(node[i].neighbourTetrahedron.begin(), node[i].neighbourTetrahedron.end());
		};
	};
}

complexNumber Computation3D::complexAddition(const complexNumber & complexNumber_1, const complexNumber & complexNumber_2) const
{
    complexNumber newComplexNumber;

    newComplexNumber.real = complexNumber_1.real + complexNumber_2.real;
    newComplexNumber.imaginary = complexNumber_1.imaginary + complexNumber_2.imaginary;

    return newComplexNumber;
}

complexNumber Computation3D::complexMultiplication(const complexNumber & complexNumber_1,
                                    const complexNumber & complexNumber_2) const
{
    complexNumber newComplexNumber;

    newComplexNumber.real = complexNumber_1.real * complexNumber_2.real
                              - complexNumber_1.imaginary * complexNumber_2.imaginary;
    newComplexNumber.imaginary = complexNumber_1.real * complexNumber_2.imaginary
                              + complexNumber_1.imaginary * complexNumber_2.real;

    return newComplexNumber;
}

complexNumber Computation3D::scalarMultiplication(const float & scalar,
                                    const complexNumber & complexNumber_1) const
{
    complexNumber newComplexNumber;

    newComplexNumber.real = scalar * complexNumber_1.real;
    newComplexNumber.imaginary = scalar * complexNumber_1.imaginary;

    return newComplexNumber;
}

complexNumber Computation3D::complexConjugate(const complexNumber & complexNumber_1) const
{
    complexNumber newComplexNumber;

    float denominator = (complexNumber_1.real * complexNumber_1.real + complexNumber_1.imaginary * complexNumber_1.imaginary);

    newComplexNumber.real = complexNumber_1.real / denominator;
    newComplexNumber.imaginary = - complexNumber_1.imaginary / denominator;

    return newComplexNumber;
}

float Computation3D::radius(const float & xCoordinate, const float & yCoordinate, const float & zCoordinate) const
{
    return sqrt(pow(xCoordinate, 2) + pow(yCoordinate, 2) + pow(zCoordinate, 2));
}

float Computation3D::theta(const float & xCoordinate, const float & yCoordinate, const float & zCoordinate) const
{
    double Theta;
    if (zCoordinate >= 0)
       Theta = atan2(sqrt(pow(xCoordinate, 2) + pow(yCoordinate, 2)), abs(zCoordinate));
    else if (zCoordinate < 0)
       Theta = M_PI - atan2(sqrt(pow(xCoordinate, 2) + pow(yCoordinate, 2)), abs(zCoordinate));

    return Theta;
}

float Computation3D::phi(const float & xCoordinate, const float & yCoordinate, const float & zCoordinate) const
{
    double Phi;
    if (xCoordinate >= 0 && yCoordinate >= 0)
       Phi = atan2(abs(yCoordinate), abs(xCoordinate));
    else if (xCoordinate < 0 && yCoordinate >= 0)
       Phi = M_PI - atan2(abs(yCoordinate), abs(xCoordinate));
    else if (xCoordinate < 0 && yCoordinate < 0)
       Phi = M_PI + atan2(abs(yCoordinate), abs(xCoordinate));
    else if (xCoordinate >= 0 && yCoordinate < 0)
       Phi = 2.0 * M_PI - atan2(abs(yCoordinate), abs(xCoordinate));

    return Phi;
}

void Computation3D::TMmodeSolution()
{
    for (int i = 0; i < nodeNumber; i++)
    {
        double xCoordinate = node[i].coordinate[0];
        double yCoordinate = node[i].coordinate[1];
        double zCoordinate = node[i].coordinate[2];

        complexNumber Ex;
        complexNumber Ey;
        complexNumber Ez;
        complexNumber Hx;
        complexNumber Hy;
        complexNumber Hz;

        Ex.real = 0.0;
        Ex.imaginary = 0.0;
        Ey.real = 0.0;
        Ey.imaginary = 0.0;
        Ez.real = 0.0;
        Ez.imaginary = 0.0;
        Hx.real = 0.0;
        Hx.imaginary = 0.0;
        Hy.real = 0.0;
        Hy.imaginary = 0.0;
        Hz.real = 0.0;
        Hz.imaginary = 0.0;

        Ez.real = cos(- propagationCoefficient * xCoordinate);
        Ez.imaginary = sin(- propagationCoefficient * xCoordinate);
        Hy.real = (- propagationCoefficient / (angularFrequency * permeability)) * cos(- propagationCoefficient * xCoordinate);
        Hy.imaginary = (- propagationCoefficient / (angularFrequency * permeability)) * sin(- propagationCoefficient * xCoordinate);

        node[i].Ex = Ex;
        node[i].Ey = Ey;
        node[i].Ez = Ez;
        node[i].Hx = Hx;
        node[i].Hy = Hy;
        node[i].Hz = Hz;
    }
}

void Computation3D::TEmodeSolution()
{
    for (int i = 0; i < nodeNumber; i++)
    {
        double xCoordinate = node[i].coordinate[0];
        double yCoordinate = node[i].coordinate[1];
        double zCoordinate = node[i].coordinate[2];

        complexNumber Ex;
        complexNumber Ey;
        complexNumber Ez;
        complexNumber Hx;
        complexNumber Hy;
        complexNumber Hz;

        Ex.real = 0.0;
        Ex.imaginary = 0.0;
        Ey.real = 0.0;
        Ey.imaginary = 0.0;
        Ez.real = 0.0;
        Ez.imaginary = 0.0;
        Hx.real = 0.0;
        Hx.imaginary = 0.0;
        Hy.real = 0.0;
        Hy.imaginary = 0.0;
        Hz.real = 0.0;
        Hz.imaginary = 0.0;

        Hz.real = cos(- propagationCoefficient * xCoordinate);
        Hz.imaginary = sin(- propagationCoefficient * xCoordinate);
        Ey.real = (propagationCoefficient / (angularFrequency * permittivity)) * cos(- propagationCoefficient * xCoordinate);
        Ey.imaginary = (propagationCoefficient / (angularFrequency * permittivity)) * sin(- propagationCoefficient * xCoordinate);

        node[i].Ex = Ex;
        node[i].Ey = Ey;
        node[i].Ez = Ez;
        node[i].Hx = Hx;
        node[i].Hy = Hy;
        node[i].Hz = Hz;
    }
}

void Computation3D::spatialSolution()
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

void Computation3D::timeDependentSolution(const int & I, const int & UVARIABLE_LEVEL, const double & Time) const
{
    complexNumber timeHarmonic;
    complexNumber Ex, Ey, Ez;
    complexNumber Hx, Hy, Hz;

    timeHarmonic.real = cos(angularFrequency * Time);
    timeHarmonic.imaginary = sin(angularFrequency * Time);

    Ex = complexMultiplication(timeHarmonic, node[I].Ex);
    Ey = complexMultiplication(timeHarmonic, node[I].Ey);
    Ez = complexMultiplication(timeHarmonic, node[I].Ez);
    Hx = complexMultiplication(timeHarmonic, node[I].Hx);
    Hy = complexMultiplication(timeHarmonic, node[I].Hy);
    Hz = complexMultiplication(timeHarmonic, node[I].Hz);

    node[I].conservedVariable[0][UVARIABLE_LEVEL] = Ex.real;
    node[I].conservedVariable[1][UVARIABLE_LEVEL] = Ey.real;
    node[I].conservedVariable[2][UVARIABLE_LEVEL] = Ez.real;
    node[I].conservedVariable[3][UVARIABLE_LEVEL] = Hx.real;
    node[I].conservedVariable[4][UVARIABLE_LEVEL] = Hy.real;
    node[I].conservedVariable[5][UVARIABLE_LEVEL] = Hz.real;
}

void Computation3D::fieldInitialization() const
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
        for (int ki = 0; ki < 6; ki++)
            node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][0];

        // ANALYTICAL
        timeDependentSolution(i, 4, timeLast / 4.0);
        timeDependentSolution(i, 5, timeLast / 2.0);
        timeDependentSolution(i, 6, 3.0 * timeLast / 4.0);
        timeDependentSolution(i, 7, timeLast);
    };

}

double Computation3D::larger(const double & value1, const double & value2) const
{
	if (value1 >= value2)
		return value1;
	else
		return value2;
}

double Computation3D::smaller(const double & value1, const double & value2) const
{
	if (value1 <= value2)
		return value1;
	else
		return value2;
}

int Computation3D::KroneckerDelta(const int & value1, const int & value2) const
{
	if (value1 == value2)
		return 1;
	else
		return 0;
}

double Computation3D::dotProduct(double* const & vector1, double* const & vector2) const
{
	double product = 0;
	for (int ix = 0; ix < 3; ix++)
	{
		product += vector1[ix] * vector2[ix];
	};
	return product;
}

double Computation3D::crossProduct(double* const & vector1, double* const & vector2, const int & coordinate) const
{
	double product;
	switch (coordinate)
	{
		case 0:
			product = (vector1[1] * vector2[2] - vector2[1] * vector1[2]);
			break;
		case 1:
			product = (vector1[2] * vector2[0] - vector2[2] * vector1[0]);
			break;
		case 2:
			product = (vector1[0] * vector2[1] - vector2[0] * vector1[1]);
			break;
	};
	return product;
}

void Computation3D::constructNodeNeighbour()
{
    for (int h = 0; h < tetrahedronNumber; h++)
    {
        for (int j = 0; j < 4; j++)
        {
            int i = tetrahedron[h].globalNode[j];
            int size = node[i].neighbourTetrahedron.size();

            if (size == 0)
            {
                node[i].neighbourTetrahedron.resize(1);
                node[i].neighbourTetrahedron[0] = h;
            }
            else
            {
                int neighbourCell [] = { h };
                int n = 0;
                while(n < size)
                {
                    if (neighbourCell[0] == node[i].neighbourTetrahedron[n])
                        break;
                    if (neighbourCell[0] < node[i].neighbourTetrahedron[n])
                    {
                        node[i].neighbourTetrahedron.insert(node[i].neighbourTetrahedron.begin() + n, neighbourCell, neighbourCell + 1);
                        break;
                    };
                    n++;
                };
                if (n == size)
                    node[i].neighbourTetrahedron.push_back(neighbourCell[0]);
            };
        }
    }

    for (int h = 0; h < tetrahedronNumber; h++)
    {
        for (int j = 0; j < 4; j++)
        {
            int i = tetrahedron[h].globalNode[j];

            for (int local = 0; local < 4; local++)
            {
                int global = tetrahedron[h].globalNode[local];
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

void Computation3D::constructTetrahedronNeighbour()
{
	vector<int>* neighbourCell;
    neighbourCell = new vector<int> [nodeNumber];

	for (int h = 0; h < tetrahedronNumber; h++)
	{
		for (int j = 0; j < 4; j++)
		{
			int i = tetrahedron[h].globalNode[j];
            int size = neighbourCell[i].size();

		    if (size == 0)
            {
                neighbourCell[i].resize(1);
                neighbourCell[i][0] = h;
            }
            else
            {
                int n = 0;
                while(n < size)
                {
					if (h == neighbourCell[i][n])
                        break;
                    if (h < neighbourCell[i][n])
                    {
                        neighbourCell[i].insert(neighbourCell[i].begin() + n, h);
                        break;
                    }
                    n++;
                };
                if (n == size)
                    neighbourCell[i].push_back(h);
            }
		}	
	};

	for (int h = 0; h < tetrahedronNumber; h++)
	{
		vector<int> tetrahedronNeighbour;

		for (int j = 0; j < 4; j++)
		{
			int i = tetrahedron[h].globalNode[j];

			int sizeNode = neighbourCell[i].size();
			for (int m = 0; m < sizeNode; m++)
			{
               	int n = 0;
				int sizeTetrahedron = tetrahedronNeighbour.size();	// the sizeTetrahedron will increase by 1 for every m-iteration

                while(n < sizeTetrahedron)
                {
					if (neighbourCell[i][m] == tetrahedronNeighbour[n])
                    {
                       	tetrahedronNeighbour.insert(tetrahedronNeighbour.begin() + n, neighbourCell[i][m]);
						if (h != neighbourCell[i][m])
							tetrahedron[h].neighbourTetrahedron.push_back(neighbourCell[i][m]);
                       	break;
                    }
					else if (neighbourCell[i][m] < tetrahedronNeighbour[n])
                    {
                       	tetrahedronNeighbour.insert(tetrahedronNeighbour.begin() + n, neighbourCell[i][m]);
                       	break;
                    };
					n++;
				};
				if (n == sizeTetrahedron)
                   	tetrahedronNeighbour.push_back(neighbourCell[i][m]);
            };
		};

		tetrahedronNeighbour.erase(tetrahedronNeighbour.begin(), tetrahedronNeighbour.end());	
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		neighbourCell[i].erase(neighbourCell[i].begin(), neighbourCell[i].end());
	};
	delete [] neighbourCell; neighbourCell = NULL;
}

// for FINITE-volume, but will be created in FLUX-DIFFERENCE AND RESIDUAL-DISTRIBUTION to calculate the CELL-volume.
void Computation3D::constructMedianDualNormal()
{
	double* positionVector1;
	double* positionVector2;
	positionVector1 = new double [3];
	positionVector2 = new double [3];

	for (int h = 0; h < tetrahedronNumber; h++)
	{
		for (int mediandual = 0; mediandual < 6; mediandual++)
			for (int coordinate = 0; coordinate < 3; coordinate++)
				tetrahedron[h].medianDualNormal[mediandual][coordinate] = 0.0;

		for (int mediandual = 0; mediandual < 6; mediandual++)
			for (int subplane = 0; subplane < 2; subplane++)
			{
				switch (mediandual)
				{
					case 0:
						switch (subplane)
						{
							case 0:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
								};
								break;
							case 1:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0;
								};
								break;
						};
						break;
					case 1:
						switch (subplane)
						{
							case 0:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
								};
								break;
							case 1:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
								};
								break;
						};
						break;
					case 2:
						switch (subplane)
						{
							case 0:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
								};
								break;
							case 1:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
								};
								break;
						};
						break;
					case 3:
						switch (subplane)
						{
							case 0:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
								};
								break;
							case 1:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
								};
								break;
						};
						break;
					case 4:
						switch (subplane)
						{
							case 0:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
								};
								break;
							case 1:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
								};
								break;
						};
						break;
					case 5:
						switch (subplane)
						{
							case 0:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0;
								};
								break;
							case 1:
								for (int coordinate = 0; coordinate < 3; coordinate++)
								{
									positionVector1[coordinate] = (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
									positionVector2[coordinate] = (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
																	- (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
								};
								break;
						};
						break;
				};
				for (int coordinate = 0; coordinate < 3; coordinate++)
					tetrahedron[h].medianDualNormal[mediandual][coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate] + 0.5 * crossProduct(positionVector1, positionVector2, coordinate);
			};
	};
	delete [] positionVector1; positionVector1 = NULL;
	delete [] positionVector2; positionVector2 = NULL;
}

// for FLUX-DIFFERENCE AND RESIDUAL-DISTRIBUTION, but will be needed in FINITE-volume just on the boundary elements.
void Computation3D::calculateInwardNormal()
{
    double* positionVector1;
	double* positionVector2;
	positionVector1 = new double [3];
	positionVector2 = new double [3];

    // Configuration B
    for (int h = 0; h < tetrahedronNumber; h++)
    {
        for (int j = 0; j < 4; j++)
        {
            switch (j)
            {
                case 0:
                    for (int coordinate = 0; coordinate < 3; coordinate++)
                    {
                        positionVector1[coordinate] = node[tetrahedron[h].globalNode[3]].coordinate[coordinate] - node[tetrahedron[h].globalNode[1]].coordinate[coordinate];
                        positionVector2[coordinate] = node[tetrahedron[h].globalNode[2]].coordinate[coordinate] - node[tetrahedron[h].globalNode[1]].coordinate[coordinate];
                    };
                    break;
                case 1:
                    for (int coordinate = 0; coordinate < 3; coordinate++)
                    {
                        positionVector1[coordinate] = node[tetrahedron[h].globalNode[2]].coordinate[coordinate] - node[tetrahedron[h].globalNode[0]].coordinate[coordinate];
                        positionVector2[coordinate] = node[tetrahedron[h].globalNode[3]].coordinate[coordinate] - node[tetrahedron[h].globalNode[0]].coordinate[coordinate];
                    };
                    break;
                case 2:
                    for (int coordinate = 0; coordinate < 3; coordinate++)
                    {
                        positionVector1[coordinate] = node[tetrahedron[h].globalNode[3]].coordinate[coordinate] - node[tetrahedron[h].globalNode[0]].coordinate[coordinate];
                        positionVector2[coordinate] = node[tetrahedron[h].globalNode[1]].coordinate[coordinate] - node[tetrahedron[h].globalNode[0]].coordinate[coordinate];
                    };
                    break;
                case 3:
                    for (int coordinate = 0; coordinate < 3; coordinate++)
                    {
                        positionVector1[coordinate] = node[tetrahedron[h].globalNode[1]].coordinate[coordinate] - node[tetrahedron[h].globalNode[0]].coordinate[coordinate];
                        positionVector2[coordinate] = node[tetrahedron[h].globalNode[2]].coordinate[coordinate] - node[tetrahedron[h].globalNode[0]].coordinate[coordinate];
                    };
                    break;
            };
            for (int coordinate = 0; coordinate < 3; coordinate++)
                tetrahedron[h].inwardNormal[j][coordinate] = 0.5 * crossProduct(positionVector1, positionVector2, coordinate);
        }
    };

    delete [] positionVector1; positionVector1 = NULL;
	delete [] positionVector2; positionVector2 = NULL;
}

void Computation3D::calculateInwardNormalFromMedianDualNormal()
{
    for (int h = 0; h < tetrahedronNumber; h++)
    {
        if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
        {
            tetrahedron[h].inwardNormal = new double* [4];
            for (int j = 0; j < 4; j++)
                tetrahedron[h].inwardNormal[j] = new double [3];

            for (int j = 0; j < 4; j++)
            {
                switch (j)
                {
                    case 0:
                        for (int coordinate = 0; coordinate < 3; coordinate++)
                        {
                            double sum = (+ tetrahedron[h].medianDualNormal[0][coordinate]
                                            + tetrahedron[h].medianDualNormal[1][coordinate]
                                            + tetrahedron[h].medianDualNormal[2][coordinate]);

                            tetrahedron[h].inwardNormal[j][coordinate] = 3.0 * sum;
                        };
                        break;
                    case 1:
                        for (int coordinate = 0; coordinate < 3; coordinate++)
                        {
                            double sum = (- tetrahedron[h].medianDualNormal[0][coordinate]
                                            + tetrahedron[h].medianDualNormal[3][coordinate]
                                            - tetrahedron[h].medianDualNormal[5][coordinate]);

                            tetrahedron[h].inwardNormal[j][coordinate] = 3.0 * sum;
                        };
                        break;
                    case 2:
                        for (int coordinate = 0; coordinate < 3; coordinate++)
                        {
                            double sum = (- tetrahedron[h].medianDualNormal[1][coordinate]
                                            - tetrahedron[h].medianDualNormal[3][coordinate]
                                            + tetrahedron[h].medianDualNormal[4][coordinate]);

                            tetrahedron[h].inwardNormal[j][coordinate] = 3.0 * sum;
                        };
                        break;
                    case 3:
                        for (int coordinate = 0; coordinate < 3; coordinate++)
                        {
                            double sum = (- tetrahedron[h].medianDualNormal[2][coordinate]
                                            - tetrahedron[h].medianDualNormal[4][coordinate]
                                            + tetrahedron[h].medianDualNormal[5][coordinate]);

                            tetrahedron[h].inwardNormal[j][coordinate] = 3.0 * sum;
                        };
                        break;
                }
            }
        };
    }
}

void Computation3D::calculateSubVolume()
{
	double* baseArea;
	double* height;
	double sum;
	baseArea = new double [3];
	height = new double [3];

	for (int h = 0; h < tetrahedronNumber; h++)
	{
		for (int j = 0; j < 4; j++)
		{
			int mediandual;
			tetrahedron[h].subVolume[j] = 0.0;

			switch (j)
			{
				case 0:
					mediandual = 0;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[0]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 1;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[0]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 2;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[0]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);
					break;
				case 1:
					mediandual = 0;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = - tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[1]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 3;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[1]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 5;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = - tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[1]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[3]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);
					break;
				case 2:
					mediandual = 1;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = - tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[2]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 3;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = - tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[2]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 4;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[2]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);
					break;
				case 3:
					mediandual = 2;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = - tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[3]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 4;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = - tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[3]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);

					mediandual = 5;
					for (int coordinate = 0; coordinate < 3; coordinate++)
					{
						baseArea[coordinate] = tetrahedron[h].medianDualNormal[mediandual][coordinate];
						height[coordinate] = node[tetrahedron[h].globalNode[3]].coordinate[coordinate] - (node[tetrahedron[h].globalNode[3]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0;
					};
					tetrahedron[h].subVolume[j] = tetrahedron[h].subVolume[j] + (1.0 / 3.0) * dotProduct(height, baseArea);
					break;
			};
		}
	};

	delete [] baseArea; baseArea = NULL;
	delete [] height; height = NULL;
}

void Computation3D::medianCellVolume()
{
	for (int i = 0; i < nodeNumber; i++)
        node[i].volume = 0.0;

	for (int h = 0; h < tetrahedronNumber; h++)
    {
        for (int j = 0; j < 4; j++)
        {
            int i = tetrahedron[h].globalNode[j];
            node[i].volume = node[i].volume + tetrahedron[h].subVolume[j];
        }
    };
}

void Computation3D::calculateCentroid()
{
	for (int h = 0; h < tetrahedronNumber; h++)
    {
		for (int coor = 0; coor < 3; coor++)
		{
			tetrahedron[h].centroid[coor] = 0.0;
			for (int j = 0; j < 4; j++)
				tetrahedron[h].centroid[coor] = tetrahedron[h].centroid[coor] + node[tetrahedron[h].globalNode[j]].coordinate[coor];
			tetrahedron[h].centroid[coor] = tetrahedron[h].centroid[coor] / 4.0;
		}
    }
}

// for ALL
void Computation3D::findCrossSection()
{
    double* vector1;
    double* vector2;
    double* vector3;
    double* vector;

    vector1 = new double [3];
    vector2 = new double [3];
    vector3 = new double [3];
    vector = new double [3];

    for (int ix = 0; ix < crossSectionNumber; ix++)
    {
        for (int h = 0; h < tetrahedronNumber; h++)
        {
            bool bound1 = false;
            bool bound2 = false;
            bool bound3 = false;
            bool bound4 = false;

            int i0 = tetrahedron[h].globalNode[0];
            int i1 = tetrahedron[h].globalNode[1];
            int i2 = tetrahedron[h].globalNode[2];
            int i3 = tetrahedron[h].globalNode[3];

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector1[coordinate] = (crossSection[ix].coordinate[coordinate] - node[i0].coordinate[coordinate]);
                vector2[coordinate] = (node[i1].coordinate[coordinate] - node[i0].coordinate[coordinate]);
                vector3[coordinate] = (node[i2].coordinate[coordinate] - node[i0].coordinate[coordinate]);
            };

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector[coordinate] = crossProduct(vector2, vector3, coordinate);
            };

            if ( dotProduct(vector1, vector) >= 0.0 )
                bound1 = true;

            /////////////////////////////////////////////////////////////////////////

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector1[coordinate] = (crossSection[ix].coordinate[coordinate] - node[i0].coordinate[coordinate]);
                vector2[coordinate] = (node[i2].coordinate[coordinate] - node[i0].coordinate[coordinate]);
                vector3[coordinate] = (node[i3].coordinate[coordinate] - node[i0].coordinate[coordinate]);
            };

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector[coordinate] = crossProduct(vector2, vector3, coordinate);
            };

            if ( dotProduct(vector1, vector) >= 0.0 )
                bound2 = true;


            /////////////////////////////////////////////////////////////////////////

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector1[coordinate] = (crossSection[ix].coordinate[coordinate] - node[i0].coordinate[coordinate]);
                vector2[coordinate] = (node[i3].coordinate[coordinate] - node[i0].coordinate[coordinate]);
                vector3[coordinate] = (node[i1].coordinate[coordinate] - node[i0].coordinate[coordinate]);
            };

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector[coordinate] = crossProduct(vector2, vector3, coordinate);
            };

            if ( dotProduct(vector1, vector) >= 0.0 )
                bound3 = true;


            /////////////////////////////////////////////////////////////////////////

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector1[coordinate] = (crossSection[ix].coordinate[coordinate] - node[i1].coordinate[coordinate]);
                vector2[coordinate] = (node[i3].coordinate[coordinate] - node[i1].coordinate[coordinate]);
                vector3[coordinate] = (node[i2].coordinate[coordinate] - node[i1].coordinate[coordinate]);
            };

            for (int coordinate = 0; coordinate < 3; coordinate++)
            {
                vector[coordinate] = crossProduct(vector2, vector3, coordinate);
            };

            if ( dotProduct(vector1, vector) >= 0.0 )
                bound4 = true;

            ///////////////////////////////////////////////////////////////////////

            if (bound1 && bound2 && bound3 && bound4)
            {
                crossSection[ix].tetrahedron = h;
                break;
            };
        };
        cout << ix << endl;
    };

    delete [] vector1; vector1 = NULL;
    delete [] vector2; vector2 = NULL;
    delete [] vector3; vector3 = NULL;
    delete [] vector; vector = NULL;
}

// for ALL
double Computation3D::globalTimeStep()
{
    double* tetrahedralTimeStep;
	tetrahedralTimeStep = new double [tetrahedronNumber];


	for (int h = 0; h < tetrahedronNumber; h++)
        tetrahedralTimeStep[h] = (2.0 / 3.0) * tetrahedron[h].volume / speed;


	double minimum = 10000.00;
	double CFL = 0.80;
	for (int h = 0; h < tetrahedronNumber; h++)
		if (tetrahedralTimeStep[h] < minimum)
			minimum = tetrahedralTimeStep[h];
	minimum = CFL * minimum;

	delete [] tetrahedralTimeStep; tetrahedralTimeStep = NULL;

	timeNumber = static_cast <int> (timeLast / minimum) + 1;

	return static_cast <double> (timeLast / timeNumber);
}

void Computation3D::GaussElimination(float** const & matrix, float** & inverse, const int & size) const
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

void Computation3D::matrixVectorMultiplication(float** const & matrix, float* const & vector, float* & result, const int & size) const
{
    for (int mi = 0; mi < size; mi++)
    {
        float sum = 0.0;
        for (int mj = 0; mj < size; mj++)
            sum += (matrix[mi][mj] * vector[mj]);
        result[mi] = sum;
    }
}

void Computation3D::findTetrahedronOrientation()
{
	for (int h = 0; h < tetrahedronNumber; h++)
	{
		double* positionVector1;
		positionVector1 = new double [3];
		for (int ix = 0; ix < 3; ix++)
		{
			positionVector1[ix] = node[tetrahedron[h].globalNode[1]].coordinate[ix] - node[tetrahedron[h].globalNode[0]].coordinate[ix];
		};

		double* positionVector2;
		positionVector2 = new double [3];
		for (int ix = 0; ix < 3; ix++)
		{
			positionVector2[ix] = node[tetrahedron[h].globalNode[2]].coordinate[ix] - node[tetrahedron[h].globalNode[0]].coordinate[ix];
		};

		double* crossPlane;
		crossPlane = new double [3];
		for (int ix = 0; ix < 3; ix++)
		{
			crossPlane[ix] = crossProduct(positionVector1, positionVector2, ix);
		};

		double* positionVector3;
		positionVector3 = new double [3];
		for (int ix = 0; ix < 3; ix++)
		{
			positionVector3[ix] = node[tetrahedron[h].globalNode[3]].coordinate[ix] - node[tetrahedron[h].globalNode[0]].coordinate[ix];
		};

		double dotPlane = dotProduct(crossPlane, positionVector3);

		if (dotPlane < 0)
			tetrahedron[h].orientation = 'A';
		else if (dotPlane >= 0)
			tetrahedron[h].orientation = 'B';

		delete [] positionVector1; positionVector1 = NULL;
		delete [] positionVector2; positionVector2 = NULL;
		delete [] crossPlane; crossPlane = NULL;
		delete [] positionVector3; positionVector3 = NULL;
	}
}

// for FINITE-volume
double Computation3D::LagrangeInterpolationFunction(double* const & median_Dual_Centroid, const int & H, const int & vertex) const
{
	double value;
	double x_0 = node[tetrahedron[H].globalNode[0]].coordinate[0];
	double y_0 = node[tetrahedron[H].globalNode[0]].coordinate[1];
	double z_0 = node[tetrahedron[H].globalNode[0]].coordinate[2];
	double x_1 = node[tetrahedron[H].globalNode[1]].coordinate[0];
	double y_1 = node[tetrahedron[H].globalNode[1]].coordinate[1];
	double z_1 = node[tetrahedron[H].globalNode[1]].coordinate[2];
	double x_2 = node[tetrahedron[H].globalNode[2]].coordinate[0];
	double y_2 = node[tetrahedron[H].globalNode[2]].coordinate[1];
	double z_2 = node[tetrahedron[H].globalNode[2]].coordinate[2];
	double x_3 = node[tetrahedron[H].globalNode[3]].coordinate[0];
	double y_3 = node[tetrahedron[H].globalNode[3]].coordinate[1];
	double z_3 = node[tetrahedron[H].globalNode[3]].coordinate[2];

	switch (tetrahedron[H].orientation)
	{
		case 'A':
			switch (vertex)
			{
				case 0:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_1) * ((y_2 - y_1) * (z_3 - z_1) - (y_3 - y_1) * (z_2 - z_1))
																	- (median_Dual_Centroid[1] - y_1) * ((x_2 - x_1) * (z_3 - z_1) - (x_3 - x_1) * (z_2 - z_1))
																	+ (median_Dual_Centroid[2] - z_1) * ((x_2 - x_1) * (y_3 - y_1) - (x_3 - x_1) * (y_2 - y_1)));
					break;
				case 1:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_0) * ((y_3 - y_0) * (z_2 - z_0) - (y_2 - y_0) * (z_3 - z_0))
																	- (median_Dual_Centroid[1] - y_0) * ((x_3 - x_0) * (z_2 - z_0) - (x_2 - x_0) * (z_3 - z_0))
																	+ (median_Dual_Centroid[2] - z_0) * ((x_3 - x_0) * (y_2 - y_0) - (x_2 - x_0) * (y_3 - y_0)));
					break;
				case 2:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_0) * ((y_1 - y_0) * (z_3 - z_0) - (y_3 - y_0) * (z_1 - z_0))
																	- (median_Dual_Centroid[1] - y_0) * ((x_1 - x_0) * (z_3 - z_0) - (x_3 - x_0) * (z_1 - z_0))
																	+ (median_Dual_Centroid[2] - z_0) * ((x_1 - x_0) * (y_3 - y_0) - (x_3 - x_0) * (y_1 - y_0)));
					break;
				case 3:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_0) * ((y_2 - y_0) * (z_1 - z_0) - (y_1 - y_0) * (z_2 - z_0))
																	- (median_Dual_Centroid[1] - y_0) * ((x_2 - x_0) * (z_1 - z_0) - (x_1 - x_0) * (z_2 - z_0))
																	+ (median_Dual_Centroid[2] - z_0) * ((x_2 - x_0) * (y_1 - y_0) - (x_1 - x_0) * (y_2 - y_0)));
					break;
			};
			break;
		case 'B':
			switch (vertex)
			{
				case 0:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_1) * ((y_3 - y_1) * (z_2 - z_1) - (y_2 - y_1) * (z_3 - z_1))
																	- (median_Dual_Centroid[1] - y_1) * ((x_3 - x_1) * (z_2 - z_1) - (x_2 - x_1) * (z_3 - z_1))
																	+ (median_Dual_Centroid[2] - z_1) * ((x_3 - x_1) * (y_2 - y_1) - (x_2 - x_1) * (y_3 - y_1)));
					break;
				case 1:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_0) * ((y_2 - y_0) * (z_3 - z_0) - (y_3 - y_0) * (z_2 - z_0))
																	- (median_Dual_Centroid[1] - y_0) * ((x_2 - x_0) * (z_3 - z_0) - (x_3 - x_0) * (z_2 - z_0))
																	+ (median_Dual_Centroid[2] - z_0) * ((x_2 - x_0) * (y_3 - y_0) - (x_3 - x_0) * (y_2 - y_0)));
					break;
				case 2:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_0) * ((y_3 - y_0) * (z_1 - z_0) - (y_1 - y_0) * (z_3 - z_0))
																	- (median_Dual_Centroid[1] - y_0) * ((x_3 - x_0) * (z_1 - z_0) - (x_1 - x_0) * (z_3 - z_0))
																	+ (median_Dual_Centroid[2] - z_0) * ((x_3 - x_0) * (y_1 - y_0) - (x_1 - x_0) * (y_3 - y_0)));
					break;
				case 3:
					value = (1.0 / (6.0 * tetrahedron[H].volume)) * ((median_Dual_Centroid[0] - x_0) * ((y_1 - y_0) * (z_2 - z_0) - (y_2 - y_0) * (z_1 - z_0))
																	- (median_Dual_Centroid[1] - y_0) * ((x_1 - x_0) * (z_2 - z_0) - (x_2 - x_0) * (z_1 - z_0))
																	+ (median_Dual_Centroid[2] - z_0) * ((x_1 - x_0) * (y_2 - y_0) - (x_2 - x_0) * (y_1 - y_0)));
					break;
			};
			break;
	}

  	return value;
}

// only for RD-LW
void Computation3D::constructDistributionMatrix()
{
	for (int h = 0; h < tetrahedronNumber; h++)
	{
		double speed = 1.0 / sqrt(permittivity * permeability);
		double impedance = sqrt(permeability / permittivity);

		double*** inflowMatrix;
		inflowMatrix = new double** [4];
		for (int j = 0; j < 4; j++)
		{
			inflowMatrix[j] = new double* [6];
			for (int ki = 0; ki < 6; ki++)
				inflowMatrix[j][ki] = new double [6];
		};

		for (int j = 0; j < 4; j++)
		{
			// Eigenvalues & Eigenvectors
			double xNormal = tetrahedron[h].inwardNormal[j][0];
			double yNormal = tetrahedron[h].inwardNormal[j][1];
			double zNormal = tetrahedron[h].inwardNormal[j][2];

			// INFLOW MATRIX
			inflowMatrix[j][0][0] = 0.0;
			inflowMatrix[j][0][1] = 0.0;
			inflowMatrix[j][0][2] = 0.0;
			inflowMatrix[j][0][3] = 0.0;
			inflowMatrix[j][0][4] = (1.0 / 3.0) * zNormal / permittivity;
			inflowMatrix[j][0][5] = - (1.0 / 3.0) * yNormal / permittivity;

			inflowMatrix[j][1][0] = 0.0;
			inflowMatrix[j][1][1] = 0.0;
			inflowMatrix[j][1][2] = 0.0;
			inflowMatrix[j][1][3] = - (1.0 / 3.0) * zNormal / permittivity;
			inflowMatrix[j][1][4] = 0.0;
			inflowMatrix[j][1][5] = (1.0 / 3.0) * xNormal / permittivity;

			inflowMatrix[j][2][0] = 0.0;
			inflowMatrix[j][2][1] = 0.0;
			inflowMatrix[j][2][2] = 0.0;
			inflowMatrix[j][2][3] = (1.0 / 3.0) * yNormal / permittivity;
			inflowMatrix[j][2][4] = - (1.0 / 3.0) * xNormal / permittivity;
			inflowMatrix[j][2][5] = 0.0;

			inflowMatrix[j][3][0] = 0.0;
			inflowMatrix[j][3][1] = - (1.0 / 3.0) * zNormal / permeability;
			inflowMatrix[j][3][2] = (1.0 / 3.0) * yNormal / permeability;
			inflowMatrix[j][3][3] = 0.0;
			inflowMatrix[j][3][4] = 0.0;
			inflowMatrix[j][3][5] = 0.0;

			inflowMatrix[j][4][0] = (1.0 / 3.0) * zNormal / permeability;
			inflowMatrix[j][4][1] = 0.0;
			inflowMatrix[j][4][2] = - (1.0 / 3.0) * xNormal / permeability;
			inflowMatrix[j][4][3] = 0.0;
			inflowMatrix[j][4][4] = 0.0;
			inflowMatrix[j][4][5] = 0.0;

			inflowMatrix[j][5][0] = - (1.0 / 3.0) * yNormal / permeability;
			inflowMatrix[j][5][1] = (1.0 / 3.0) * xNormal / permeability;
			inflowMatrix[j][5][2] = 0.0;
			inflowMatrix[j][5][3] = 0.0;
			inflowMatrix[j][5][4] = 0.0;
			inflowMatrix[j][5][5] = 0.0;
		};


		// Distribution Matrix
		for (int j = 0; j < 4; j++)
			for (int ki = 0; ki < 6; ki++)
				for (int kj = 0; kj < 6; kj++)
				{
					tetrahedron[h].distributionMatrix[j][ki][kj] = (1.0 / 4.0) * KroneckerDelta(ki, kj) + (timeDelta / (2.0 * tetrahedron[h].volume)) * inflowMatrix[j][ki][kj];
				};


		for (int j = 0; j < 4; j++)
		{
			for (int ki = 0; ki < 6; ki++)
			{
				delete [] inflowMatrix[j][ki]; inflowMatrix[j][ki] = NULL;
			};
			delete [] inflowMatrix[j]; inflowMatrix[j] = NULL;
		};
		delete [] inflowMatrix; inflowMatrix = NULL;
	}
}

void Computation3D::calculateFiniteVolume(const int & UVARIABLE_LEVEL)
{
	// int ipositive, inegative;
    double* fluxResidual;
    double* interpolateVariable;
    fluxResidual = new double [6];
	interpolateVariable = new double [6];
	double* displacement = new double [3];
	float* medianDualNormal = new float [3];
	double* medianDualCentroid = new double [3];

    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 6; ki++)
			node[i].fluxResidual[ki] = 0.0;

    for (int h = 0; h < tetrahedronNumber; h++)
    {
		for (int mediandual = 0; mediandual < 6; mediandual++)
        {
			int ipositive, inegative;
            switch (mediandual)
            {
                case 0:
                    ipositive = tetrahedron[h].globalNode[1];
                    inegative = tetrahedron[h].globalNode[0];
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0) / 4.0;
                    break;
                case 1:
                    ipositive = tetrahedron[h].globalNode[2];
                    inegative = tetrahedron[h].globalNode[0];
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0) / 4.0;
                    break;
                case 2:
                    ipositive = tetrahedron[h].globalNode[3];
                    inegative = tetrahedron[h].globalNode[0];
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0) / 4.0;
                    break;
                case 3:
                    ipositive = tetrahedron[h].globalNode[2];
                    inegative = tetrahedron[h].globalNode[1];
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0) / 4.0;
                    break;
                case 4:
                    ipositive = tetrahedron[h].globalNode[3];
                    inegative = tetrahedron[h].globalNode[2];
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0) / 4.0;
                    break;
                case 5:
                    ipositive = tetrahedron[h].globalNode[1];
                    inegative = tetrahedron[h].globalNode[3];
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0) / 4.0;
                    break;
            };

			////////////////////////////////////////////////////////////////////////////////////////////////
			// positive side of median dual normal
			for (int coor = 0; coor < 3; coor++)
			{
				displacement[coor] = medianDualCentroid[coor] - node[ipositive].coordinate[coor];
				medianDualNormal[coor] = tetrahedron[h].medianDualNormal[mediandual][coor];
			};

			for (int ki = 0; ki < 6; ki++)
			{
				interpolateVariable[ki] = node[ipositive].conservedVariable[ki][UVARIABLE_LEVEL] 
											+ dotProduct(displacement, node[ipositive].firstDerivative[ki]);
			};

            fluxResidual[0] = (interpolateVariable[4] * medianDualNormal[2] / permittivity - interpolateVariable[5] * medianDualNormal[1] / permittivity);
            fluxResidual[1] = (interpolateVariable[5] * medianDualNormal[0] / permittivity - interpolateVariable[3] * medianDualNormal[2] / permittivity);
            fluxResidual[2] = (interpolateVariable[3] * medianDualNormal[1] / permittivity - interpolateVariable[4] * medianDualNormal[0] / permittivity);
            fluxResidual[3] = (interpolateVariable[2] * medianDualNormal[1] / permeability - interpolateVariable[1] * medianDualNormal[2] / permeability);
            fluxResidual[4] = (interpolateVariable[0] * medianDualNormal[2] / permeability - interpolateVariable[2] * medianDualNormal[0] / permeability);
            fluxResidual[5] = (interpolateVariable[1] * medianDualNormal[0] / permeability - interpolateVariable[0] * medianDualNormal[1] / permeability);

			// flux out from local node ipositive
            for (int ki = 0; ki < 6; ki++)
                node[ipositive].fluxResidual[ki] = node[ipositive].fluxResidual[ki] + 0.5 * fluxResidual[ki];

			// flux into local node inegative
			for (int ki = 0; ki < 6; ki++)
                node[inegative].fluxResidual[ki] = node[inegative].fluxResidual[ki] - 0.5 * fluxResidual[ki];
			
			////////////////////////////////////////////////////////////////////////////////////////////////
			// negative side of median dual normal
			for (int coor = 0; coor < 3; coor++)
			{
				displacement[coor] = medianDualCentroid[coor] - node[inegative].coordinate[coor];
				medianDualNormal[coor] = - tetrahedron[h].medianDualNormal[mediandual][coor];
			};

			for (int ki = 0; ki < 6; ki++)
			{
				interpolateVariable[ki] = node[inegative].conservedVariable[ki][UVARIABLE_LEVEL] 
											+ dotProduct(displacement, node[inegative].firstDerivative[ki]);
			};

            fluxResidual[0] = (interpolateVariable[4] * medianDualNormal[2] / permittivity - interpolateVariable[5] * medianDualNormal[1] / permittivity);
            fluxResidual[1] = (interpolateVariable[5] * medianDualNormal[0] / permittivity - interpolateVariable[3] * medianDualNormal[2] / permittivity);
            fluxResidual[2] = (interpolateVariable[3] * medianDualNormal[1] / permittivity - interpolateVariable[4] * medianDualNormal[0] / permittivity);
            fluxResidual[3] = (interpolateVariable[2] * medianDualNormal[1] / permeability - interpolateVariable[1] * medianDualNormal[2] / permeability);
            fluxResidual[4] = (interpolateVariable[0] * medianDualNormal[2] / permeability - interpolateVariable[2] * medianDualNormal[0] / permeability);
            fluxResidual[5] = (interpolateVariable[1] * medianDualNormal[0] / permeability - interpolateVariable[0] * medianDualNormal[1] / permeability);

			// flux into local node ipositive
            for (int ki = 0; ki < 6; ki++)
                node[ipositive].fluxResidual[ki] = node[ipositive].fluxResidual[ki] - 0.5 * fluxResidual[ki];
			
			// flux out from local node inegative
			for (int ki = 0; ki < 6; ki++)
                node[inegative].fluxResidual[ki] = node[inegative].fluxResidual[ki] + 0.5 * fluxResidual[ki];
        };
    };

	delete [] fluxResidual; fluxResidual = NULL;
	delete [] interpolateVariable; interpolateVariable = NULL;
	delete [] displacement; displacement = NULL;
	delete [] medianDualNormal; medianDualNormal = NULL;
	delete [] medianDualCentroid; medianDualCentroid = NULL;
}

void Computation3D::calculateGradient(const int & UVARIABLE_LEVEL)
{
	if (method == 'A')
	{
		for (int i = 0; i < nodeNumber; i++)
    	{
			float Sxx = 0.0;
			float Syy = 0.0;
			float Szz = 0.0;
			float Sxy = 0.0;
			float Sxz = 0.0;
			float Syz = 0.0;
			float Sxu;
			float Syu;
			float Szu;

			float* xk = new float [node[i].neighbourNode.size()];
			float* yk = new float [node[i].neighbourNode.size()];
			float* zk = new float [node[i].neighbourNode.size()];

			float** matrix = new float* [3];
			float** inverse = new float* [3];
			for (int mi = 0; mi < 3; mi++)
			{
				matrix[mi] = new float [3];
				inverse[mi] = new float [3];
			};

			for (int n = 0; n < node[i].neighbourNode.size(); n++)
			{
				xk[n] = node[node[i].neighbourNode[n]].coordinate[0] - node[i].coordinate[0];
				yk[n] = node[node[i].neighbourNode[n]].coordinate[1] - node[i].coordinate[1];
				zk[n] = node[node[i].neighbourNode[n]].coordinate[2] - node[i].coordinate[2];
			};

			// calculations for Sxx, Sxy, Syy
			for (int n = 0; n < node[i].neighbourNode.size(); n++)
			{
				if (node[i].neighbourNode[n] != i)
				{
					Sxx += (xk[n] * xk[n]);
					Syy += (yk[n] * yk[n]);
					Szz += (zk[n] * zk[n]);
					Sxy += (xk[n] * yk[n]);
					Sxz += (xk[n] * zk[n]);
					Syz += (yk[n] * zk[n]);
				};
			};

			matrix[0][0] = Sxx;
			matrix[0][1] = Sxy;
			matrix[0][2] = Sxz;

			matrix[1][0] = Sxy;
			matrix[1][1] = Syy;
			matrix[1][2] = Syz;

			matrix[2][0] = Sxz;
			matrix[2][1] = Syz;
			matrix[2][2] = Szz;

			GaussElimination(matrix, inverse, 3);

			// calculations for first derivative
			float* vector = new float [3];
			float* result = new float [3];

			for (int ki = 0; ki < 6; ki++)
			{
				Sxu = 0.0;
				Syu = 0.0;
				Szu = 0.0;
			
				for (int n = 0; n < node[i].neighbourNode.size(); n++)
				{
					float uk = node[node[i].neighbourNode[n]].conservedVariable[ki][UVARIABLE_LEVEL] - node[i].conservedVariable[ki][UVARIABLE_LEVEL];

					if (node[i].neighbourNode[n] != i)
					{
						Sxu += (xk[n] * uk);
						Syu += (yk[n] * uk);
						Szu += (zk[n] * uk);
					};
				};

				vector[0] = Sxu;
				vector[1] = Syu;
				vector[2] = Szu;

				matrixVectorMultiplication(inverse, vector, result, 3);

				// inverse estimator matrix mulply with Sxu, Syu, Szu
				// node[i].firstDerivative[ki][0] = result[0];
				// node[i].firstDerivative[ki][1] = result[1];
				// node[i].firstDerivative[ki][2] = result[2];

				// gradient limiter
				if (ki == 2 || ki == 4)
				{
					node[i].firstDerivative[ki][0] = result[0];
				}
				else
				{
					node[i].firstDerivative[ki][0] = 0.0;
				}
				node[i].firstDerivative[ki][1] = 0.0;
				node[i].firstDerivative[ki][2] = 0.0;
			};

			for (int mi = 0; mi < 3; mi++)
			{
				delete [] matrix[mi]; matrix[mi] = NULL;
				delete [] inverse[mi]; inverse[mi] = NULL;
			};
			delete [] matrix; matrix = NULL;
			delete [] inverse; inverse = NULL;
			delete [] vector; vector = NULL;
			delete [] result; result = NULL;

			delete [] xk; xk = NULL;
			delete [] yk; yk = NULL;
			delete [] zk; zk = NULL;

			///////////////////////////////////////////////////////////////////////////////////
			// analytical first derivative
			// complexNumber Hyx, Ezx;

			// Ezx.real = propagationCoefficient * sin(- propagationCoefficient * node[i].coordinate[0]);
			// Ezx.imaginary = - propagationCoefficient * cos(- propagationCoefficient * node[i].coordinate[0]);

			// Hyx.real = (- pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * sin(- propagationCoefficient * node[i].coordinate[0]);
			// Hyx.imaginary = (pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * cos(- propagationCoefficient * node[i].coordinate[0]);

			// complexNumber timeHarmonic;
			// timeHarmonic.real = cos(angularFrequency * time);
			// timeHarmonic.imaginary = sin(angularFrequency * time);

			// for (int ki = 0; ki < 6; ki++)
			// 	for (int coor = 0; coor < 3; coor++)
			// 		node[i].firstDerivative[ki][coor] = 0.0;

			// node[i].firstDerivative[2][0] = complexMultiplication(timeHarmonic, Ezx).real;
			// node[i].firstDerivative[4][0] = complexMultiplication(timeHarmonic, Hyx).real;

		};
	}
	else if (method == 'B')
	{
		for (int h = 0; h < tetrahedronNumber; h++)
    	{
			if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
			{
			   	for (int ki = 0; ki < 6; ki++)
        		{
            		tetrahedron[h].gradient[ki][0] = 0.0;
            		tetrahedron[h].gradient[ki][1] = 0.0;
					tetrahedron[h].gradient[ki][2] = 0.0;
            		for (int j = 0; j < 4; j++)
            		{
                		int i = tetrahedron[h].globalNode[j];
                		tetrahedron[h].gradient[ki][0] = tetrahedron[h].gradient[ki][0] 
															+ (1.0 / 3.0) * (tetrahedron[h].inwardNormal[j][0] / tetrahedron[h].volume) * node[i].conservedVariable[ki][UVARIABLE_LEVEL];
                		tetrahedron[h].gradient[ki][1] = tetrahedron[h].gradient[ki][1] 
															+ (1.0 / 3.0) * (tetrahedron[h].inwardNormal[j][1] / tetrahedron[h].volume) * node[i].conservedVariable[ki][UVARIABLE_LEVEL];
						tetrahedron[h].gradient[ki][2] = tetrahedron[h].gradient[ki][2] 
															+ (1.0 / 3.0) * (tetrahedron[h].inwardNormal[j][2] / tetrahedron[h].volume) * node[i].conservedVariable[ki][UVARIABLE_LEVEL];
            		};
        		}
			}
    	};

		// for (int i = 0; i < nodeNumber; i++)
    	// {
		// 	float Sxx = 0.0;
		// 	float Syy = 0.0;
		// 	float Szz = 0.0;
		// 	float Sxy = 0.0;
		// 	float Sxz = 0.0;
		// 	float Syz = 0.0;
		// 	float Sxu;
		// 	float Syu;
		// 	float Szu;

		// 	float* xk = new float [node[i].neighbourNode.size()];
		// 	float* yk = new float [node[i].neighbourNode.size()];
		// 	float* zk = new float [node[i].neighbourNode.size()];

		// 	float** matrix = new float* [3];
		// 	float** inverse = new float* [3];
		// 	for (int mi = 0; mi < 3; mi++)
		// 	{
		// 		matrix[mi] = new float [3];
		// 		inverse[mi] = new float [3];
		// 	};

		// 	for (int n = 0; n < node[i].neighbourNode.size(); n++)
		// 	{
		// 		xk[n] = node[node[i].neighbourNode[n]].coordinate[0] - node[i].coordinate[0];
		// 		yk[n] = node[node[i].neighbourNode[n]].coordinate[1] - node[i].coordinate[1];
		// 		zk[n] = node[node[i].neighbourNode[n]].coordinate[2] - node[i].coordinate[2];
		// 	};

		// 	// calculations for Sxx, Sxy, Syy
		// 	for (int n = 0; n < node[i].neighbourNode.size(); n++)
		// 	{
		// 		if (node[i].neighbourNode[n] != i)
		// 		{
		// 			Sxx += (xk[n] * xk[n]);
		// 			Syy += (yk[n] * yk[n]);
		// 			Szz += (zk[n] * zk[n]);
		// 			Sxy += (xk[n] * yk[n]);
		// 			Sxz += (xk[n] * zk[n]);
		// 			Syz += (yk[n] * zk[n]);
		// 		};
		// 	};

		// 	matrix[0][0] = Sxx;
		// 	matrix[0][1] = Sxy;
		// 	matrix[0][2] = Sxz;

		// 	matrix[1][0] = Sxy;
		// 	matrix[1][1] = Syy;
		// 	matrix[1][2] = Syz;

		// 	matrix[2][0] = Sxz;
		// 	matrix[2][1] = Syz;
		// 	matrix[2][2] = Szz;

		// 	GaussElimination(matrix, inverse, 3);

		// 	// calculations for first derivative
		// 	float* vector = new float [3];
		// 	float* result = new float [3];

		// 	for (int ki = 0; ki < 6; ki++)
		// 	{
		// 		Sxu = 0.0;
		// 		Syu = 0.0;
		// 		Szu = 0.0;
			
		// 		for (int n = 0; n < node[i].neighbourNode.size(); n++)
		// 		{
		// 			float uk = node[node[i].neighbourNode[n]].conservedVariable[ki][UVARIABLE_LEVEL] - node[i].conservedVariable[ki][UVARIABLE_LEVEL];

		// 			if (node[i].neighbourNode[n] != i)
		// 			{
		// 				Sxu += (xk[n] * uk);
		// 				Syu += (yk[n] * uk);
		// 				Szu += (zk[n] * uk);
		// 			};
		// 		};

		// 		vector[0] = Sxu;
		// 		vector[1] = Syu;
		// 		vector[2] = Szu;

		// 		matrixVectorMultiplication(inverse, vector, result, 3);

		// 		// inverse estimator matrix mulply with Sxu, Syu, Szu
		// 		node[i].firstDerivative[ki][0] = result[0];
		// 		node[i].firstDerivative[ki][1] = result[1];
		// 		node[i].firstDerivative[ki][2] = result[2];
		// 	};

		// 	for (int mi = 0; mi < 3; mi++)
		// 	{
		// 		delete [] matrix[mi]; matrix[mi] = NULL;
		// 		delete [] inverse[mi]; inverse[mi] = NULL;
		// 	};
		// 	delete [] matrix; matrix = NULL;
		// 	delete [] inverse; inverse = NULL;
		// 	delete [] vector; vector = NULL;
		// 	delete [] result; result = NULL;

		// 	delete [] xk; xk = NULL;
		// 	delete [] yk; yk = NULL;
		// 	delete [] zk; zk = NULL;

		// 	for (int ki = 0; ki < 6; ki++)
        // 	    for (int coor = 0; coor < 3; coor++)
        // 	    {
        // 	        node[i].firstDerivative[ki][coor] = 0.0; 

        // 	        for (int n = 0; n < node[i].neighbourTetrahedron.size(); n++)
        // 	        {
        // 	            node[i].firstDerivative[ki][coor] = node[i].firstDerivative[ki][coor]
        // 	                                                    + tetrahedron[node[i].neighbourTetrahedron[n]].gradient[ki][coor];
        // 	        }

        // 	        if (node[i].neighbourTetrahedron.size() != 0)
        // 	            node[i].firstDerivative[ki][coor] = node[i].firstDerivative[ki][coor] / (node[i].neighbourTetrahedron.size());
        // 	    };

		// 	///////////////////////////////////////////////////////////////////////////////////
		// 	// analytical first derivative
		// 	complexNumber Hyx, Ezx;

		// 	Ezx.real = propagationCoefficient * sin(- propagationCoefficient * node[i].coordinate[0]);
		// 	Ezx.imaginary = - propagationCoefficient * cos(- propagationCoefficient * node[i].coordinate[0]);

		// 	Hyx.real = (- pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * sin(- propagationCoefficient * node[i].coordinate[0]);
		// 	Hyx.imaginary = (pow(propagationCoefficient, 2) / (angularFrequency * permeability)) * cos(- propagationCoefficient * node[i].coordinate[0]);

		// 	complexNumber timeHarmonic;
		// 	timeHarmonic.real = cos(angularFrequency * time);
		// 	timeHarmonic.imaginary = sin(angularFrequency * time);

		// 	for (int ki = 0; ki < 6; ki++)
		// 		for (int coor = 0; coor < 3; coor++)
		// 			node[i].firstDerivative[ki][coor] = 0.0;

		// 	node[i].firstDerivative[2][0] = complexMultiplication(timeHarmonic, Ezx).real;
		// 	node[i].firstDerivative[4][0] = complexMultiplication(timeHarmonic, Hyx).real;
		// };
	}
}

void Computation3D::fluxDifference(const int & H, const int & UVARIABLE_LEVEL)
{
	///////////////////////////////////////////////////////////
    // calculate flux residual
	double xElectricField = 0.0;
    double yElectricField = 0.0;
    double zElectricField = 0.0;
    double xMagneticField = 0.0;
    double yMagneticField = 0.0;
    double zMagneticField = 0.0;

	for (int j = 0; j < 4; j++)
	{
		int i = tetrahedron[H].globalNode[j];

		xElectricField = xElectricField + node[i].conservedVariable[0][UVARIABLE_LEVEL];
		yElectricField = yElectricField + node[i].conservedVariable[1][UVARIABLE_LEVEL];
		zElectricField = zElectricField + node[i].conservedVariable[2][UVARIABLE_LEVEL];
		xMagneticField = xMagneticField + node[i].conservedVariable[3][UVARIABLE_LEVEL];
		yMagneticField = yMagneticField + node[i].conservedVariable[4][UVARIABLE_LEVEL];
		zMagneticField = zMagneticField + node[i].conservedVariable[5][UVARIABLE_LEVEL];
	};
	xElectricField = xElectricField / 4.0;
	yElectricField = yElectricField / 4.0;
	zElectricField = zElectricField / 4.0;
	xMagneticField = xMagneticField / 4.0;
	yMagneticField = yMagneticField / 4.0;
	zMagneticField = zMagneticField / 4.0;

	for (int j = 0; j < 4; j++)
	{
		int i = tetrahedron[H].globalNode[j];

		double xMedianDualNormal = - (1.0 / 3.0) * tetrahedron[H].inwardNormal[j][0];
        double yMedianDualNormal = - (1.0 / 3.0) * tetrahedron[H].inwardNormal[j][1];
        double zMedianDualNormal = - (1.0 / 3.0) * tetrahedron[H].inwardNormal[j][2];

		node[i].fluxResidual[0] = node[i].fluxResidual[0] + ((yMagneticField / permittivity) * zMedianDualNormal - (zMagneticField / permittivity) * yMedianDualNormal);
	    node[i].fluxResidual[1] = node[i].fluxResidual[1] + ((zMagneticField / permittivity) * xMedianDualNormal - (xMagneticField / permittivity) * zMedianDualNormal);
	    node[i].fluxResidual[2] = node[i].fluxResidual[2] + ((xMagneticField / permittivity) * yMedianDualNormal - (yMagneticField / permittivity) * xMedianDualNormal);
	    node[i].fluxResidual[3] = node[i].fluxResidual[3] + ((zElectricField / permittivity) * yMedianDualNormal - (yElectricField / permittivity) * zMedianDualNormal);
	    node[i].fluxResidual[4] = node[i].fluxResidual[4] + ((xElectricField / permittivity) * zMedianDualNormal - (zElectricField / permittivity) * xMedianDualNormal);
	    node[i].fluxResidual[5] = node[i].fluxResidual[5] + ((yElectricField / permittivity) * xMedianDualNormal - (xElectricField / permittivity) * yMedianDualNormal);
	};

	// ///////////////////////////////////////////////////////////
	// // calculate numerical dissipation
	// double* secondDerivative = new double [6];
	// for (int j = 0; j < 4; j++)
	// {
	// 	int i = tetrahedron[H].globalNode[j];
	// 	int i1 = tetrahedron[H].globalNode[(j + 1) % 4];
	// 	int i2 = tetrahedron[H].globalNode[(j + 2) % 4];
	// 	int i3 = tetrahedron[H].globalNode[(j + 3) % 4];

	// 	for (int ki = 0; ki < 6; ki++)
	// 	{
	// 		secondDerivative[ki] = (1.0 / 4.0) * (0.5 * pow(node[i1].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
	// 												+ 0.5 * (node[i1].coordinate[0] - node[i].coordinate[0]) * (node[i1].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
	// 												+ 0.5 * (node[i1].coordinate[0] - node[i].coordinate[0]) * (node[i1].coordinate[2] - node[i].coordinate[2]) * node[i].secondDerivative[ki][0][2]
	// 												+ 0.5 * (node[i1].coordinate[1] - node[i].coordinate[1]) * (node[i1].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
	// 												+ 0.5 * pow(node[i1].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]
	// 												+ 0.5 * (node[i1].coordinate[1] - node[i].coordinate[1]) * (node[i1].coordinate[2] - node[i].coordinate[2]) * node[i].secondDerivative[ki][1][2]
	// 												+ 0.5 * (node[i1].coordinate[2] - node[i].coordinate[2]) * (node[i1].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][2][0]
	// 												+ 0.5 * (node[i1].coordinate[2] - node[i].coordinate[2]) * (node[i1].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][2][1]
	// 												+ 0.5 * pow(node[i1].coordinate[2] - node[i].coordinate[2], 2) * node[i].secondDerivative[ki][2][2]

	// 												+ 0.5 * pow(node[i2].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
	// 												+ 0.5 * (node[i2].coordinate[0] - node[i].coordinate[0]) * (node[i2].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
	// 												+ 0.5 * (node[i2].coordinate[0] - node[i].coordinate[0]) * (node[i2].coordinate[2] - node[i].coordinate[2]) * node[i].secondDerivative[ki][0][2]
	// 												+ 0.5 * (node[i2].coordinate[1] - node[i].coordinate[1]) * (node[i2].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
	// 												+ 0.5 * pow(node[i2].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]
	// 												+ 0.5 * (node[i2].coordinate[1] - node[i].coordinate[1]) * (node[i2].coordinate[2] - node[i].coordinate[2]) * node[i].secondDerivative[ki][1][2]
	// 												+ 0.5 * (node[i2].coordinate[2] - node[i].coordinate[2]) * (node[i2].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][2][0]
	// 												+ 0.5 * (node[i2].coordinate[2] - node[i].coordinate[2]) * (node[i2].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][2][1]
	// 												+ 0.5 * pow(node[i2].coordinate[2] - node[i].coordinate[2], 2) * node[i].secondDerivative[ki][2][2]
									
	// 												+ 0.5 * pow(node[i3].coordinate[0] - node[i].coordinate[0], 2) * node[i].secondDerivative[ki][0][0]
	// 												+ 0.5 * (node[i3].coordinate[0] - node[i].coordinate[0]) * (node[i3].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][0][1]
	// 												+ 0.5 * (node[i3].coordinate[0] - node[i].coordinate[0]) * (node[i3].coordinate[2] - node[i].coordinate[2]) * node[i].secondDerivative[ki][0][2]
	// 												+ 0.5 * (node[i3].coordinate[1] - node[i].coordinate[1]) * (node[i3].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][1][0]
	// 												+ 0.5 * pow(node[i3].coordinate[1] - node[i].coordinate[1], 2) * node[i].secondDerivative[ki][1][1]
	// 												+ 0.5 * (node[i3].coordinate[1] - node[i].coordinate[1]) * (node[i3].coordinate[2] - node[i].coordinate[2]) * node[i].secondDerivative[ki][1][2]
	// 												+ 0.5 * (node[i3].coordinate[2] - node[i].coordinate[2]) * (node[i3].coordinate[0] - node[i].coordinate[0]) * node[i].secondDerivative[ki][2][0]
	// 												+ 0.5 * (node[i3].coordinate[2] - node[i].coordinate[2]) * (node[i3].coordinate[1] - node[i].coordinate[1]) * node[i].secondDerivative[ki][2][1]
	// 												+ 0.5 * pow(node[i3].coordinate[2] - node[i].coordinate[2], 2) * node[i].secondDerivative[ki][2][2]);
	// 	};

	// 	double xMedianDualNormal = - (1.0 / 3.0) * tetrahedron[H].inwardNormal[j][0];
    //     double yMedianDualNormal = - (1.0 / 3.0) * tetrahedron[H].inwardNormal[j][1];
    //     double zMedianDualNormal = - (1.0 / 3.0) * tetrahedron[H].inwardNormal[j][2];

	// 	node[i].dissipation[0] = node[i].dissipation[0] + ((secondDerivative[4] / permittivity) * zMedianDualNormal - (secondDerivative[5] / permittivity) * yMedianDualNormal);
	//     node[i].dissipation[1] = node[i].dissipation[1] + ((secondDerivative[5] / permittivity) * xMedianDualNormal - (secondDerivative[3] / permittivity) * zMedianDualNormal);
	//     node[i].dissipation[2] = node[i].dissipation[2] + ((secondDerivative[3] / permittivity) * yMedianDualNormal - (secondDerivative[4] / permittivity) * xMedianDualNormal);
	//     node[i].dissipation[3] = node[i].dissipation[3] + ((secondDerivative[2] / permittivity) * yMedianDualNormal - (secondDerivative[1] / permittivity) * zMedianDualNormal);
	//     node[i].dissipation[4] = node[i].dissipation[4] + ((secondDerivative[0] / permittivity) * zMedianDualNormal - (secondDerivative[2] / permittivity) * xMedianDualNormal);
	//     node[i].dissipation[5] = node[i].dissipation[5] + ((secondDerivative[1] / permittivity) * xMedianDualNormal - (secondDerivative[0] / permittivity) * yMedianDualNormal);
	// };
	// delete [] secondDerivative; secondDerivative = NULL;
}

void Computation3D::calculateFluxDifference(const int & UVARIABLE_LEVEL)
{
    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 6; ki++)
		{
			node[i].fluxResidual[ki] = 0.0;
			node[i].dissipation[ki] = 0.0;
		}
	
    for (int h = 0; h < tetrahedronNumber; h++)
	{
        fluxDifference(h, UVARIABLE_LEVEL);
	};
}

void Computation3D::boundaryFluxResidual(double* & flux_Residual, double* const & interpolate_Variable,
                                        const double & x_Normal, const double & y_Normal, const double & z_Normal) const
{
    flux_Residual[0] = ((interpolate_Variable[4] / permittivity) * z_Normal - (interpolate_Variable[5] / permittivity) * y_Normal);
    flux_Residual[1] = ((interpolate_Variable[5] / permittivity) * x_Normal - (interpolate_Variable[3] / permittivity) * z_Normal);
    flux_Residual[2] = ((interpolate_Variable[3] / permittivity) * y_Normal - (interpolate_Variable[4] / permittivity) * x_Normal);
    flux_Residual[3] = ((interpolate_Variable[2] / permeability) * y_Normal - (interpolate_Variable[1] / permeability) * z_Normal);
    flux_Residual[4] = ((interpolate_Variable[0] / permeability) * z_Normal - (interpolate_Variable[2] / permeability) * x_Normal);
    flux_Residual[5] = ((interpolate_Variable[1] / permeability) * x_Normal - (interpolate_Variable[0] / permeability) * y_Normal);
}

void Computation3D::calculateFiniteVolumeBoundaryFlux(const int & UVARIABLE_LEVEL)
{
    double* fluxResidual;
    fluxResidual = new double [6];

    double* coordinate;
    coordinate = new double [3];

    double* centroid;
    centroid = new double [3];

    double* interpolateVariable;
    interpolateVariable = new double [6];

	double* displacement;
	displacement = new double [3];

	double* firstDerivative;
	firstDerivative = new double [6];

    switch (TMTEMode)
    {
      case 'A':
      case 'B':
        for (int h = 0; h < tetrahedronNumber; h++)
        {
          	if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
          	{
              	for (int ki = 0; ki < 6; ki++)
                	fluxResidual[ki] = 0.0;

              	int j = tetrahedron[h].boundaryVertex;
				int i1 = tetrahedron[h].globalNode[(j + 1) % 4];
				int i2 = tetrahedron[h].globalNode[(j + 2) % 4];
				int i3 = tetrahedron[h].globalNode[(j + 3) % 4];

              	double xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][0];
              	double yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][1];
              	double zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][2];

              	for (int ix = 0; ix < 3; ix++)
                	centroid[ix] = (node[i1].coordinate[ix]
                                    + node[i2].coordinate[ix]
                                    + node[i3].coordinate[ix]) / 3.0;

              	/////////////////////////////////////////////////////
              	// (j + 1) % 4
              	for (int ix = 0; ix < 3; ix++)
                	coordinate[ix] = (node[i1].coordinate[ix]
                    	                + (node[i1].coordinate[ix] + node[i2].coordinate[ix]) / 2.0
                                        + (node[i1].coordinate[ix] + node[i3].coordinate[ix]) / 2.0
                                        + centroid[ix]) / 4.0;

				for (int ix = 0; ix < 3; ix++)
					displacement[ix] = coordinate[ix] - node[tetrahedron[h].globalNode[j]].coordinate[ix];

				for (int ki = 0; ki < 6; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[tetrahedron[h].globalNode[j]].firstDerivative[ki]);

              	for (int ki = 0; ki < 6; ki++)
                	interpolateVariable[ki] = node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
												+ firstDerivative[ki];

              	boundaryFluxResidual(fluxResidual, interpolateVariable,
                	                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              	if (tetrahedron[h].boundaryType == 'B')
              	{
                	for (int ki = 3; ki < 6; ki++)
                    	node[i1].fluxResidual[ki] = node[i1].fluxResidual[ki] + fluxResidual[ki];
              	}
              	else if (tetrahedron[h].boundaryType == 'C')
              	{
                	for (int ki = 0; ki < 6; ki++)
                    	node[i1].fluxResidual[ki] = node[i1].fluxResidual[ki] + fluxResidual[ki];
              	};

              	/////////////////////////////////////////////////////
              	// (j + 2) % 4
              	for (int ix = 0; ix < 3; ix++)
                	coordinate[ix] = (node[i2].coordinate[ix]
                                        + (node[i2].coordinate[ix] + node[i3].coordinate[ix]) / 2.0
                                        + (node[i2].coordinate[ix] + node[i1].coordinate[ix]) / 2.0
                                        + centroid[ix]) / 4.0;

				for (int ix = 0; ix < 3; ix++)
					displacement[ix] = coordinate[ix] - node[tetrahedron[h].globalNode[j]].coordinate[ix];

				for (int ki = 0; ki < 6; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[tetrahedron[h].globalNode[j]].firstDerivative[ki]);

              	for (int ki = 0; ki < 6; ki++)
                	interpolateVariable[ki] = node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
												+ firstDerivative[ki];

              	boundaryFluxResidual(fluxResidual, interpolateVariable,
                	                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              	if (tetrahedron[h].boundaryType == 'B')
              	{
                	for (int ki = 3; ki < 6; ki++)
                    	node[i2].fluxResidual[ki] = node[i2].fluxResidual[ki] + fluxResidual[ki];
              	}
              	else if (tetrahedron[h].boundaryType == 'C')
              	{
                	for (int ki = 0; ki < 6; ki++)
                    	node[i2].fluxResidual[ki] = node[i2].fluxResidual[ki] + fluxResidual[ki];
              	};

              	/////////////////////////////////////////////////////
              	// (j + 3) % 4
              	for (int ix = 0; ix < 3; ix++)
                	coordinate[ix] = (node[i3].coordinate[ix]
                                        + (node[i3].coordinate[ix] + node[i1].coordinate[ix]) / 2.0
                                        + (node[i3].coordinate[ix] + node[i2].coordinate[ix]) / 2.0
                                        + centroid[ix]) / 4.0;

				for (int ix = 0; ix < 3; ix++)
					displacement[ix] = coordinate[ix] - node[tetrahedron[h].globalNode[j]].coordinate[ix];

				for (int ki = 0; ki < 6; ki++)
					firstDerivative[ki] = dotProduct(displacement, node[tetrahedron[h].globalNode[j]].firstDerivative[ki]);

				for (int ki = 0; ki < 6; ki++)
                  	interpolateVariable[ki] = node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
												+ firstDerivative[ki];

              	boundaryFluxResidual(fluxResidual, interpolateVariable,
                	                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              	if (tetrahedron[h].boundaryType == 'B')
              	{
                	for (int ki = 3; ki < 6; ki++)
                    	node[i3].fluxResidual[ki] = node[i3].fluxResidual[ki] + fluxResidual[ki];
              	}
              	else if (tetrahedron[h].boundaryType == 'C')
              	{
                	for (int ki = 0; ki < 6; ki++)
                    	node[i3].fluxResidual[ki] = node[i3].fluxResidual[ki] + fluxResidual[ki];
              	};
          	}
        };
        break;
    };
    delete [] fluxResidual; fluxResidual = NULL;
    delete [] coordinate; coordinate = NULL;
    delete [] centroid; centroid = NULL;
    delete [] interpolateVariable; interpolateVariable = NULL;
	delete [] displacement; displacement = NULL;
	delete [] firstDerivative; firstDerivative = NULL;
}

void Computation3D::calculateFluxDifferenceBoundaryFlux(const int & UVARIABLE_LEVEL)
{
    double* fluxResidual;
    fluxResidual = new double [6];

    double* coordinate;
    coordinate = new double [3];

    double* centroid;
    centroid = new double [3];

    double* interpolateVariable;
    interpolateVariable = new double [6];

	double* displacement;
	displacement = new double [3];

	double* firstDerivative;
	firstDerivative = new double [6];

    switch (TMTEMode)
    {
    	case 'A':
      	case 'B':
        	for (int h = 0; h < tetrahedronNumber; h++)
        	{
          		if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
          		{
              		for (int ki = 0; ki < 6; ki++)
                  		fluxResidual[ki] = 0.0;

              		int j = tetrahedron[h].boundaryVertex;
					int i1 = tetrahedron[h].globalNode[(j + 1) % 4];
					int i2 = tetrahedron[h].globalNode[(j + 2) % 4];
					int i3 = tetrahedron[h].globalNode[(j + 3) % 4];

              		double xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][0];
              		double yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][1];
              		double zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][2];

              		for (int ix = 0; ix < 3; ix++)
                  		centroid[ix] = (node[i1].coordinate[ix] + node[i2].coordinate[ix] + node[i3].coordinate[ix]) / 3.0;

              		/////////////////////////////////////////////////////
              		// (j + 1) % 4
              		for (int ix = 0; ix < 3; ix++)
                  		coordinate[ix] = (node[i1].coordinate[ix]
                                              + (node[i1].coordinate[ix] + node[i2].coordinate[ix]) / 2.0
                                              + (node[i1].coordinate[ix] + node[i3].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

					for (int ix = 0; ix < 3; ix++)
						displacement[ix] = coordinate[ix] - node[tetrahedron[h].globalNode[j]].coordinate[ix];

					for (int ki = 0; ki < 6; ki++)
						firstDerivative[ki] = dotProduct(displacement, tetrahedron[h].gradient[ki]);

              		for (int ki = 0; ki < 6; ki++)
                  		interpolateVariable[ki] = node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
													+ firstDerivative[ki];

              		boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    		xOutwardNormal, yOutwardNormal, zOutwardNormal);

              		if (tetrahedron[h].boundaryType == 'B')
              		{
                  		for (int ki = 3; ki < 6; ki++)
                      		node[i1].fluxResidual[ki] = node[i1].fluxResidual[ki] + fluxResidual[ki];
              		}
              		else if (tetrahedron[h].boundaryType == 'C')
              		{
                  		for (int ki = 0; ki < 6; ki++)
                      		node[i1].fluxResidual[ki] = node[i1].fluxResidual[ki] + fluxResidual[ki];
              		};

              		/////////////////////////////////////////////////////
              		// (j + 2) % 4
              		for (int ix = 0; ix < 3; ix++)
                  		coordinate[ix] = (node[i2].coordinate[ix]
                                              + (node[i2].coordinate[ix] + node[i3].coordinate[ix]) / 2.0
                                              + (node[i2].coordinate[ix] + node[i1].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

					for (int ix = 0; ix < 3; ix++)
						displacement[ix] = coordinate[ix] - node[tetrahedron[h].globalNode[j]].coordinate[ix];

					for (int ki = 0; ki < 6; ki++)
						firstDerivative[ki] = dotProduct(displacement, tetrahedron[h].gradient[ki]);

					for (int ki = 0; ki < 6; ki++)
	                  	interpolateVariable[ki] = node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
													+ firstDerivative[ki];

              		boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    		xOutwardNormal, yOutwardNormal, zOutwardNormal);

              		if (tetrahedron[h].boundaryType == 'B')
              		{
                  		for (int ki = 3; ki < 6; ki++)
                      		node[i2].fluxResidual[ki] = node[i2].fluxResidual[ki] + fluxResidual[ki];
              		}
              		else if (tetrahedron[h].boundaryType == 'C')
              		{
                  		for (int ki = 0; ki < 6; ki++)
                      		node[i2].fluxResidual[ki] = node[i2].fluxResidual[ki] + fluxResidual[ki];
              		};

              		/////////////////////////////////////////////////////
              		// (j + 3) % 4
              		for (int ix = 0; ix < 3; ix++)
                  		coordinate[ix] = (node[i3].coordinate[ix]
                                              + (node[i3].coordinate[ix] + node[i1].coordinate[ix]) / 2.0
                                              + (node[i3].coordinate[ix] + node[i2].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

					for (int ix = 0; ix < 3; ix++)
						displacement[ix] = coordinate[ix] - node[tetrahedron[h].globalNode[j]].coordinate[ix];

					for (int ki = 0; ki < 6; ki++)
						firstDerivative[ki] = dotProduct(displacement, tetrahedron[h].gradient[ki]);

					for (int ki = 0; ki < 6; ki++)
	                  	interpolateVariable[ki] = node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
													+ firstDerivative[ki];

              		boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    		xOutwardNormal, yOutwardNormal, zOutwardNormal);
			

              		if (tetrahedron[h].boundaryType == 'B')
              		{
                  		for (int ki = 3; ki < 6; ki++)
                      		node[i3].fluxResidual[ki] = node[i3].fluxResidual[ki] + fluxResidual[ki];
              		}
              		else if (tetrahedron[h].boundaryType == 'C')
              		{
                  		for (int ki = 0; ki < 6; ki++)
                      		node[i3].fluxResidual[ki] = node[i3].fluxResidual[ki] + fluxResidual[ki];
              		};
          		}
        	};
        	break;
    };
	delete [] fluxResidual; fluxResidual = NULL;
    delete [] coordinate; coordinate = NULL;
    delete [] centroid; centroid = NULL;
    delete [] interpolateVariable; interpolateVariable = NULL;
	delete [] displacement; displacement = NULL;
	delete [] firstDerivative; firstDerivative = NULL;
}

void Computation3D::boundaryCondition(const int & UVARIABLE_LEVEL)
{
	for (int i = 0; i < nodeNumber; i++)
	{
		if (node[i].boundary == 'A' || node[i].boundary == 'C' || node[i].boundary == 'D' || node[i].boundary == 'B')
		{
            timeDependentSolution(i, UVARIABLE_LEVEL, time);
		};
	};
}

void Computation3D::finiteVolumeNodalUpdate()
{
	/*/ (1/3) Simpson's Rule /**/
    double** k1 = new double* [nodeNumber];
    double** k2 = new double* [nodeNumber];
    double** k3 = new double* [nodeNumber];
    for (int i = 0; i < nodeNumber; i++) {
        k1[i] = new double [6];
        k2[i] = new double [6];
        k3[i] = new double [6];
    };
    
    /*/ Simpson's Rule Stage-1: (tn, yn) /**/
	calculateGradient(1);
	calculateFiniteVolume(1);
	calculateFiniteVolumeBoundaryFlux(1);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++) {
    		k1[i][ki] = node[i].fluxResidual[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].volume) * (k1[i][ki]);
    	};
	time -= (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-2: (t + tdelta / 2, y + tdelta * k1) /**/
	calculateGradient(2);
	calculateFiniteVolume(2);
	calculateFiniteVolumeBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++) {
    		k2[i][ki] = node[i].fluxResidual[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * (k2[i][ki]);
    	};
	time += (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-3: (t + tdelta, y + tdelta * k2) /**/
	calculateGradient(2);
	calculateFiniteVolume(2);
	calculateFiniteVolumeBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++) {
    		k3[i][ki] = node[i].fluxResidual[ki];
    	};
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++)
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * (k1[i][ki] + 4.0 * k2[i][ki] + k3[i][ki]) / 6.0;
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
    	for (int ki = 0; ki < 6; ki++)
    		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];


	// /*/ Runge-Kutta Stage-1 /**/
	// calculateGradient(1);
	// calculateFiniteVolume(1);
	// // calculateFiniteVolumeBoundaryFlux(1);
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 6; ki++)
	// 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].volume) * node[i].fluxResidual[ki];
	// time -= timeDelta / 2.0;
	// boundaryCondition(2);

	// /*/ Runge-Kutta Stage-2 /**/
	// calculateGradient(2);
	// calculateFiniteVolume(2);
	// // calculateFiniteVolumeBoundaryFlux(2);
	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 6; ki++)
	// 		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * node[i].fluxResidual[ki];
	// time += timeDelta / 2.0;
	// boundaryCondition(2);

	// for (int i = 0; i < nodeNumber; i++)
	// 	for (int ki = 0; ki < 6; ki++)
	// 		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation3D::fluxDifferenceNodalUpdate(const int & timeStep)
{
	/*/ (1/3) Simpson's Rule /**/
    double** k1 = new double* [nodeNumber];
    double** k2 = new double* [nodeNumber];
    double** k3 = new double* [nodeNumber];
    for (int i = 0; i < nodeNumber; i++) {
        k1[i] = new double [6];
        k2[i] = new double [6];
        k3[i] = new double [6];
    };
    
    /*/ Simpson's Rule Stage-1: (tn, yn) /**/
    calculateFluxDifference(1);
    calculateFluxDifferenceBoundaryFlux(1);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++) {
    		k1[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].volume) * (k1[i][ki]);
    	};
	time -= (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-2: (t + tdelta / 2, y + tdelta * k1) /**/
    calculateFluxDifference(2);
    calculateFluxDifferenceBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++) {
    		k2[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * (k2[i][ki]);
    	};
	time += (timeDelta / 2.0);
    boundaryCondition(2);
    
    /*/ Simpson's Rule Stage-3: (t + tdelta, y + tdelta * k2) /**/
    calculateFluxDifference(2);
    calculateFluxDifferenceBoundaryFlux(2);
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++) {
    		k3[i][ki] = node[i].fluxResidual[ki] - node[i].dissipation[ki];
    	};
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++)
    		node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * (k1[i][ki] + 4.0 * k2[i][ki] + k3[i][ki]) / 6.0;
    boundaryCondition(2);

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// compute the amplificationFactor
	for (int ki = 0; ki < 6; ki++)
		amplificationFactor[timeStep][ki] = 0.0;
	int count = 0;
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++)
		{
			float amplification_numerator = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * ((k1[i][ki] + 2.0 * k2[i][ki]) / 6.0);
			float amplification_denominator = node[i].conservedVariable[ki][2] + (timeDelta / node[i].volume) * ((2.0 * k2[i][ki] + k3[i][ki]) / 6.0);
        	if (amplification_denominator != 0.0)
			{
				amplificationFactor[timeStep][ki] = amplificationFactor[timeStep][ki] + amplification_numerator / amplification_denominator;
				count++;
			}
		};
	if (count != 0)
		for (int ki = 0; ki < 6; ki++)
			amplificationFactor[timeStep][ki] = amplificationFactor[timeStep][ki] / nodeNumber;
    
    for (int i = 0; i < nodeNumber; i++) {
        delete [] k1[i]; k1[i] = NULL;
        delete [] k2[i]; k2[i] = NULL;
        delete [] k3[i]; k3[i] = NULL;
    };
    delete [] k1; k1 = NULL;
    delete [] k2; k2 = NULL;
    delete [] k3; k3 = NULL;    
    
    for (int i = 0; i < nodeNumber; i++)
    	for (int ki = 0; ki < 6; ki++)
    		node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation3D::errorsCalculation()
{
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][3] = node[i].conservedVariable[ki][2] - node[i].conservedVariable[ki][7];
}

void Computation3D::interpolateCrossSection()
{
    if (randomization)
    {
        for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
        {
            double sum = 0.0;
            for (int j = 0; j < 4; j++)
            {
                int h = crossSection[crosssection].tetrahedron;
                int i = tetrahedron[h].globalNode[j];
                sum += LagrangeInterpolationFunction(crossSection[crosssection].coordinate, h, j) * node[i].conservedVariable[5][2];
            };
            crossSection[crosssection].Hz = sum;
        }
    }
    else
    {
        for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
        {
            int i = crossSection[crosssection].node;
            crossSection[crosssection].Hz = node[i].conservedVariable[5][2];
        };
    };
}

void Computation3D::printCrossSection() const
{
	// For cross-sectional plots in Asymptote
	///////////////////////////////////////////
	string stringTime;
	if (time == 0.5)
        stringTime = "t05";
	else if (time == 1.0)
        stringTime = "t1";
    else if (time == 1.5)
        stringTime = "t15";
    else if (time == 2.0)
        stringTime = "t2";


    ofstream results;
    if (method == 'A')
    {
        results.open("Spherical_Scattering_Cross_Section_" + stringTime + "_Exact.txt");
        results << crossSectionNumber << endl;
        for (int ix = 0; ix < crossSectionNumber; ix++)
        {
            results << setw(20) << crossSection[ix].coordinate[0]
                    << setw(20) << crossSection[ix].coordinate[2]
                    << setw(20) << crossSection[ix].Hz << endl;
        };
        results.close();
    };


    switch (method)
	{
        case 'A':
            results.open("Spherical_Scattering_Cross_Section_" + stringTime + "_SecondOrderFV.txt");
            break;
        case 'B':
            results.open("Spherical_Scattering_Cross_Section_" + stringTime + "_FluxDifference.txt");
            break;
        case 'C':
            results.open("Spherical_Scattering_Cross_Section_" + stringTime + "_RDGalerkin.txt");
            break;
        case 'D':
            results.open("Spherical_Scattering_Cross_Section_" + stringTime + "_RDLW.txt");
            break;
	};
    results << crossSectionNumber << endl;
    for (int ix = 0; ix < crossSectionNumber; ix++)
    {
        results << setw(20) << crossSection[ix].coordinate[0]
                << setw(20) << crossSection[ix].coordinate[2]
                << setw(20) << crossSection[ix].Hz << endl;
    };
    results.close();
}

void Computation3D::printResults() const
{
    ofstream results;

    // For cross-sectional plots in Asymptote
	///////////////////////////////////////////
	string stringTime;
	if (time == 0.5)
        stringTime = "t05";
	else if (time == 1.0)
        stringTime = "t1";
    else if (time == 1.5)
        stringTime = "t15";
    else if (time == 2.0)
        stringTime = "t2";



    if (method == 'A')
    {
        results.open("./Output/Spherical_Scattering_" + stringTime + "_Exact.txt");
        results << nodeNumber;
        for (int i = 0; i < nodeNumber; i++)
        {
            results << setw(20) << node[i].coordinate[0]
                    << setw(20) << node[i].coordinate[1]
                    << setw(20) << node[i].coordinate[2]
                    << setw(20) << node[i].conservedVariable[5][4] << endl;
        };
        results.close();
    };


    switch (method)
	{
        case 'A':
            results.open("Spherical_Scattering_" + stringTime + "_SecondOrderFV.txt");
            break;
        case 'B':
            results.open("Spherical_Scattering_" + stringTime + "_FluxDifference.txt");
            break;
        case 'C':
            results.open("Spherical_Scattering_" + stringTime + "_RDGalerkin.txt");
            break;
        case 'D':
            results.open("Spherical_Scattering_" + stringTime + "_RDLW.txt");
            break;
	};
    results << nodeNumber;
    for (int i = 0; i < nodeNumber; i++)
    {
        results << setw(20) << node[i].coordinate[0]
                << setw(20) << node[i].coordinate[1]
                << setw(20) << node[i].coordinate[2]
                << setw(20) << node[i].conservedVariable[5][2] << endl;
    };
    results.close();



    // to calculate L2-errors
    ///////////////////////////////////
	double averageVolume = 0.0;
	for (int h = 0; h < tetrahedronNumber; h++)
		averageVolume += tetrahedron[h].volume;
	averageVolume = averageVolume / tetrahedronNumber;

    if (abs(time - 0.5) <= timeDelta / 50.0 || abs(time - 1.0) <= timeDelta / 50.0
        || abs(time - 1.5) <= timeDelta / 50.0 || abs(time - 2.0) <= timeDelta / 50.0)
    {
        ofstream L2Errors;
        switch (method)
        {
            case 'A':
                L2Errors.open("L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_SecondOrderFV.txt");
                break;
            case 'B':
                L2Errors.open("L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_FluxDifference.txt");
                break;
            case 'C':
                L2Errors.open("L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_RDGalerkin.txt");
                break;
            case 'D':
                L2Errors.open("L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_RDLW.txt");
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
        L2Errors << setw(20) << log10(averageVolume) << setw(20) << log10(sum) << endl;

        L2Errors.close();
    };
}

void Computation3D::intervalResults(bool & TIME_05, bool & TIME_1, bool & TIME_15, bool & TIME_2) const
{
    if (time >= 0.5 && TIME_05)
    {
        printCrossSection();
        printResults();
        TIME_05 = false;
    }
    else if (time >= 1.0 && TIME_1)
    {
        printCrossSection();
        printResults();
        TIME_1 = false;
    }
    else if (time >= 1.5 && TIME_15)
    {
        printCrossSection();
        printResults();
        TIME_15 = false;
    }
    else if (time >= 2.0 && TIME_2)
    {
        printCrossSection();
        printResults();
        TIME_2 = false;
    }
}

void Computation3D::initializeArray()
{
	for (int h = 0; h < tetrahedronNumber; h++)
	{
		tetrahedron[h].medianDualNormal = new double* [6];
		for (int mediandual = 0; mediandual < 6; mediandual++)
		{
			tetrahedron[h].medianDualNormal[mediandual] = new double [3];
		};

		if (method == 'A')
		{
			if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
			{
				tetrahedron[h].inwardNormal = new double* [4];
				for (int j = 0; j < 4; j++)
				{
					tetrahedron[h].inwardNormal[j] = new double [3];
				};
			};
			tetrahedron[h].centroid = new float [3];
		}
		else if (method == 'B')
		{
			tetrahedron[h].inwardNormal = new double* [4];
			for (int j = 0; j < 4; j++)
			{
				tetrahedron[h].inwardNormal[j] = new double [3];
			};
		};

		if (method == 'B')
		{
			tetrahedron[h].centroid = new float [3];
			if (tetrahedron[h].boundaryType == 'B' || tetrahedron[h].boundaryType == 'C')
			{
				tetrahedron[h].gradient = new double* [6];
				for (int ki = 0; ki < 6; ki++)
					tetrahedron[h].gradient[ki] = new double [3];
			}
		};

		tetrahedron[h].subVolume = new double [4];
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		node[i].fluxResidual = new double [6];
		node[i].dissipation = new double [6];

		node[i].conservedVariable = new double* [6];
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki] = new double [8];

		if (method == 'A')
		{
			node[i].firstDerivative = new double* [6];
			for (int ki = 0; ki < 6; ki++)
				node[i].firstDerivative[ki] = new double [3];
		};
	};


	crossSection = new crossSectionArray [crossSectionNumber];
	for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
	{
		crossSection[crosssection].coordinate = new double [3];
	};
}

void Computation3D::timeComputations(const char & TM_TE_Mode, const char & Method, const double & Time_Last, const double & Time_Delta)
{
	TMTEMode = TM_TE_Mode;
	method = Method;
	initializeArray();

	constructMedianDualNormal();
	calculateSubVolume();
	medianCellVolume();
	findTetrahedronOrientation();

    if (method == 'A')
    {
		constructNodeNeighbour();
        calculateInwardNormalFromMedianDualNormal();
		calculateCentroid();
    }
	else if (method == 'B')
	{
		constructNodeNeighbour();
		calculateInwardNormal();
		calculateCentroid();
	}

	interpolateCrossSection();

    /*/
    if (randomization)
        findCrossSection();
    /**/

    if (randomization)
    {
        int number;
        ifstream inputCrossSection;
        inputCrossSection.open("./Output/Cross Section (Randomized).txt");
        inputCrossSection >> number;
        for (int ix = 0; ix < crossSectionNumber; ix++)
            inputCrossSection >> crossSection[ix].coordinate[0]
                    			>> crossSection[ix].coordinate[1]
                    			>> crossSection[ix].coordinate[2]
                    			>> crossSection[ix].node
                    			>> crossSection[ix].tetrahedron;
        inputCrossSection.close();
    };


	timeLast = Time_Last;

	if (Time_Delta == 0.0)
	{
		timeDelta = globalTimeStep();
	}
	else
	{
		timeDelta = Time_Delta;
	}

	// to correlate the timeNumber exactly halt at timeLast / 4, timeLast / 2, 0.75 * timeLast, timeLast
	// timeNumber = static_cast <int> ((timeLast / 4.0) / timeDelta) * 4;
	// timeDelta = timeLast / timeNumber;

	timeNumber = static_cast <int> (timeLast / timeDelta);
	time = 0.0;

	// construct the **amplificationFactor pointer
	amplificationFactor = new float* [timeNumber];
	for (int t = 0; t < timeNumber; t++)
		amplificationFactor[t] = new float [6];

    spatialSolution();
    fieldInitialization();

    /**/
    bool time05 = true;
    bool time1 = true;
    bool time15 = true;
    bool time2 = true;

    ofstream outputTime;
	clock_t STARTTIME;

	switch (method)
	{
        case 'A':
            outputTime.open("Time Record (Second-Order FV) " + to_string(tetrahedronNumber) + ".txt");
            break;
        case 'B':
            outputTime.open("Time Record (Flux-Difference) " + to_string(tetrahedronNumber) + ".txt");
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
		L2Errors.open("L2-Errors-3D-RK2.txt");
		for (int i = 0; i < nodeNumber; i++)
			sum += pow(node[i].conservedVariable[2][7] - node[i].conservedVariable[2][2], 2);
		L2Errors << showpoint << setprecision(10) << log10(sqrt(sum / nodeNumber)) << endl;
		L2Errors.close();
    }
    else if (method == 'B')
    {
		for (int t = 1; t <= timeNumber; t++)
        {
            time = t * timeDelta;
            cout << "Time is " << time << endl;
            fluxDifferenceNodalUpdate(t - 1);
            intervalResults(time05, time1, time15, time2);
        };
		double sum = 0.0;
		ofstream L2Errors;
		L2Errors.open("L2-Errors-3D-RK2.txt");
		for (int i = 0; i < nodeNumber; i++)
			sum += pow(node[i].conservedVariable[2][7] - node[i].conservedVariable[2][2], 2);
		L2Errors << showpoint << setprecision(10) << log10(sqrt(sum / nodeNumber)) << endl;
		L2Errors.close();
    };
	// output the amplificationFactor to .csv file
	ofstream amplificationFactorFile;
	amplificationFactorFile.open("AmplificationFactor_3D.csv");
	amplificationFactorFile << "time,U1,U2,U3,U4,U5,U6\n";
	for (int t = 0; t < timeNumber; t++)
	{
		amplificationFactorFile << (t + 1) * timeDelta;
		for (int ki = 0; ki < 6; ki++)
			amplificationFactorFile << "," << amplificationFactor[t][ki];
		amplificationFactorFile << "\n";		
	};
	amplificationFactorFile.close();
	// destruct the **amplificationFactor pointer
	for (int t = 0; t < timeNumber; t++)
    {
        delete [] amplificationFactor[t]; amplificationFactor[t] = NULL;
    };
    delete [] amplificationFactor; amplificationFactor = NULL;

	ofstream results;
	// print 3D results
	switch (method)
	{
		case 'A':
			results.open("Second_Order_Finite_Volume_3D.vtk");
			break;
		case 'B':
			results.open("Flux_Difference_3D.vtk");
			break;
	};
	results << "# vtk DataFile Version 2.0" << endl;
	results << "3D Sample" << endl;
	results << "ASCII" << endl << endl;

	results << "DATASET UNSTRUCTURED_GRID" << endl;
	results << "POINTS " << nodeNumber << " float" << endl;
	for (int i = 0; i < nodeNumber; i++)
		results << showpoint << setprecision(10)
							<< setw(20) << node[i].coordinate[0]
							<< setw(20) << node[i].coordinate[1]
							<< setw(20) << node[i].coordinate[2] << endl;
	results << endl;

	results << "CELLS " << tetrahedronNumber << setw(10) << tetrahedronNumber * (4 + 1) << endl;
	for (int h = 0; h < tetrahedronNumber; h++)
		results << setw(12) << "4"
							<< setw(10) << tetrahedron[h].globalNode[0]
							<< setw(10) << tetrahedron[h].globalNode[1]
							<< setw(10) << tetrahedron[h].globalNode[2]
							<< setw(10) << tetrahedron[h].globalNode[3] << endl;
	results << endl;

	results << "CELL_TYPES " << tetrahedronNumber << endl;
	for (int h = 0; h < tetrahedronNumber; h++)
		results << setw(12) << "10" << endl;
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

	// print cross-section y = -0.9
	switch (method)
	{
        case 'A':
            results.open("Second_Order_Finite_Volume_3D_iy=1.vtk");
            break;
        case 'B':
            results.open("Flux_Difference_3D_iy=1.vtk");
            break;
	};
    results << "# vtk DataFile Version 2.0" << endl;
    results << "Spherical Scattering" << endl;
    results << "ASCII" << endl << endl;

	int xNumber = 20;
	int zNumber = 20;
    results << "DATASET UNSTRUCTURED_GRID" << endl;
    results << "POINTS " << (xNumber + 1) * (zNumber + 1) << " float" << endl;
	int crossSectionNumber = 0;
    for (int i = 0; i < nodeNumber; i++)
		if (abs(node[i].coordinate[1] - (- 0.9)) < 0.0001)
        {
			crossSectionNumber++;
			results << showpoint << setprecision(10)
                           		<< setw(20) << node[i].coordinate[0]
                           		<< setw(20) << node[i].coordinate[2]
                           		<< setw(20) << 0.0 << endl;
		};
	results << endl;

	int elementNumber = xNumber * zNumber;
    results << "CELLS " << elementNumber << setw(10) << elementNumber * (4 + 1) << endl;
    for (int iz = 0; iz < zNumber; iz++)
		for (int ix = 0; ix < xNumber; ix++)
	    	results << setw(12) << "4"
                    << setw(10) << iz * (xNumber + 1) + ix
                    << setw(10) << (iz + 1) * (xNumber + 1) + ix
					<< setw(10) << (iz + 1) * (xNumber + 1) + (ix + 1)
                    << setw(10) << iz * (xNumber + 1) + (ix + 1) << endl;
    results << endl;

    results << "CELL_TYPES " << elementNumber << endl;
    for (int e = 0; e < elementNumber; e++)
        results << setw(12) << "9" << endl;
    results << endl;

    results << "POINT_DATA " << crossSectionNumber << endl;
    results << "SCALARS Ez_Numerical float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < crossSectionNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][2] << endl;
    results << endl;

    results << "SCALARS Ez_Analytical_Q1 float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < crossSectionNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][4] << endl;
    results << endl;

	results << "SCALARS Ez_Analytical_Q2 float 1" << endl;
	results << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < crossSectionNumber; i++)
		results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][5] << endl;
	results << endl;

	results << "SCALARS Ez_Analytical_Q3 float 1" << endl;
	results << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < crossSectionNumber; i++)
		results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][6] << endl;
	results << endl;

	results << "SCALARS Ez_Analytical_Q4 float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < crossSectionNumber; i++)
        results << showpoint << setprecision(10) << setw(20) << node[i].conservedVariable[2][7] << endl;
    results << endl;
	results.close();

	outputTime << "The global time step is " << timeDelta << endl;
	outputTime << "The execution time is " << static_cast <double> ((clock() - STARTTIME) / static_cast <double> (CLOCKS_PER_SEC)) << ". " << endl;
	outputTime.close();
	errorsCalculation();
  	/**/
}

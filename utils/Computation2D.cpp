# include "Computation2D.h"

Computation2D::~Computation2D()
{
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

            for (int j = 0; j < 3; j++)
            {
                delete [] element[e].interpolateVariable[j]; element[e].interpolateVariable[j] = NULL;
            };
            delete [] element[e].interpolateVariable; element[e].interpolateVariable = NULL;
        }
        else if (method == 'B' || method == 'C' || method == 'D')
        {
            for (int j = 0; j < 3; j++)
            {
                delete [] element[e].inwardNormal[j]; element[e].inwardNormal[j] = NULL;
            };
            delete [] element[e].inwardNormal; element[e].inwardNormal = NULL;
        };

        if (method == 'D')
        {
            for (int j = 0; j < 3; j++)
            {
                for (int ki = 0; ki < 3; ki++)
                {
                    delete [] element[e].distributionMatrix[j][ki]; element[e].distributionMatrix[j][ki] = NULL;
                };
                delete [] element[e].distributionMatrix[j]; element[e].distributionMatrix[j] = NULL;
            };
            delete [] element[e].distributionMatrix; element[e].distributionMatrix = NULL;
        };
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		for (int ki = 0; ki < 3; ki++)
		{
			delete [] node[i].conservedVariable[ki]; node[i].conservedVariable[ki] = NULL;
		};
		delete [] node[i].conservedVariable; node[i].conservedVariable = NULL;
		delete [] node[i].fluxResidual; node[i].fluxResidual = NULL;
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

void Computation2D::spatialSolution()
{
		for (int i = 0; i < nodeNumber; i++)
		{
			complexNumber Hr, Hphi, Ez;

			Hr.real = 0.0;
			Hr.imaginary = 0.0;
			Hphi.real = (- propagationCoefficient / (angularFrequency * permeability)) * cos(- propagationCoefficient * node[i].coordinate[0]);
			Hphi.imaginary = (- propagationCoefficient / (angularFrequency * permeability)) * sin(- propagationCoefficient * node[i].coordinate[0]);
			Ez.real = cos(- propagationCoefficient * node[i].coordinate[0]);
			Ez.imaginary = sin(- propagationCoefficient * node[i].coordinate[0]);

			node[i].Hr = Hr;
			node[i].Hphi = Hphi;
			node[i].Ez = Ez;
		}
}

void Computation2D::timeDependentSolution(const int & I, const int & UVARIABLE_LEVEL, const double & Time)
{
		complexNumber timeHarmonic;
		complexNumber HR, HPhi, EZ;

		timeHarmonic.real = cos(angularFrequency * Time);
		timeHarmonic.imaginary = sin(angularFrequency * Time);

		HR = complexComplexMultiplication(timeHarmonic, node[I].Hr);
		HPhi = complexComplexMultiplication(timeHarmonic, node[I].Hphi);
		EZ = complexComplexMultiplication(timeHarmonic, node[I].Ez);

		node[I].conservedVariable[0][UVARIABLE_LEVEL] = HR.real;
		node[I].conservedVariable[1][UVARIABLE_LEVEL] = HPhi.real;
		node[I].conservedVariable[2][UVARIABLE_LEVEL] = EZ.real;
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

void Computation2D::interpolateMedianDualCenter(const int & UVARIABLE_LEVEL) const
{
	for (int e = 0; e < elementNumber; e++)
	{
		for (int j = 0; j < 3; j++)
		{
			double* vector;
			vector = new double [2];

			int vertex;
			double xCoordinate;
			double yCoordinate;

			/*/	side 0 /**/
			vector[0] = element[e].edgeTangent[j][0] / 2.0 - element[e].medianDualNormal[j][1] / 2.0;
			vector[1] = element[e].edgeTangent[j][1] / 2.0 + element[e].medianDualNormal[j][0] / 2.0;
			switch (j)
			{
				case 0:
					vertex = element[e].globalNode[1];
					xCoordinate = node[element[e].globalNode[1]].coordinate[0] + vector[0];
					yCoordinate = node[element[e].globalNode[1]].coordinate[1] + vector[1];
					break;
				case 1:
					vertex = element[e].globalNode[2];
					xCoordinate = node[element[e].globalNode[2]].coordinate[0] + vector[0];
					yCoordinate = node[element[e].globalNode[2]].coordinate[1] + vector[1];
					break;
				case 2:
					vertex = element[e].globalNode[0];
					xCoordinate = node[element[e].globalNode[0]].coordinate[0] + vector[0];
					yCoordinate = node[element[e].globalNode[0]].coordinate[1] + vector[1];
					break;
			};

			/*/	side 1 /**/
			/*/
			vector[0] = - element[e].edgeTangent[j][0] / 2.0 - element[e].medianDualNormal[j][1] / 2.0;
			vector[1] = - element[e].edgeTangent[j][1] / 2.0 + element[e].medianDualNormal[j][0] / 2.0;
			switch (j)
			{
				case 0:
					vertex = element[e].globalNode[1];
					xCoordinate = node[element[e].globalNode[2]].coordinate[0] + vector[0];
					yCoordinate = node[element[e].globalNode[2]].coordinate[1] + vector[1];
					break;
				case 1:
					vertex = element[e].globalNode[2];
					xCoordinate = node[element[e].globalNode[0]].coordinate[0] + vector[0];
					yCoordinate = node[element[e].globalNode[0]].coordinate[1] + vector[1];
					break;
				case 2:
					vertex = element[e].globalNode[0];
					xCoordinate = node[element[e].globalNode[1]].coordinate[0] + vector[0];
					yCoordinate = node[element[e].globalNode[1]].coordinate[1] + vector[1];
					break;
			};
			/**/

			for (int ki = 0; ki < 3; ki++)
			{
				element[e].interpolateVariable[j][ki] = 0.0;
				for (int vertex = 0; vertex < 3; vertex++)
					element[e].interpolateVariable[j][ki] = element[e].interpolateVariable[j][ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
			};

			delete [] vector; vector = NULL;
		}
	};
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

	switch (TMTEMode)
	{
        case 'A':
            for (int e = 0; e < elementNumber; e++)
            {
                for (int j = 0; j < 3; j++)
                {
                    int i = element[e].globalNode[j];
                    double sum = 0.0;
                    for (int side = 0; side < 2; side++)
                    {
                        double xMedianDualNormal;
                        double yMedianDualNormal;
                        double xMagneticField;
                        double yMagneticField;
                        double zElectricField;
                        switch (side)
                        {
                            case 0:
                                xMedianDualNormal = element[e].medianDualNormal[(j + 2) % 3][0];
                                yMedianDualNormal = element[e].medianDualNormal[(j + 2) % 3][1];
                                xMagneticField = element[e].interpolateVariable[(j + 2) % 3][0];
                                yMagneticField = element[e].interpolateVariable[(j + 2) % 3][1];
                                zElectricField = element[e].interpolateVariable[(j + 2) % 3][2];
                                break;
                            case 1:
                                xMedianDualNormal = - element[e].medianDualNormal[(j + 1) % 3][0];
                                yMedianDualNormal = - element[e].medianDualNormal[(j + 1) % 3][1];
                                xMagneticField = element[e].interpolateVariable[(j + 1) % 3][0];
                                yMagneticField = element[e].interpolateVariable[(j + 1) % 3][1];
                                zElectricField = element[e].interpolateVariable[(j + 1) % 3][2];
                                break;
                        };
                        node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- 0.0 * xMedianDualNormal + (yMedianDualNormal / permeability) * zElectricField);
                        node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- (xMedianDualNormal / permeability) * zElectricField + 0.0 * yMedianDualNormal);
                        node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- (xMedianDualNormal / permittivity) * yMagneticField
                                                                                + (yMedianDualNormal / permittivity) * xMagneticField);
                    }
                }
            };
            break;
        case 'B':
            for (int e = 0; e < elementNumber; e++)
            {
                for (int j = 0; j < 3; j++)
                {
                    int i = element[e].globalNode[j];
                    double sum = 0.0;
                    for (int side = 0; side < 2; side++)
                    {
                        double xMedianDualNormal;
                        double yMedianDualNormal;
                        double xElectricField;
                        double yElectricField;
                        double zMagneticField;
                        switch (side)
                        {
                            case 0:
                                xMedianDualNormal = element[e].medianDualNormal[(j + 2) % 3][0];
                                yMedianDualNormal = element[e].medianDualNormal[(j + 2) % 3][1];
                                xElectricField = element[e].interpolateVariable[(j + 2) % 3][0];
                                yElectricField = element[e].interpolateVariable[(j + 2) % 3][1];
                                zMagneticField = element[e].interpolateVariable[(j + 2) % 3][2];
                                break;
                            case 1:
                                xMedianDualNormal = - element[e].medianDualNormal[(j + 1) % 3][0];
                                yMedianDualNormal = - element[e].medianDualNormal[(j + 1) % 3][1];
                                xElectricField = element[e].interpolateVariable[(j + 1) % 3][0];
                                yElectricField = element[e].interpolateVariable[(j + 1) % 3][1];
                                zMagneticField = element[e].interpolateVariable[(j + 1) % 3][2];
                                break;
                        };
                        node[i].fluxResidual[0] = node[i].fluxResidual[0] + (0.0 * xMedianDualNormal - (yMedianDualNormal / permittivity) * zMagneticField);
                        node[i].fluxResidual[1] = node[i].fluxResidual[1] + ((xMedianDualNormal / permittivity) * zMagneticField - 0.0 * yMedianDualNormal);
                        node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- (yMedianDualNormal / permeability) * xElectricField
                                                                                + (xMedianDualNormal / permeability) * yElectricField);
                    }
                }
            };
            break;
	}
}

void Computation2D::TMFluxDifference(const int & E, const int & UVARIABLE_LEVEL)
{
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];
				double xMagneticField = 0.0;
				double yMagneticField = 0.0;
				double zElectricField = 0.0;

				for (int vertex = 0; vertex < 3; vertex++)
				{
						xMagneticField += node[element[E].globalNode[vertex]].conservedVariable[0][UVARIABLE_LEVEL];
						yMagneticField += node[element[E].globalNode[vertex]].conservedVariable[1][UVARIABLE_LEVEL];
						zElectricField += node[element[E].globalNode[vertex]].conservedVariable[2][UVARIABLE_LEVEL];
				};

				xMagneticField = xMagneticField / 3.0;
				yMagneticField = yMagneticField / 3.0;
				zElectricField = zElectricField / 3.0;

				double xMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][0];
				double yMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][1];


				node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- 0.0 * xMedianDualNormal + (zElectricField / permeability) * yMedianDualNormal);
				node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- (zElectricField / permeability) * xMedianDualNormal + 0.0 * yMedianDualNormal);
				node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- (yMagneticField / permittivity) * xMedianDualNormal + (xMagneticField / permittivity) * yMedianDualNormal);
		};
}

void Computation2D::TEFluxDifference(const int & E, const int & UVARIABLE_LEVEL)
{
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];
				double xElectricField = 0.0;
				double yElectricField = 0.0;
				double zMagneticField = 0.0;

				for (int vertex = 0; vertex < 3; vertex++)
				{
						xElectricField += node[element[E].globalNode[vertex]].conservedVariable[0][UVARIABLE_LEVEL];
						yElectricField += node[element[E].globalNode[vertex]].conservedVariable[1][UVARIABLE_LEVEL];
						zMagneticField += node[element[E].globalNode[vertex]].conservedVariable[2][UVARIABLE_LEVEL];
				};

				xElectricField = xElectricField / 3.0;
				yElectricField = yElectricField / 3.0;
				zMagneticField = zMagneticField / 3.0;

				double xMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][0];
				double yMedianDualNormal = - (1.0 / 2.0) * element[E].inwardNormal[j][1];


				node[i].fluxResidual[0] = node[i].fluxResidual[0] + (0.0 * xMedianDualNormal - (zMagneticField / permittivity) * yMedianDualNormal);
				node[i].fluxResidual[1] = node[i].fluxResidual[1] + ((zMagneticField / permittivity) * xMedianDualNormal - 0.0 * yMedianDualNormal);
				node[i].fluxResidual[2] = node[i].fluxResidual[2] + ((yElectricField / permeability) * xMedianDualNormal - (xElectricField / permeability) * yMedianDualNormal);
		};
}

void Computation2D::calculateFluxDifference(const int & UVARIABLE_LEVEL)
{
	for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 3; ki++)
            node[i].fluxResidual[ki] = 0.0;

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

void Computation2D::TMLocalFluxResidual(const int & E, const int & UVARIABLE_LEVEL, double* const & fluxResidual)
{
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];
				double xMagneticField = node[i].conservedVariable[0][UVARIABLE_LEVEL];
				double yMagneticField = node[i].conservedVariable[1][UVARIABLE_LEVEL];
				double zElectricField = node[i].conservedVariable[2][UVARIABLE_LEVEL];
				double xInwardNormal = element[E].inwardNormal[j][0];
				double yInwardNormal = element[E].inwardNormal[j][1];

				if (element[E].boundaryType == 'A'
						|| element[E].boundaryType == 'C'
						|| element[E].boundaryType == 'D'
						|| element[E].boundaryType == 'N')
				{
						fluxResidual[0] = fluxResidual[0] + (1.0 / 2.0) * (- 0.0 * xInwardNormal + (zElectricField / permeability) * yInwardNormal);
						fluxResidual[1] = fluxResidual[1] + (1.0 / 2.0) * (- (zElectricField / permeability) * xInwardNormal + 0.0 * yInwardNormal);
						fluxResidual[2] = fluxResidual[2] + (1.0 / 2.0) * (- (yMagneticField / permittivity) * xInwardNormal + (xMagneticField / permittivity) * yInwardNormal);
				};
		}
}

void Computation2D::TELocalFluxResidual(const int & E, const int & UVARIABLE_LEVEL, double* const & fluxResidual)
{
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];
				double xElectricField = node[i].conservedVariable[0][UVARIABLE_LEVEL];
				double yElectricField = node[i].conservedVariable[1][UVARIABLE_LEVEL];
				double zMagneticField = node[i].conservedVariable[2][UVARIABLE_LEVEL];
				double xInwardNormal = element[E].inwardNormal[j][0];
				double yInwardNormal = element[E].inwardNormal[j][1];

				if (element[E].boundaryType == 'A'
						|| element[E].boundaryType == 'C'
						|| element[E].boundaryType == 'D'
						|| element[E].boundaryType == 'N')
				{
						fluxResidual[0] = fluxResidual[0] + (1.0 / 2.0) * (0.0 * xInwardNormal - (zMagneticField / permittivity) * yInwardNormal);
						fluxResidual[1] = fluxResidual[1] + (1.0 / 2.0) * ((zMagneticField / permittivity) * xInwardNormal - 0.0 * yInwardNormal);
						fluxResidual[2] = fluxResidual[2] + (1.0 / 2.0) * ((yElectricField / permeability) * xInwardNormal - (xElectricField / permeability) * yInwardNormal);
				};
		};
}

void Computation2D::TMRDGalerkin(const int & E, const int & UVARIABLE_LEVEL)
{
		double* fluxResidual;
		fluxResidual = new double [3];

		for (int ki = 0; ki < 3; ki++)
				fluxResidual[ki] = 0.0;

		TMLocalFluxResidual(E, UVARIABLE_LEVEL, fluxResidual);

		// Distributed Residual
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];

				for (int ki = 0; ki < 3; ki++)
				{
						node[i].fluxResidual[ki] = node[i].fluxResidual[ki] + (1.0 / 3.0) * fluxResidual[ki];
				};
		};

		delete [] fluxResidual; fluxResidual = NULL;
}

void Computation2D::TERDGalerkin(const int & E, const int & UVARIABLE_LEVEL)
{
		double* fluxResidual;
		fluxResidual = new double [3];

		for (int ki = 0; ki < 3; ki++)
				fluxResidual[ki] = 0.0;

		TELocalFluxResidual(E, UVARIABLE_LEVEL, fluxResidual);

		// Distributed Residual
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];

				for (int ki = 0; ki < 3; ki++)
				{
						node[i].fluxResidual[ki] = node[i].fluxResidual[ki] + (1.0 / 3.0) * fluxResidual[ki];
				};
		};

		delete [] fluxResidual; fluxResidual = NULL;
}

void Computation2D::calculateRDGalerkin(const int & UVARIABLE_LEVEL)
{
    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 3; ki++)
            node[i].fluxResidual[ki] = 0.0;

        switch (TMTEMode)
        {
						case 'A':
								for (int e = 0; e < elementNumber; e++)
								{
										if (element[e].boundaryType == 'A'
												|| element[e].boundaryType == 'C'
												|| element[e].boundaryType == 'D'
												|| element[e].boundaryType == 'N')
										{
												TMRDGalerkin(e, UVARIABLE_LEVEL);
										}
										else if (element[e].boundaryType == 'B')
										{
												// TMRDGALERKIN(e, UVARIABLE_LEVEL);
												TMFluxDifference(e, UVARIABLE_LEVEL);
										};
								};
								break;
						case 'B':
								for (int e = 0; e < elementNumber; e++)
								{
										if (element[e].boundaryType == 'A'
												|| element[e].boundaryType == 'C'
												|| element[e].boundaryType == 'D'
												|| element[e].boundaryType == 'N')
										{
												TERDGalerkin(e, UVARIABLE_LEVEL);
										}
										else if (element[e].boundaryType == 'B')
										{
												// TERDGALERKIN(e, UVARIABLE_LEVEL);
												TEFluxDifference(e, UVARIABLE_LEVEL);
										};
								};
                break;
        };
}

void Computation2D::TMLaxWendroff(const int & E, const int & UVARIABLE_LEVEL)
{
		double* fluxResidual;
		fluxResidual = new double [3];

		for (int ki = 0; ki < 3; ki++)
				fluxResidual[ki] = 0.0;

		TMLocalFluxResidual(E, UVARIABLE_LEVEL, fluxResidual);

		// Distributed Residual
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];

				for (int ki = 0; ki < 3; ki++)
				{
						double sum = 0.0;
						for (int kj = 0; kj < 3; kj++)
								sum += element[E].distributionMatrix[j][ki][kj] * fluxResidual[kj];
						node[i].fluxResidual[ki] = node[i].fluxResidual[ki] + sum;
				};
		};

		delete [] fluxResidual; fluxResidual = NULL;
}

void Computation2D::TELaxWendroff(const int & E, const int & UVARIABLE_LEVEL)
{
		double* fluxResidual;
		fluxResidual = new double [3];

		for (int ki = 0; ki < 3; ki++)
				fluxResidual[ki] = 0.0;

		TELocalFluxResidual(E, UVARIABLE_LEVEL, fluxResidual);

		// Distributed Residual
		for (int j = 0; j < 3; j++)
		{
				int i = element[E].globalNode[j];

				for (int ki = 0; ki < 3; ki++)
				{
						double sum = 0.0;
						for (int kj = 0; kj < 3; kj++)
								sum += element[E].distributionMatrix[j][ki][kj] * fluxResidual[kj];
						node[i].fluxResidual[ki] = node[i].fluxResidual[ki] + sum;
				};
		};

		delete [] fluxResidual; fluxResidual = NULL;
}

void Computation2D::calculateLaxWendroff(const int & UVARIABLE_LEVEL)
{
    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 3; ki++)
            node[i].fluxResidual[ki] = 0.0;


		switch (TMTEMode)
		{
				case 'A':
						for (int e = 0; e < elementNumber; e++)
						{
								if (element[e].boundaryType == 'A'
										|| element[e].boundaryType == 'C'
										|| element[e].boundaryType == 'D'
										|| element[e].boundaryType == 'N')
								{
										TMLaxWendroff(e, UVARIABLE_LEVEL);
								}
								else if (element[e].boundaryType == 'B')
								{
										// TMLAXWENDROFF(e, UVARIABLE_LEVEL);
										TMFluxDifference(e, UVARIABLE_LEVEL);
								};
						};
						break;
				case 'B':
						for (int e = 0; e < elementNumber; e++)
						{
								if (element[e].boundaryType == 'A'
										|| element[e].boundaryType == 'C'
										|| element[e].boundaryType == 'D'
										|| element[e].boundaryType == 'N')
								{
										TELaxWendroff(e, UVARIABLE_LEVEL);
								}
								else if (element[e].boundaryType == 'B')
								{
										// TELAXWENDROFF(e, UVARIABLE_LEVEL);
										TEFluxDifference(e, UVARIABLE_LEVEL);
								};
						};
						break;
		};
}

void Computation2D::calculateFiniteVolumeBoundaryFlux(const int & UVARIABLE_LEVEL)
{
	switch (TMTEMode)
	{
		case 'A':
		    for (int e = 0; e < elementNumber; e++)
            {
				if (element[e].boundaryType == 'B')
				{
										int j = element[e].boundaryVertex;

										double* interpolateVariable;
										interpolateVariable = new double [3];

										/*********************************************/
										// Boundary on vertex (j + 2) % 3
										int i = element[e].globalNode[(j + 2) % 3];
										double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										// node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
										// node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];


										/*********************************************/
										// Boundary on vertex (j + 1) % 3
										i = element[e].globalNode[(j + 1) % 3];
										xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};


										xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


										// node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
										// node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];


										delete [] interpolateVariable; interpolateVariable = NULL;
								}
				else if (element[e].boundaryType == 'C')
								{
										int j = element[e].boundaryVertex;

										double* interpolateVariable;
										interpolateVariable = new double [3];

										/*********************************************/
										// Boundary on vertex (j + 2) % 3
										int i = element[e].globalNode[(j + 2) % 3];
										double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
										node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];


										/*********************************************/
										// Boundary on vertex (j + 1) % 3
										i = element[e].globalNode[(j + 1) % 3];
										xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};


										xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


										node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
										node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];


										delete [] interpolateVariable; interpolateVariable = NULL;
								}
            };
		    break;
        case 'B':
            for (int e = 0; e < elementNumber; e++)
            {
								if (element[e].boundaryType == 'B')
								{
										int j = element[e].boundaryVertex;

										double* interpolateVariable;
										interpolateVariable = new double [3];

										/*********************************************/
										// Boundary on vertex (j + 2) % 3
										int i = element[e].globalNode[(j + 2) % 3];
										double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										//node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										//node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];


										/*********************************************/
										// Boundary on vertex (j + 1) % 3
										i = element[e].globalNode[(j + 1) % 3];
										xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


										//node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										//node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

										delete [] interpolateVariable; interpolateVariable = NULL;
								}
								else if (element[e].boundaryType == 'C')
								{
										int j = element[e].boundaryVertex;

										double* interpolateVariable;
										interpolateVariable = new double [3];

										/*********************************************/
										// Boundary on vertex (j + 2) % 3
										int i = element[e].globalNode[(j + 2) % 3];
										double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];


										/*********************************************/
										// Boundary on vertex (j + 1) % 3
										i = element[e].globalNode[(j + 1) % 3];
										xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


										node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

										delete [] interpolateVariable; interpolateVariable = NULL;
								};
            };
            break;
	};
}

void Computation2D::calculateFluxDifferenceBoundaryFlux(const int & UVARIABLE_LEVEL)
{
	switch (TMTEMode)
	{
			case 'A':
		    	for (int e = 0; e < elementNumber; e++)
          {
			  if (element[e].boundaryType == 'B')
							{
									int j = element[e].boundaryVertex;

									double* interpolateVariable;
									interpolateVariable = new double [3];

									/*********************************************/
									// Boundary on vertex (j + 2) % 3
									int i = element[e].globalNode[(j + 2) % 3];

									// for (int ki = 0; ki < 3; ki++)
									// {
									// 		interpolateVariable[ki] = 0.0;
									// 		for (int vertex = 0; vertex < 3; vertex++)
									// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
									// };

									double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
									double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

									for (int ki = 0; ki < 3; ki++)
									{
											interpolateVariable[ki] = 0.0;
											for (int vertex = 0; vertex < 3; vertex++)
													interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									};

									double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
									double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

									// node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
									// node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
									node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

									/*********************************************/
									// Boundary on vertex (j + 1) % 3
									i = element[e].globalNode[(j + 1) % 3];

									// for (int ki = 0; ki < 3; ki++)
									// {
									// 		interpolateVariable[ki] = 0.0;
									// 		for (int vertex = 0; vertex < 3; vertex++)
									// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
									// };

									xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
									yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

									for (int ki = 0; ki < 3; ki++)
									{
											interpolateVariable[ki] = 0.0;
											for (int vertex = 0; vertex < 3; vertex++)
													interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									};


									xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
									yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


									// node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
									// node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
									node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

									delete [] interpolateVariable; interpolateVariable = NULL;
							}
			  else if (element[e].boundaryType == 'C')
							{
									int j = element[e].boundaryVertex;

									double* interpolateVariable;
									interpolateVariable = new double [3];

									/*********************************************/
									// Boundary on vertex (j + 2) % 3
									int i = element[e].globalNode[(j + 2) % 3];

									// for (int ki = 0; ki < 3; ki++)
									// {
									// 		interpolateVariable[ki] = 0.0;
									// 		for (int vertex = 0; vertex < 3; vertex++)
									// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
									// };

									double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
									double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

									for (int ki = 0; ki < 3; ki++)
									{
											interpolateVariable[ki] = 0.0;
											for (int vertex = 0; vertex < 3; vertex++)
													interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									};

									double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
									double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

									node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
									node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
									node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

									/*********************************************/
									// Boundary on vertex 1
									i = element[e].globalNode[(j + 1) % 3];

									// for (int ki = 0; ki < 3; ki++)
									// {
									// 		interpolateVariable[ki] = 0.0;
									// 		for (int vertex = 0; vertex < 3; vertex++)
									// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
									// };

									xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
									yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

									for (int ki = 0; ki < 3; ki++)
									{
											interpolateVariable[ki] = 0.0;
											for (int vertex = 0; vertex < 3; vertex++)
													interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
									};


									xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
									yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


									node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
									node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
									node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

									delete [] interpolateVariable; interpolateVariable = NULL;
							}
          };
		    	break;
        case 'B':
            for (int e = 0; e < elementNumber; e++)
            {
								if (element[e].boundaryType == 'B')
								{
										int j = element[e].boundaryVertex;

										double* interpolateVariable;
										interpolateVariable = new double [3];

										/*********************************************/
										// Boundary on vertex (j + 2) % 3
										int i = element[e].globalNode[(j + 2) % 3];

										// for (int ki = 0; ki < 3; ki++)
										// {
										// 		interpolateVariable[ki] = 0.0;
										// 		for (int vertex = 0; vertex < 3; vertex++)
										// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										// 				interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
										// };

										double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										//node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										//node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

										/*********************************************/
										// Boundary on vertex (j + 1) % 3
										i = element[e].globalNode[(j + 1) % 3];

										// for (int ki = 0; ki < 3; ki++)
										// {
										// 		interpolateVariable[ki] = 0.0;
										// 		for (int vertex = 0; vertex < 3; vertex++)
										// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
										// };

										xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										//node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										//node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

										delete [] interpolateVariable; interpolateVariable = NULL;
								}
								else if (element[e].boundaryType == 'C')
								{
										int j = element[e].boundaryVertex;

										double* interpolateVariable;
										interpolateVariable = new double [3];

										/*********************************************/
										// Boundary on vertex (j + 2) % 3
										int i = element[e].globalNode[(j + 2) % 3];

										// for (int ki = 0; ki < 3; ki++)
										// {
										// 		interpolateVariable[ki] = 0.0;
										// 		for (int vertex = 0; vertex < 3; vertex++)
										// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
										// };

										double xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										double yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										double xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										double yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;

										node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

										/*********************************************/
										// Boundary on vertex (j + 1) % 3
										i = element[e].globalNode[(j + 1) % 3];

										// for (int ki = 0; ki < 3; ki++)
										// {
										// 		interpolateVariable[ki] = 0.0;
										// 		for (int vertex = 0; vertex < 3; vertex++)
										// 				interpolateVariable[ki] = interpolateVariable[ki] + node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										// 		interpolateVariable[ki] = interpolateVariable[ki] / 3.0;
										// };

										xCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
										yCoordinate = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + 3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

										for (int ki = 0; ki < 3; ki++)
										{
												interpolateVariable[ki] = 0.0;
												for (int vertex = 0; vertex < 3; vertex++)
														interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
										};

										xEdgeNormal = (node[element[e].globalNode[(j + 2) % 3]].coordinate[1] - node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 2.0;
										yEdgeNormal = (node[element[e].globalNode[(j + 1) % 3]].coordinate[0] - node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 2.0;


										node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
										node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

										delete [] interpolateVariable; interpolateVariable = NULL;
								};
            };
            break;
	};
}

void Computation2D::addResidualDistributionBoundaryFlux(const int & UVARIABLE_LEVEL)
{
	switch (TMTEMode)
	{
		case 'A':
	   		for (int e = 0; e < elementNumber; e++)
            {
				if (element[e].boundaryType == 'B')
				{
					int j = element[e].boundaryVertex;

					double* interpolateVariable;
					interpolateVariable = new double [3];

					/*********************************************/
					// Boundary on vertex j and (j + 1) % 3
					int i = element[e].globalNode[j];
					double xCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
					double yCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					double xEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][0] / 2.0;
					double yEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

					////////////////////////
					// Boundary on vertex j and (j + 2) % 3
					xCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[0] + node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 4.0;
					yCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[1] + node[element[e].globalNode[(j + 2) % 3]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					xEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][0] / 2.0;
					yEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

					/*********************************************/
					// Boundary on vertex (j + 1) % 3
					i = element[e].globalNode[(j + 1) % 3];
					xCoordinate = (3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0] + node[element[e].globalNode[j]].coordinate[0]) / 4.0;
					yCoordinate = (3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1] + node[element[e].globalNode[j]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					xEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][0] / 2.0;
					yEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];


					/*********************************************/
					// Boundary on vertex (j + 2) % 3
					i = element[e].globalNode[(j + 2) % 3];
					xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[j]].coordinate[0]) / 4.0;
					yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[j]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					xEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][0] / 2.0;
					yEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (yEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (- xEdgeNormal / permeability) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- xEdgeNormal / permittivity) * interpolateVariable[1] + (yEdgeNormal / permittivity) * interpolateVariable[0];

					delete [] interpolateVariable; interpolateVariable = NULL;
				}
            };
		    break;
        case 'B':
			for (int e = 0; e < elementNumber; e++)
	        {
				if (element[e].boundaryType == 'B')
				{
					int j = element[e].boundaryVertex;

					double* interpolateVariable;
					interpolateVariable = new double [3];

					/*********************************************/
					// Boundary on vertex j and (j + 1) % 3
					int i = element[e].globalNode[j];
					double xCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[0] + node[element[e].globalNode[(j + 1) % 3]].coordinate[0]) / 4.0;
					double yCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[1] + node[element[e].globalNode[(j + 1) % 3]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					double xEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][0] / 2.0;
					double yEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

					////////////////////////
					// Boundary on vertex j and (j + 2) % 3
					xCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[0] + node[element[e].globalNode[(j + 2) % 3]].coordinate[0]) / 4.0;
					yCoordinate = (3.0 * node[element[e].globalNode[j]].coordinate[1] + node[element[e].globalNode[(j + 2) % 3]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					xEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][0] / 2.0;
					yEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

					/*********************************************/
					// Boundary on vertex (j + 1) % 3
					i = element[e].globalNode[(j + 1) % 3];
					xCoordinate = (3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[0] + node[element[e].globalNode[j]].coordinate[0]) / 4.0;
					yCoordinate = (3.0 * node[element[e].globalNode[(j + 1) % 3]].coordinate[1] + node[element[e].globalNode[j]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					xEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][0] / 2.0;
					yEdgeNormal = - element[e].inwardNormal[(j + 2) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

					/*********************************************/
					// Boundary on vertex (j + 2) % 3
					i = element[e].globalNode[(j + 2) % 3];
					xCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[0] + node[element[e].globalNode[j]].coordinate[0]) / 4.0;
					yCoordinate = (3.0 * node[element[e].globalNode[(j + 2) % 3]].coordinate[1] + node[element[e].globalNode[j]].coordinate[1]) / 4.0;

					for (int ki = 0; ki < 3; ki++)
					{
						interpolateVariable[ki] = 0.0;
						for (int vertex = 0; vertex < 3; vertex++)
							interpolateVariable[ki] = interpolateVariable[ki] + LagrangeInterpolation(xCoordinate, yCoordinate, e, vertex) * node[element[e].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
					};

					xEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][0] / 2.0;
					yEdgeNormal = - element[e].inwardNormal[(j + 1) % 3][1] / 2.0;

					node[i].fluxResidual[0] = node[i].fluxResidual[0] + (- yEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[1] = node[i].fluxResidual[1] + (xEdgeNormal / permittivity) * interpolateVariable[2];
					node[i].fluxResidual[2] = node[i].fluxResidual[2] + (- yEdgeNormal / permeability) * interpolateVariable[0] + (xEdgeNormal / permeability) * interpolateVariable[1];

					delete [] interpolateVariable; interpolateVariable = NULL;
				}
            };
            break;
	};
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
		{
			for (int ki = 0; ki < 3; ki++)
				timeDependentSolution(i, 2, time);
		}
}

void Computation2D::finiteVolumeNodalUpdate()
{
	/*/ Runge-Kutta Method Stage-1 /**/
	interpolateMedianDualCenter(1);
	calculateFiniteVolume(1);
	calculateFiniteVolumeBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * node[i].fluxResidual[ki];
	time -= timeDelta / 2.0;
	boundaryCondition(2);

	/*/ Runge-Kutta Method Stage-2 /**/
	interpolateMedianDualCenter(2);
	calculateFiniteVolume(2);
	calculateFiniteVolumeBoundaryFlux(2);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * node[i].fluxResidual[ki];
	time += timeDelta / 2.0;
	boundaryCondition(2);

	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation2D::fluxDifferenceNodalUpdate()
{
	/*/ Runge-Kutta Method Stage-1 /**/
	calculateFluxDifference(1);
	calculateFluxDifferenceBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * node[i].fluxResidual[ki];
	time -= timeDelta / 2.0;
	boundaryCondition(2);

	/*/ Runge-Kutta Method Stage-2 /**/
	calculateFluxDifference(2);
	calculateFluxDifferenceBoundaryFlux(2);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * node[i].fluxResidual[ki];
	time += timeDelta / 2.0;
	boundaryCondition(2);

	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation2D::RDGalerkinNodalUpdate()
{
	/*/ Runge-Kutta Method Stage-1 /**/
	calculateRDGalerkin(1);
	addResidualDistributionBoundaryFlux(1);
	calculateFluxDifferenceBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].nodeArea) * node[i].fluxResidual[ki];
	time -= timeDelta / 2.0;
	boundaryCondition(2);

	/*/ Runge-Kutta Method Stage-2 /**/
	calculateRDGalerkin(2);
	addResidualDistributionBoundaryFlux(2);
	calculateFluxDifferenceBoundaryFlux(2);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * node[i].fluxResidual[ki];
	time += timeDelta / 2.0;
	boundaryCondition(2);

	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation2D::LaxWendroffNodalUpdate()
{
	/*/ Runge-Kutta Method Stage-1 /**/
	calculateLaxWendroff(1);
	addResidualDistributionBoundaryFlux(1);
	calculateFluxDifferenceBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 3; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].nodeArea) * node[i].fluxResidual[ki];

	boundaryCondition(2);

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
            results.open("./Output/Continuous_Finite_Volume_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
        case 'B':
            results.open("./Output/Flux_Difference_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
        case 'C':
            results.open("./Output/RD_Galerkin_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
            break;
        case 'D':
            results.open("./Output/RD_Lax_Wendroff_" + to_string(elementNumber) + "_" + stringTime + ".vtk");
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
        results.open("./Output/Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_Exact.txt");
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
            results.open("./Output/Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_ContinuousFV.txt");
            break;
        case 'B':
            results.open("./Output/Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_FluxDifference.txt");
            break;
        case 'C':
            results.open("./Output/Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_RDGalerkin.txt");
            break;
        case 'D':
            results.open("./Output/Wedge_Scattering_" + to_string(elementNumber) + "_" + stringTime + "_RDLW.txt");
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
                L2Errors.open("./Output/L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_ContinuousFV.txt");
                break;
            case 'B':
                L2Errors.open("./Output/L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_FluxDifference.txt");
                break;
            case 'C':
                L2Errors.open("./Output/L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_RDGalerkin.txt");
                break;
            case 'D':
                L2Errors.open("./Output/L2_Errors_" + to_string(elementNumber) + "_" + stringTime + "_RDLW.txt");
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
	element = new elementArray [elementNumber];
	for (int e = 0; e < elementNumber; e++)
	{
		element[e].globalNode = new int [3];

		if (method == 'A')
		{
			element[e].medianDualNormal = new double* [3];
			element[e].edgeTangent = new double* [3];
			for (int j = 0; j < 3; j++)
			{
				element[e].medianDualNormal[j] = new double [2];
				element[e].edgeTangent[j] = new double [2];
			};

			element[e].interpolateVariable = new double* [3];
			for (int j = 0; j < 3; j++)
				element[e].interpolateVariable[j] = new double [3];
		}
		else if (method == 'B' || method == 'C' || method == 'D')
		{
			element[e].inwardNormal = new double* [3];
			for (int j = 0; j < 3; j++)
				element[e].inwardNormal[j] = new double [2];
		};

		if (method == 'D')
		{
			element[e].distributionMatrix = new double** [3];
			for (int j = 0; j < 3; j++)
			{
				element[e].distributionMatrix[j] = new double* [3];
				for (int ki = 0; ki < 3; ki++)
					element[e].distributionMatrix[j][ki] = new double [3];
			}
		};
	};

	node = new nodeArray [nodeNumber];
	for (int i = 0; i < nodeNumber; i++)
	{
		node[i].coordinate = new double [2];

		node[i].conservedVariable = new double* [3];
		for (int ki = 0; ki < 3; ki++)
		{
			node[i].conservedVariable[ki] = new double [8];
			//node[i].conservedVariable[ki] = new double [6];
		};

		node[i].fluxResidual = new double [3];
	};
}

void Computation2D::timeIterations()
{
	initializeArray();

	if (method == 'A')
    {
        calculateEdgeTangent();
        calculateMedianDualNormal();
    }
    else if (method == 'B' || method == 'C' || method == 'D')
    {
        constructInwardNormal();
    };
	medianCellArea();

	timeLast = 2.0;
	timeDelta = globalTimeStep();
    timeDelta = 0.0001;
	timeNumber = timeLast / timeDelta;

	// timeNumber = static_cast <int> ((timeLast / 4.0) / timeDelta + 1) * 4;
	// timeDelta = static_cast<double> (timeLast / timeNumber);
	time = 0.0;


	// must be after timeDelta being defined
	if (method == 'D')
    {
        constructDistributionMatrix();
    };

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
            outputTime.open("./Output/Time Record (Continuous FV) " + to_string(elementNumber) + ".txt");
            break;
        case 'B':
            outputTime.open("./Output/Time Record (Flux-Difference) " + to_string(elementNumber) + ".txt");
            break;
        case 'C':
            outputTime.open("./Output/Time Record (RD-Galerkin) " + to_string(elementNumber) + ".txt");
            break;
        case 'D':
            outputTime.open("./Output/Time Record (RDLW) " + to_string(elementNumber) + ".txt");
            break;
	};
	STARTTIME = clock();

	if (method == 'A')
	{
		double* sum;
        sum = new double [3];

		for (int t = 1; t <= timeNumber; t++)
        {
            time = t * timeDelta;
			cout << "Time is " << time << endl;
            finiteVolumeNodalUpdate();

            intervalResults(time05, time1, time15, time2);
        };

		delete [] sum; sum = NULL;
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
	}
    else if (method == 'C')
	{
        for (int t = 1; t <= timeNumber; t++)
        {
            time = t * timeDelta;
            cout << "Time is " << time << endl;
            RDGalerkinNodalUpdate();

            intervalResults(time05, time1, time15, time2);
        };
	}
    else if (method == 'D')
	{
        for (int t = 1; t <= timeNumber; t++)
        {
            time = t * timeDelta;
            cout << "Time is " << time << endl;
            LaxWendroffNodalUpdate();

            intervalResults(time05, time1, time15, time2);
        };
	};
	outputTime << "The global time step is " << timeDelta << endl;
	outputTime << "The execution time is " << static_cast <double> ((clock() - STARTTIME) / static_cast <double> (CLOCKS_PER_SEC)) << ". " << endl;
	outputTime.close();
	errorsCalculation();
}

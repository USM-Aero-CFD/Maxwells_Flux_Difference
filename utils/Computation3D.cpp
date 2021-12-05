# include "Computation3D.h"

Computation3D::~Computation3D()
{
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
				for (int j = 0; j < 3; j++)
				{
					delete [] tetrahedron[h].inwardNormal[j]; tetrahedron[h].inwardNormal[j] = NULL;
				};
				delete [] tetrahedron[h].inwardNormal; tetrahedron[h].inwardNormal = NULL;
			};

			for (int mediandual = 0; mediandual < 6; mediandual++)
			{
				delete [] tetrahedron[h].LagrangeInterpolation[mediandual]; tetrahedron[h].LagrangeInterpolation[mediandual] = NULL;
			};
			delete [] tetrahedron[h].LagrangeInterpolation; tetrahedron[h].LagrangeInterpolation = NULL;
		}
		else if (method == 'B' || method == 'C' || method == 'D')
		{
			for (int j = 0; j < 4; j++)
			{
				delete [] tetrahedron[h].inwardNormal[j]; tetrahedron[h].inwardNormal[j] = NULL;
			};
			delete [] tetrahedron[h].inwardNormal; tetrahedron[h].inwardNormal = NULL;
		};

		if (method == 'D')
		{
			for (int j = 0; j < 4; j++)
			{
				for (int ki = 0; ki < 6; ki++)
				{
					delete [] tetrahedron[h].distributionMatrix[j][ki]; tetrahedron[h].distributionMatrix[j][ki] = NULL;
				};
				delete [] tetrahedron[h].distributionMatrix[j]; tetrahedron[h].distributionMatrix[j] = NULL;
			};
			delete [] tetrahedron[h].distributionMatrix; tetrahedron[h].distributionMatrix = NULL;
		};

		delete [] tetrahedron[h].subVolume; tetrahedron[h].subVolume = NULL;
	};

	for (int i = 0; i < nodeNumber; i++)
	{
		delete [] node[i].fluxResidual; node[i].fluxResidual = NULL;

		for (int ki = 0; ki < 6; ki++)
		{
			delete [] node[i].conservedVariable[ki]; node[i].conservedVariable[ki] = NULL;
		};
		delete [] node[i].conservedVariable; node[i].conservedVariable = NULL;
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

        complexNumber ERIncident;
        complexNumber ERScattered;
        complexNumber EThetaIncident;
        complexNumber EThetaScattered;
        complexNumber EPhiIncident;
        complexNumber EPhiScattered;
        complexNumber HRIncident;
        complexNumber HRScattered;
        complexNumber HThetaIncident;
        complexNumber HThetaScattered;
        complexNumber HPhiIncident;
        complexNumber HPhiScattered;

        ERIncident.real = 0.0;
        ERIncident.imaginary = 0.0;
        ERScattered.real = 0.0;
        ERScattered.imaginary = 0.0;
        EThetaIncident.real = 0.0;
        EThetaIncident.imaginary = 0.0;
        EThetaScattered.real = 0.0;
        EThetaScattered.imaginary = 0.0;
        EPhiIncident.real = 0.0;
        EPhiIncident.imaginary = 0.0;
        EPhiScattered.real = 0.0;
        EPhiScattered.imaginary = 0.0;
        HRIncident.real = 0.0;
        HRIncident.imaginary = 0.0;
        HRScattered.real = 0.0;
        HRScattered.imaginary = 0.0;
        HThetaIncident.real = 0.0;
        HThetaIncident.imaginary = 0.0;
        HThetaScattered.real = 0.0;
        HThetaScattered.imaginary = 0.0;
        HPhiIncident.real = 0.0;
        HPhiIncident.imaginary = 0.0;
        HPhiScattered.real = 0.0;
        HPhiScattered.imaginary = 0.0;

        EPhiIncident.real = cos(- propagationCoefficient * xCoordinate);
        EPhiIncident.imaginary = sin(- propagationCoefficient * xCoordinate);
        HThetaIncident.real = (- propagationCoefficient / (angularFrequency * permeability)) * cos(- propagationCoefficient * xCoordinate);
        HThetaIncident.imaginary = (- propagationCoefficient / (angularFrequency * permeability)) * sin(- propagationCoefficient * xCoordinate);

        node[i].ERIncident = ERIncident;
        node[i].ERScattered = ERScattered;
        node[i].EThetaIncident = EThetaIncident;
        node[i].EThetaScattered = EThetaScattered;
        node[i].EPhiIncident = EPhiIncident;
        node[i].EPhiScattered = EPhiScattered;
        node[i].HRIncident = HRIncident;
        node[i].HRScattered = HRScattered;
        node[i].HThetaIncident = HThetaIncident;
        node[i].HThetaScattered = HThetaScattered;
        node[i].HPhiIncident = HPhiIncident;
        node[i].HPhiScattered = HPhiScattered;
    }
}

void Computation3D::spatialSolution()
{
    switch (TMTEMode)
    {
        case 'A':
            TMmodeSolution();
            break;
        // case 'B':
        //     TEmodeSolution();
        //     break;
    };
}

void Computation3D::timeDependentSolution(const int & I, const int & UVARIABLE_LEVEL, const double & Time) const
{
    complexNumber timeHarmonic;
    complexNumber ER, ETheta, EPhi;
    complexNumber HR, HTheta, HPhi;

    timeHarmonic.real = cos(angularFrequency * Time);
    timeHarmonic.imaginary = sin(angularFrequency * Time);

    ER = complexMultiplication(timeHarmonic, node[I].ERIncident);
    ETheta = complexMultiplication(timeHarmonic, node[I].EThetaIncident);
    EPhi = complexMultiplication(timeHarmonic, node[I].EPhiIncident);
    HR = complexMultiplication(timeHarmonic, node[I].HRIncident);
    HTheta = complexMultiplication(timeHarmonic, node[I].HThetaIncident);
    HPhi = complexMultiplication(timeHarmonic, node[I].HPhiIncident);

    node[I].conservedVariable[0][UVARIABLE_LEVEL] = ER.real;
    node[I].conservedVariable[1][UVARIABLE_LEVEL] = ETheta.real;
    node[I].conservedVariable[2][UVARIABLE_LEVEL] = EPhi.real;
    node[I].conservedVariable[3][UVARIABLE_LEVEL] = HR.real;
    node[I].conservedVariable[4][UVARIABLE_LEVEL] = HTheta.real;
    node[I].conservedVariable[5][UVARIABLE_LEVEL] = HPhi.real;
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

void Computation3D::calculateLagrangeInterpolation()
{
	for (int h = 0; h < tetrahedronNumber; h++)
		for (int mediandual = 0; mediandual < 6; mediandual++)
		{
			double* medianDualCentroid;
			medianDualCentroid = new double [3];

			switch (mediandual)
			{
				case 0:
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate]) / 2.0) / 4.0;
					break;
				case 1:
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0) / 4.0;
					break;
				case 2:
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0) / 4.0;
					break;
				case 3:
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate]) / 2.0) / 4.0;
					break;
				case 4:
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0) / 4.0;
					break;
				case 5:
					for (int coordinate = 0; coordinate < 3; coordinate++)
						medianDualCentroid[coordinate] = ((node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[2]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 4.0
															+ (node[tetrahedron[h].globalNode[0]].coordinate[coordinate] + node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 3.0
															+ (node[tetrahedron[h].globalNode[1]].coordinate[coordinate] + node[tetrahedron[h].globalNode[3]].coordinate[coordinate]) / 2.0) / 4.0;
					break;
			};

			for (int j = 0; j < 4; j++)
				tetrahedron[h].LagrangeInterpolation[mediandual][j] = LagrangeInterpolationFunction(medianDualCentroid, h, j);

			delete [] medianDualCentroid; medianDualCentroid = NULL;
		}
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
	int ipositive, inegative;
    double* fluxResidual;
    double* interpolateVariable;
    fluxResidual = new double [6];
	interpolateVariable = new double [6];


    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 6; ki++)
			node[i].fluxResidual[ki] = 0.0;


    for (int h = 0; h < tetrahedronNumber; h++)
    {
        for (int mediandual = 0; mediandual < 6; mediandual++)
        {
            for (int ki = 0; ki < 6; ki++)
                fluxResidual[ki] = 0.0;

            for (int ki = 0; ki < 6; ki++)
            {
                interpolateVariable[ki] = 0.0;
                for (int vertex = 0; vertex < 4; vertex++)
                    interpolateVariable[ki] = interpolateVariable[ki] + tetrahedron[h].LagrangeInterpolation[mediandual][vertex] * node[tetrahedron[h].globalNode[vertex]].conservedVariable[ki][UVARIABLE_LEVEL];
            };

            switch (mediandual)
            {
                case 0:
                    ipositive = tetrahedron[h].globalNode[1];
                    inegative = tetrahedron[h].globalNode[0];
                    break;
                case 1:
                    ipositive = tetrahedron[h].globalNode[2];
                    inegative = tetrahedron[h].globalNode[0];
                    break;
                case 2:
                    ipositive = tetrahedron[h].globalNode[3];
                    inegative = tetrahedron[h].globalNode[0];
                    break;
                case 3:
                    ipositive = tetrahedron[h].globalNode[2];
                    inegative = tetrahedron[h].globalNode[1];
                    break;
                case 4:
                    ipositive = tetrahedron[h].globalNode[3];
                    inegative = tetrahedron[h].globalNode[2];
                    break;
                case 5:
                    ipositive = tetrahedron[h].globalNode[1];
                    inegative = tetrahedron[h].globalNode[3];
                    break;
            };

            fluxResidual[0] = (interpolateVariable[4] * tetrahedron[h].medianDualNormal[mediandual][2] / permittivity - interpolateVariable[5] * tetrahedron[h].medianDualNormal[mediandual][1] / permittivity);
            fluxResidual[1] = (interpolateVariable[5] * tetrahedron[h].medianDualNormal[mediandual][0] / permittivity - interpolateVariable[3] * tetrahedron[h].medianDualNormal[mediandual][2] / permittivity);
            fluxResidual[2] = (interpolateVariable[3] * tetrahedron[h].medianDualNormal[mediandual][1] / permittivity - interpolateVariable[4] * tetrahedron[h].medianDualNormal[mediandual][0] / permittivity);
            fluxResidual[3] = (interpolateVariable[2] * tetrahedron[h].medianDualNormal[mediandual][1] / permeability - interpolateVariable[1] * tetrahedron[h].medianDualNormal[mediandual][2] / permeability);
            fluxResidual[4] = (interpolateVariable[0] * tetrahedron[h].medianDualNormal[mediandual][2] / permeability - interpolateVariable[2] * tetrahedron[h].medianDualNormal[mediandual][0] / permeability);
            fluxResidual[5] = (interpolateVariable[1] * tetrahedron[h].medianDualNormal[mediandual][0] / permeability - interpolateVariable[0] * tetrahedron[h].medianDualNormal[mediandual][1] / permeability);

            for (int ki = 0; ki < 6; ki++)
                node[ipositive].fluxResidual[ki] = node[ipositive].fluxResidual[ki] + fluxResidual[ki];

            for (int ki = 0; ki < 6; ki++)
                node[inegative].fluxResidual[ki] = node[inegative].fluxResidual[ki] - fluxResidual[ki];
        };
    };

	delete [] fluxResidual; fluxResidual = NULL;
	delete [] interpolateVariable; interpolateVariable = NULL;
}

void Computation3D::fluxDifference(const int & H, const int & UVARIABLE_LEVEL)
{
    double xElectricField = 0.0;
    double yElectricField = 0.0;
    double zElectricField = 0.0;
    double xMagneticField = 0.0;
    double yMagneticField = 0.0;
    double zMagneticField = 0.0;

    for (int j = 0; j < 4; j++)
    {
        xElectricField += node[tetrahedron[H].globalNode[j]].conservedVariable[0][UVARIABLE_LEVEL];
        yElectricField += node[tetrahedron[H].globalNode[j]].conservedVariable[1][UVARIABLE_LEVEL];
        zElectricField += node[tetrahedron[H].globalNode[j]].conservedVariable[2][UVARIABLE_LEVEL];
        xMagneticField += node[tetrahedron[H].globalNode[j]].conservedVariable[3][UVARIABLE_LEVEL];
        yMagneticField += node[tetrahedron[H].globalNode[j]].conservedVariable[4][UVARIABLE_LEVEL];
        zMagneticField += node[tetrahedron[H].globalNode[j]].conservedVariable[5][UVARIABLE_LEVEL];
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
}

void Computation3D::calculateFluxDifference(const int & UVARIABLE_LEVEL)
{
    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 6; ki++)
            node[i].fluxResidual[ki] = 0.0;

    for (int h = 0; h < tetrahedronNumber; h++)
	  {
        fluxDifference(h, UVARIABLE_LEVEL);
	  };
}

// both RD-Galerkin and RD-LW
void Computation3D::localFluxResidual(const int & H, const int & UVARIABLE_LEVEL, double* & fluxResidual) const
{
    for (int j = 0; j < 4; j++)
    {
        int i = tetrahedron[H].globalNode[j];

        double xElectricField = node[i].conservedVariable[0][UVARIABLE_LEVEL];
        double yElectricField = node[i].conservedVariable[1][UVARIABLE_LEVEL];
        double zElectricField = node[i].conservedVariable[2][UVARIABLE_LEVEL];
        double xMagneticField = node[i].conservedVariable[3][UVARIABLE_LEVEL];
        double yMagneticField = node[i].conservedVariable[4][UVARIABLE_LEVEL];
        double zMagneticField = node[i].conservedVariable[5][UVARIABLE_LEVEL];
        double xInwardNormal = tetrahedron[H].inwardNormal[j][0];
        double yInwardNormal = tetrahedron[H].inwardNormal[j][1];
        double zInwardNormal = tetrahedron[H].inwardNormal[j][2];

        if (tetrahedron[H].boundaryType == 'A'
            || tetrahedron[H].boundaryType == 'C'
            || tetrahedron[H].boundaryType == 'D'
            || tetrahedron[H].boundaryType == 'N')
        {
            fluxResidual[0] = fluxResidual[0] + (1.0 / 3.0) * ((yMagneticField / permittivity) * zInwardNormal - (zMagneticField / permittivity) * yInwardNormal);
            fluxResidual[1] = fluxResidual[1] + (1.0 / 3.0) * ((zMagneticField / permittivity) * xInwardNormal - (xMagneticField / permittivity) * zInwardNormal);
            fluxResidual[2] = fluxResidual[2] + (1.0 / 3.0) * ((xMagneticField / permittivity) * yInwardNormal - (yMagneticField / permittivity) * xInwardNormal);
            fluxResidual[3] = fluxResidual[3] + (1.0 / 3.0) * ((zElectricField / permittivity) * yInwardNormal - (yElectricField / permittivity) * zInwardNormal);
            fluxResidual[4] = fluxResidual[4] + (1.0 / 3.0) * ((xElectricField / permittivity) * zInwardNormal - (zElectricField / permittivity) * xInwardNormal);
            fluxResidual[5] = fluxResidual[5] + (1.0 / 3.0) * ((yElectricField / permittivity) * xInwardNormal - (xElectricField / permittivity) * yInwardNormal);
        };;
    };
}

void Computation3D::RDGalerkin(const int & H, const int & UVARIABLE_LEVEL)
{
    double* fluxResidual;
    fluxResidual = new double [6];
    for (int ki = 0; ki < 6; ki++)
        fluxResidual[ki] = 0.0;

    localFluxResidual(H, UVARIABLE_LEVEL, fluxResidual);

    for (int j = 0; j < 4; j++)
    {
        int i = tetrahedron[H].globalNode[j];

        for (int ki = 0; ki < 6; ki++)
            node[i].fluxResidual[ki] = node[i].fluxResidual[ki] + (1.0 / 4.0) * fluxResidual[ki];
    };

    delete [] fluxResidual; fluxResidual = NULL;
}

void Computation3D::calculateRDGalerkin(const int & UVARIABLE_LEVEL)
{
    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 6; ki++)
            node[i].fluxResidual[ki] = 0.0;

    switch (TMTEMode)
    {
        case 'A':
        case 'B':
            for (int h = 0; h < tetrahedronNumber; h++)
            {
                if (tetrahedron[h].boundaryType == 'A'
                    || tetrahedron[h].boundaryType == 'C'
                    || tetrahedron[h].boundaryType == 'D'
                    || tetrahedron[h].boundaryType == 'N')
                {
                    RDGalerkin(h, UVARIABLE_LEVEL);
                }
                else if (tetrahedron[h].boundaryType == 'B')
                {
                    // RDGALERKIN(h, UVARIABLE_LEVEL);
                    fluxDifference(h, UVARIABLE_LEVEL);
                };
            };
            break;
    };
}

void Computation3D::LaxWendroff(const int & H, const int & UVARIABLE_LEVEL)
{
    double* fluxResidual;
    fluxResidual = new double [6];
    for (int ki = 0; ki < 6; ki++)
        fluxResidual[ki] = 0.0;

    localFluxResidual(H, UVARIABLE_LEVEL, fluxResidual);

    for (int j = 0; j < 4; j++)
    {
        int i = tetrahedron[H].globalNode[j];

        for (int ki = 0; ki < 6; ki++)
            for (int kj = 0; kj < 6; kj++)
                node[i].fluxResidual[ki] = node[i].fluxResidual[ki] + tetrahedron[H].distributionMatrix[j][ki][kj] * fluxResidual[kj];
    };

    delete [] fluxResidual; fluxResidual = NULL;
}

void Computation3D::calculateLaxWendroff(const int & UVARIABLE_LEVEL)
{
    for (int i = 0; i < nodeNumber; i++)
        for (int ki = 0; ki < 6; ki++)
            node[i].fluxResidual[ki] = 0.0;

    switch (TMTEMode)
    {
        case 'A':
        case 'B':
            for (int h = 0; h < tetrahedronNumber; h++)
            {
                if (tetrahedron[h].boundaryType == 'A'
                    || tetrahedron[h].boundaryType == 'C'
                    || tetrahedron[h].boundaryType == 'D'
                    || tetrahedron[h].boundaryType == 'N')
                {
                    LaxWendroff(h, UVARIABLE_LEVEL);
                }
                else if (tetrahedron[h].boundaryType == 'B')
                {
                    // LAXWENDROFF(h, UVARIABLE_LEVEL);
                    fluxDifference(h, UVARIABLE_LEVEL);
                };
            };
            break;
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

              double xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][0];
              double yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][1];
              double zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][2];

              for (int ix = 0; ix < 3; ix++)
                  centroid[ix] = (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                          + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                          + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 3.0;

              /////////////////////////////////////////////////////
              // (j + 1) % 4
              for (int ix = 0; ix < 3; ix++)
                  coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                              + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                              + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

              for (int ki = 0; ki < 6; ki++)
                  interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

              boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              if (tetrahedron[h].boundaryType == 'B')
              {
                  for (int ki = 3; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] + fluxResidual[ki];
              }
              else if (tetrahedron[h].boundaryType == 'C')
              {
                  for (int ki = 0; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] + fluxResidual[ki];
              };

              /////////////////////////////////////////////////////
              // (j + 2) % 4
              for (int ix = 0; ix < 3; ix++)
                  coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                              + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                              + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

              for (int ki = 0; ki < 6; ki++)
                  interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

              boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              if (tetrahedron[h].boundaryType == 'B')
              {
                  for (int ki = 3; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] + fluxResidual[ki];
              }
              else if (tetrahedron[h].boundaryType == 'C')
              {
                  for (int ki = 0; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] + fluxResidual[ki];
              };

              /////////////////////////////////////////////////////
              // (j + 3) % 4
              for (int ix = 0; ix < 3; ix++)
                  coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]
                                              + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                              + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

              for (int ki = 0; ki < 6; ki++)
                  interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

              boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              if (tetrahedron[h].boundaryType == 'B')
              {
                  for (int ki = 3; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] + fluxResidual[ki];
              }
              else if (tetrahedron[h].boundaryType == 'C')
              {
                  for (int ki = 0; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] + fluxResidual[ki];
              };
          }
        };
        break;
    };

    delete [] fluxResidual; fluxResidual = NULL;
    delete [] coordinate; coordinate = NULL;
    delete [] centroid; centroid = NULL;
    delete [] interpolateVariable; interpolateVariable = NULL;
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

              double xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][0];
              double yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][1];
              double zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[j][2];

              for (int ix = 0; ix < 3; ix++)
                  centroid[ix] = (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                          + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                          + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 3.0;

              /////////////////////////////////////////////////////
              // (j + 1) % 4
              for (int ix = 0; ix < 3; ix++)
                  coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                              + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                              + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

              for (int ki = 0; ki < 6; ki++)
                  interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

              boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              if (tetrahedron[h].boundaryType == 'B')
              {
                  for (int ki = 3; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] + fluxResidual[ki];
              }
              else if (tetrahedron[h].boundaryType == 'C')
              {
                  for (int ki = 0; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] + fluxResidual[ki];
              };

              /////////////////////////////////////////////////////
              // (j + 2) % 4
              for (int ix = 0; ix < 3; ix++)
                  coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                              + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                              + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

              for (int ki = 0; ki < 6; ki++)
                  interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

              boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              if (tetrahedron[h].boundaryType == 'B')
              {
                  for (int ki = 3; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] + fluxResidual[ki];
              }
              else if (tetrahedron[h].boundaryType == 'C')
              {
                  for (int ki = 0; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] + fluxResidual[ki];
              };

              /////////////////////////////////////////////////////
              // (j + 3) % 4
              for (int ix = 0; ix < 3; ix++)
                  coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]
                                              + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                              + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                              + centroid[ix]) / 4.0;

              for (int ki = 0; ki < 6; ki++)
                  interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                             + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

              boundaryFluxResidual(fluxResidual, interpolateVariable,
                                    xOutwardNormal, yOutwardNormal, zOutwardNormal);

              if (tetrahedron[h].boundaryType == 'B')
              {
                  for (int ki = 3; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] + fluxResidual[ki];
              }
              else if (tetrahedron[h].boundaryType == 'C')
              {
                  for (int ki = 0; ki < 6; ki++)
                      node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] + fluxResidual[ki];
              };
          }
        };
        break;
    };
	delete [] fluxResidual; fluxResidual = NULL;
    delete [] coordinate; coordinate = NULL;
    delete [] centroid; centroid = NULL;
    delete [] interpolateVariable; interpolateVariable = NULL;
}

void Computation3D::addResidualDistributionBoundaryFlux(const int & UVARIABLE_LEVEL)
{
    double* fluxResidual;
    fluxResidual = new double [6];

    double* coordinate;
    coordinate = new double [3];

    double* centroid;
    centroid = new double [3];

    double* interpolateVariable;
    interpolateVariable = new double [6];

    switch (TMTEMode)
    {
        case 'A':
        case 'B':
            for (int h = 0; h < tetrahedronNumber; h++)
            {
                if (tetrahedron[h].boundaryType == 'B')
                {
                    for (int ki = 0; ki < 6; ki++)
                        fluxResidual[ki] = 0.0;

                    int j = tetrahedron[h].boundaryVertex;


                    /////////////////////////////////////////////////////
                    // plane j, (j + 1) % 4, (j + 2) % 4
                    double xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 3) % 4][0];
                    double yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 3) % 4][1];
                    double zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 3) % 4][2];

                    for (int ix = 0; ix < 3; ix++)
                        centroid[ix] = (node[tetrahedron[h].globalNode[j]].coordinate[ix]
                                                + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                                + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 3.0;

                    /////////////////
                    // j
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[j]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[j]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[j]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[j]].fluxResidual[ki] = node[tetrahedron[h].globalNode[j]].fluxResidual[ki] + fluxResidual[ki];

                    /////////////////
                    // (j + 1) % 4
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[j]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] + fluxResidual[ki];

                    /////////////////
                    // (j + 2) % 4
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[j]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] + fluxResidual[ki];



                    /////////////////////////////////////////////////////
                    // plane j, (j + 2) % 4, (j + 3) % 4
                    xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 1) % 4][0];
                    yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 1) % 4][1];
                    zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 1) % 4][2];

                    for (int ix = 0; ix < 3; ix++)
                        centroid[ix] = (node[tetrahedron[h].globalNode[j]].coordinate[ix]
                                                + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                                + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 3.0;

                    /////////////////
                    // j
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[j]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[j]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[j]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[j]].fluxResidual[ki] = node[tetrahedron[h].globalNode[j]].fluxResidual[ki] + fluxResidual[ki];

                    /////////////////
                    // (j + 2) % 4
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[j]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 2) % 4]].fluxResidual[ki] + fluxResidual[ki];

                    /////////////////
                    // (j + 3) % 4
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[j]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 2) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 2) % 4) * node[tetrahedron[h].globalNode[(j + 2) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] + fluxResidual[ki];


                    /////////////////////////////////////////////////////
                    // plane j, (j + 3) % 4, (j + 1) % 4
                    xOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 2) % 4][0];
                    yOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 2) % 4][1];
                    zOutwardNormal = - (1.0 / 3.0) * tetrahedron[h].inwardNormal[(j + 2) % 4][2];

                    for (int ix = 0; ix < 3; ix++)
                        centroid[ix] = (node[tetrahedron[h].globalNode[j]].coordinate[ix]
                                                + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]
                                                + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 3.0;

                    /////////////////
                    // j
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[j]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[j]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[j]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[j]].fluxResidual[ki] = node[tetrahedron[h].globalNode[j]].fluxResidual[ki] + fluxResidual[ki];

                    /////////////////
                    // (j + 3) % 4
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[j]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 3) % 4]].fluxResidual[ki] + fluxResidual[ki];

                    /////////////////
                    // (j + 1) % 4
                    for (int ix = 0; ix < 3; ix++)
                        coordinate[ix] = (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix]
                                                + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[(j + 3) % 4]].coordinate[ix]) / 2.0
                                                + (node[tetrahedron[h].globalNode[(j + 1) % 4]].coordinate[ix] + node[tetrahedron[h].globalNode[j]].coordinate[ix]) / 2.0
                                                + centroid[ix]) / 4.0;

                    for (int ki = 0; ki < 6; ki++)
                        interpolateVariable[ki] = (LagrangeInterpolationFunction(coordinate, h, j) * node[tetrahedron[h].globalNode[j]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 3) % 4) * node[tetrahedron[h].globalNode[(j + 3) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]
                                                + LagrangeInterpolationFunction(coordinate, h, (j + 1) % 4) * node[tetrahedron[h].globalNode[(j + 1) % 4]].conservedVariable[ki][UVARIABLE_LEVEL]);

                    boundaryFluxResidual(fluxResidual, interpolateVariable, xOutwardNormal, yOutwardNormal, zOutwardNormal);

                    for (int ki = 0; ki < 6; ki++)
                        node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] = node[tetrahedron[h].globalNode[(j + 1) % 4]].fluxResidual[ki] + fluxResidual[ki];
                }
            };
            break;
    };

    delete [] fluxResidual; fluxResidual = NULL;
    delete [] coordinate; coordinate = NULL;
    delete [] centroid; centroid = NULL;
    delete [] interpolateVariable; interpolateVariable = NULL;
}

void Computation3D::boundaryCondition(const int & UVARIABLE_LEVEL)
{
	for (int i = 0; i < nodeNumber; i++)
	{
		if (node[i].boundary == 'A' || node[i].boundary == 'D' || node[i].boundary == 'C')
		{
            timeDependentSolution(i, UVARIABLE_LEVEL, time);
		};
	};

}

void Computation3D::finiteVolumeNodalUpdate()
{
	/*/ Runge-Kutta Stage-1 /**/
	calculateFiniteVolume(1);
	calculateFiniteVolumeBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].volume) * node[i].fluxResidual[ki];
	time -= timeDelta / 2.0;
	boundaryCondition(2);

	/*/ Runge-Kutta Stage-2 /**/
	calculateFiniteVolume(2);
	calculateFiniteVolumeBoundaryFlux(2);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * node[i].fluxResidual[ki];
	time += timeDelta / 2.0;
	boundaryCondition(2);

	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation3D::fluxDifferenceNodalUpdate()
{
	/*/ Runge-Kutta Stage-1 /**/
	calculateFluxDifference(1);
	calculateFluxDifferenceBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].volume) * node[i].fluxResidual[ki];
	time -= timeDelta / 2.0;
	boundaryCondition(2);

	/*/ Runge-Kutta Stage-2 /**/
	calculateFluxDifference(2);
	calculateFluxDifferenceBoundaryFlux(2);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * node[i].fluxResidual[ki];
	time += timeDelta / 2.0;
	boundaryCondition(2);

	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];
}

void Computation3D::RDGalerkinNodalUpdate()
{
	/*/ Runge-Kutta Stage-1 /**/
	calculateRDGalerkin(1);
	addResidualDistributionBoundaryFlux(1);
	calculateFluxDifferenceBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - ((timeDelta / 2.0) / node[i].volume) * node[i].fluxResidual[ki];
	time -= timeDelta / 2.0;
	boundaryCondition(2);

	/*/ Runge-Kutta Stage-2 /**/
	calculateRDGalerkin(2);
	addResidualDistributionBoundaryFlux(2);
	calculateFluxDifferenceBoundaryFlux(2);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * node[i].fluxResidual[ki];
	time += timeDelta / 2.0;
	boundaryCondition(2);

	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][1] = node[i].conservedVariable[ki][2];

}

void Computation3D::LaxWendroffNodalUpdate()
{
	calculateLaxWendroff(1);
    addResidualDistributionBoundaryFlux(1);
	calculateFluxDifferenceBoundaryFlux(1);
	for (int i = 0; i < nodeNumber; i++)
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki][2] = node[i].conservedVariable[ki][1] - (timeDelta / node[i].volume) * node[i].fluxResidual[ki];
	boundaryCondition(2);

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
        results.open("./Output/Spherical_Scattering_Cross_Section_" + stringTime + "_Exact.txt");
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
            results.open("./Output/Spherical_Scattering_Cross_Section_" + stringTime + "_ContinuousFV.txt");
            break;
        case 'B':
            results.open("./Output/Spherical_Scattering_Cross_Section_" + stringTime + "_FluxDifference.txt");
            break;
        case 'C':
            results.open("./Output/Spherical_Scattering_Cross_Section_" + stringTime + "_RDGalerkin.txt");
            break;
        case 'D':
            results.open("./Output/Spherical_Scattering_Cross_Section_" + stringTime + "_RDLW.txt");
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
            results.open("./Output/Spherical_Scattering_" + stringTime + "_ContinuousFV.txt");
            break;
        case 'B':
            results.open("./Output/Spherical_Scattering_" + stringTime + "_FluxDifference.txt");
            break;
        case 'C':
            results.open("./Output/Spherical_Scattering_" + stringTime + "_RDGalerkin.txt");
            break;
        case 'D':
            results.open("./Output/Spherical_Scattering_" + stringTime + "_RDLW.txt");
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
                L2Errors.open("./Output/L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_ContinuousFV.txt");
                break;
            case 'B':
                L2Errors.open("./Output/L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_FluxDifference.txt");
                break;
            case 'C':
                L2Errors.open("./Output/L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_RDGalerkin.txt");
                break;
            case 'D':
                L2Errors.open("./Output/L2_Errors_" + to_string(tetrahedronNumber) + "_" + stringTime + "_RDLW.txt");
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
	tetrahedron = new tetrahedronArray [tetrahedronNumber];
	for (int h = 0; h < tetrahedronNumber; h++)
	{
		tetrahedron[h].globalNode = new int [4];
		tetrahedron[h].medianDualNormal = new double* [6];
		for (int mediandual = 0; mediandual < 6; mediandual++)
		{
			tetrahedron[h].medianDualNormal[mediandual] = new double [3];
		};

		if (method == 'A')
		{
			tetrahedron[h].LagrangeInterpolation = new double* [6];
			for (int mediandual = 0; mediandual < 6; mediandual++)
			{
				tetrahedron[h].LagrangeInterpolation[mediandual] = new double [4];
			};
		}
		else if (method == 'B' || method == 'C' || method == 'D')
		{
			tetrahedron[h].inwardNormal = new double* [4];
			for (int j = 0; j < 4; j++)
			{
				tetrahedron[h].inwardNormal[j] = new double [3];
			};
		};

		if (method == 'D')
		{
			tetrahedron[h].distributionMatrix = new double** [4];
			for (int j = 0; j < 4; j++)
			{
				tetrahedron[h].distributionMatrix[j] = new double* [6];
				for (int ki = 0; ki < 6; ki++)
					tetrahedron[h].distributionMatrix[j][ki] = new double [6];
			};
		};

		tetrahedron[h].subVolume = new double [4];
	};

	node = new nodeArray [nodeNumber];
	for (int i = 0; i < nodeNumber; i++)
	{
		node[i].coordinate = new double [3];
		node[i].fluxResidual = new double [6];

		node[i].conservedVariable = new double* [6];
		for (int ki = 0; ki < 6; ki++)
			node[i].conservedVariable[ki] = new double [8];
	};


	crossSection = new crossSectionArray [crossSectionNumber];
	for (int crosssection = 0; crosssection < crossSectionNumber; crosssection++)
	{
		crossSection[crosssection].coordinate = new double [3];
	};
}

void Computation3D::timeComputations()
{
	initializeArray();

	constructMedianDualNormal();
	calculateSubVolume();
	medianCellVolume();
	findTetrahedronOrientation();

    if (method == 'A')
    {
        calculateLagrangeInterpolation();
        calculateInwardNormalFromMedianDualNormal();
    }
    else if (method == 'B' || method == 'C' || method == 'D')
    {
        calculateInwardNormal();
    };

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


	double timeLast = 0.2;

	if (method == 'A')
	{
	 	ofstream outputInterpolation;
	  	outputInterpolation.open("./Output/Interpolation.txt");
	  	for (int h = 0; h < tetrahedronNumber; h++)
	  	{
		  	outputInterpolation << "Tetrahedron is " << h << endl;
		  	for (int mediandual = 0; mediandual < 6; mediandual++)
		  	{
			 	for (int j = 0; j < 4; j++)
					outputInterpolation << setw(20) << tetrahedron[h].LagrangeInterpolation[mediandual][j];
			  	outputInterpolation << endl;
		  	};
		  	outputInterpolation << endl;
	  	};
	  	outputInterpolation.close();
	};


	// timeDelta = globalTimeStep();
	timeDelta = 0.0005;

	// to correlate the timeNumber exactly halt at timeLast / 4, timeLast / 2, 0.75 * timeLast, timeLast
	timeNumber = static_cast <int> ((timeLast / 4.0) / timeDelta) * 4;
	timeDelta = timeLast / timeNumber;

	time = 0.0;

    spatialSolution();


	// must be placed after timeDelta is being defined
    if (method == 'D')
        constructDistributionMatrix();
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
            outputTime.open("./Output/Time Record (Continuous FV) " + to_string(tetrahedronNumber) + ".txt");
            break;
        case 'B':
            outputTime.open("./Output/Time Record (Flux-Difference) " + to_string(tetrahedronNumber) + ".txt");
            break;
        case 'C':
            outputTime.open("./Output/Time Record (RD-Galerkin) " + to_string(tetrahedronNumber) + ".txt");
            break;
        case 'D':
            outputTime.open("./Output/Time Record (RDLW) " + to_string(tetrahedronNumber) + ".txt");
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
  /**/
}

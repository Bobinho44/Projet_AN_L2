#include "remez.hpp"


// Question 2
/*
	Create an array of nElements equidistant points in the interval [boundA, boundB].
	The array must contain at least two points (boudA and boundB).
	After the procedure, interpolationPoints is an array composed of the nElements points.
	We assume that interpolationPoints has been allocated outside the function.
*/
void createInterpolationPoints(double *interpolationPoints, double boundA, double boundB, uint64_t nElements) {

	// The array must contain at least two points (boudA and boundB)
    assert(nElements>1 && boundA < boundB);

	// Interval between two interpolation points
	double interval = (boundB - boundA) / (nElements - 1);

	// Create the array containing nElements interpolation points
    for (int i = 0; i < nElements; i++) {
        interpolationPoints[i] = boundA + i * interval;
    }
}


//  Question 3
/*
	Search for a root of the derivative of the function interpolation(x) - interpolated function(x) 
	on the [boundA, boundB] interval.
	We assume that poly has been allocated outside the function.
	Return a double:
		*The root of the derivative of the function interpolation(x) - interpolated function(x) 
		* on the [boundA, boundB] interval corresponds to an extremum of the function 
		* interpolation(x) - interpolated function(x) on that interval.
*/
double extremumLocal(const double* poly, double boundA, double boundB, fonction f, uint64_t polyDegree) {

	// Assignment of the limits of the interval to the first two terms of the Newton iteration
    double x1 = boundA;
    double x2 = boundB;
	double deriver_x1 = 1;
	double deriver_x2 = 2;

    // Search for a root of the derivative of the function interpolation(x) - interpolated function(x) 
	// such that the value of this root is accurate to within epsilon. If the algorithm does not find a root 
	// and diverges, it will stop after stop iterations.
    for (int i = 0; fabs(deriver_x2) > epsilon  && i < stop; i++) {

		// Calculates the x1 and x2 derivative of the the function interpolation(x) - interpolated function(x)
		deriver_x1 = ((horner(poly, x1 + epsilon, polyDegree) - f(x1 + epsilon)) - (horner(poly, x1 - epsilon, polyDegree) - f(x1 - epsilon))) / (2 * epsilon);
		deriver_x2 = ((horner(poly, x2 + epsilon, polyDegree) - f(x2 + epsilon)) - (horner(poly, x2 - epsilon, polyDegree) - f(x2 - epsilon))) / (2 * epsilon);

        // Calculation of the next term of the Newton iteration
        double resultat = x2 - deriver_x2 / ((deriver_x2 - deriver_x1) / (x2 - x1));
        x1 = x2;
        x2 = resultat;

		// The extremum point cannot exceed our chosen interval
        if (x2 > boundB) x2 = boundA;
    }
    return x2;
}


// Question 4
/*
	At the end of the procedure, extremumPoints is a an array containing the extremes of 
	the function interpolation(x) - interpolated function(x) on the interval 
	[interpolationPoints[0], interpolationPoints[polyDegree +1]].
	We assume that extremumPoints, poly and interpolationPoints has been allocated outside the function.
*/
void extremums(double *extremumPoints, const double *poly, double *interpolationPoints, fonction f, uint64_t polyDegree) {

	// The lower bound of interpolationPoints is the first extremum of the function
	// interpolation(x) - interpolated function(x) on the selected interval
	extremumPoints[0] = interpolationPoints[0];
	int extremumFound = 1;

	// Search for other extremes
    for (int i = 0; i < polyDegree + 1; i++) {

		// Search for an extremum on the intervals [interpolationPoints[i], interpolationPoints[i + 1]]
		double extremum = extremumLocal(poly, interpolationPoints[i], interpolationPoints[i + 1], f, polyDegree);

		// If the found extremum belongs to the chosen small interval, then it is kept. 
		// Otherwise it is that there was no extremum in this small interval.
        if (extremum >= interpolationPoints[i] && extremum <= interpolationPoints[i + 1]) {
			extremumPoints[extremumFound] = extremum;
			extremumFound++;
		}
    }

	// The upper bound of interpolationPoints is the first extremum of the function
	// interpolation(x) - interpolated function(x) on the selected interval
    extremumPoints[polyDegree + 1] = interpolationPoints[polyDegree + 1];
}


// Question 5
/*
	Create the matrix A (the matrix of the system).
	At the end of the procedure, A is a matrix of the form:
		1 x0    x0^2    ...  x0^n 1
		1 x1    x1^2    ...  x1^n -1
		 ...     ...    ...    ...
		1 xn+1  xn+1^2  ...  xn+1^n -1^n+1
	Allowing to solve the remez system allowing the interpolation of a function.
	We assume that A has been allocated outside the function.

*/
void createA(double *A, double *interpolationPoints, uint64_t polyDegree) {

	// Row Filling
	for (int i = 0; i < polyDegree + 2; i++) {
		double previousPuissance = 1;

		// Filling in the columns of the row you are on
		for (int j = 0; j < polyDegree + 2; j++) {

			// Filling the first element of a row with the value 1
			if (j == 0) {
				A[convert_ij_vers_ind(i, j, polyDegree + 2)] = previousPuissance;
			} else {
				
				// Filling the last element of a line with the value 1 (odd line) or -1 (even line)
				if (j == polyDegree + 1) {
					A[convert_ij_vers_ind(i, j, polyDegree + 2)] = -2 * (i % 2) + 1;
				
				// Filling of the rows with powers from 1 to polyDegree of the interpolation points
				} else {
					A[convert_ij_vers_ind(i, j, polyDegree + 2)] = previousPuissance * interpolationPoints[i];
					previousPuissance *= interpolationPoints[i];
				}
			}
		}
	}
}


// Question 6
/*
	Create the vector B (the right-hand vector of the system).
	At the end of the procedure, B is a vector containing the images of the 
	interpolation points by the interpolation function
	We assume that b and interpolationPoints has been allocated outside the function.
*/
void createB(double *b, double *interpolationPoints, fonction f, uint64_t n) {
	for (int i = 0; i < n + 2; i++) {
		b[i] = f(interpolationPoints[i]);
	}
}


// Question 8
/*
	Constructs the best polynomial approximation of degree polyDegree of a function f, 
	continuous over a bounded interval [boudA, boundB].
	After the procedure, a good interpolation of the function f was found.
	Return a double* variable:
		* A dynamic array (poly) containing the coefficients of the polynomial found.
		* This array can contain the extremes of the function interpolation(x) - interpolated function(x)
		* if isAuto is true.
*/
double* remez(double boundA, double boundB, fonction f, uint64_t polyDegree, bool isAuto) {
	
	// Create vectors and matrices
	double *poly = isAuto ? allocateVector(2 * polyDegree + 3) : allocateVector(polyDegree + 1);
	double *x = allocateVector(polyDegree + 2);
	double *b = allocateVector(polyDegree + 2);
	double *A = allocateMatrix(polyDegree + 2, polyDegree + 2);
	
	double *interpolationPoints = allocateVector(polyDegree + 2);
	double *extremumPoints = allocateVector(polyDegree + 2);
	double *interpolationErrors = allocateVector(polyDegree + 2);

	double maximumInterpolationError[2] = {100,200};

	// Create an array of n equidistant points
	createInterpolationPoints(interpolationPoints, boundA, boundB, polyDegree + 2);

	// Search for a polynomial interpolation until it is sufficiently precise, 
	// or until the precision does not increase any more
	for (int i = 0; fabs(maximumInterpolationError[0] - maximumInterpolationError[1]) > epsilon && maximumInterpolationError[0] > lambda; i++) {

		// Create the vector B
		createB(b, interpolationPoints, f, polyDegree);

		// Create the matrice A
		createA(A, interpolationPoints, polyDegree);

		// Solves the system Ax=B
		solveSystemGauss(x, A, b, polyDegree + 2);

		// Find all extremes of the function interpolation(x) - interpolated function(x)
		extremums(extremumPoints, x, interpolationPoints, f, polyDegree);
		
		// Find the values of all extremes of the function interpolation(x) - interpolated function(x)
		estimateError(interpolationErrors, x, f, boundA, boundB, polyDegree, polyDegree + 2, false);

		// Choosing the new interpolation points (extremes)
		for (int j = 0; j < polyDegree + 2; j++) {
			interpolationPoints[j] = extremumPoints[j];
		}

		// Find the maximum interpolation error
		maximumInterpolationError[i % 2] = getMax(interpolationErrors, polyDegree + 2);
		
	}

	// Create an array containing the coefficients of the polynomial found
	// If isAuto is true, this array also contains the extremes of the function so 
	// that the relative maximum error in the remez_plus function can easily be obtained.
	for (int k = 0; (k < polyDegree + 1) || (isAuto && k < 2 * polyDegree + 3); k++) {
		poly[k] = k < polyDegree + 1 ? x[k] : extremumPoints[k - polyDegree - 1];
	}

	// Release of dynamically allocated memory
	free(x);
	free(b);
	free(A);
	free(interpolationPoints);
	free(extremumPoints);
	free(interpolationErrors);

	return poly;
}


// Question 9
/*
	Constructs the best polynomial approximation of a function f, 
	continuous over a bounded interval [boudA, boundB].
	The polyDegree is searched by the procedure, so that the relative interpolation 
	error of the found polynomial is bounded by delta.
	After the procedure, a good interpolation of the function f was found.
	Return a double* variable:
		* A dynamic array (poly) containing the coefficients of the polynomial found.
*/
double* remezAuto(double boundA, double boundB, fonction f, double delta, int &polyDegree) {

	// Create vectors and the variable maximumInterpolationError, that will stop the 
	// algorithm if its value becomes less than delta
	double *x;
	double *interpolationErrors;
	double maximumInterpolationError = delta + 1;

	// Incrementation of the polyDegree as long as the relative interpolation error is not delta-bounded
	for (polyDegree = 1; maximumInterpolationError > delta ; polyDegree++) {

		// Search for an interpolation polynomial of the interpolated function
		x = remez(boundA, boundB, f, polyDegree, true);
		
		// Find the values of all extremes of the function interpolation(x) - interpolated function(x)
		interpolationErrors = allocateVector(polyDegree + 2);
		estimateError(interpolationErrors, x, f, boundA, boundB, polyDegree, polyDegree+2, true);

		// Find the maximum interpolation error
		maximumInterpolationError = getMax(interpolationErrors, polyDegree + 2);
		if (maximumInterpolationError < 0) { maximumInterpolationError = delta + 1; }

		// Release of dynamically allocated memory
		free(interpolationErrors);
		if (maximumInterpolationError > delta) free(x); 
	}

	// Create an array containing the coefficients of the polynomial found
	polyDegree--;
	double *poly = allocateVector(polyDegree + 1);
	for (int k = 0; k < polyDegree + 1; k++) {
		poly[k] = x[k];
	}
	
	// Release of dynamically allocated memory
	free(x);

	return poly;
}
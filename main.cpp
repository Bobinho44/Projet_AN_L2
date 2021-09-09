#include "remez.hpp"

using namespace std;

// ------- Declaration of constants -------

// Precision for the iteration of Newtion
double epsilon = 10e-12;

// Maximum number of iterations
int stop = 1000;

// Maximum precision of Remez
double lambda = 10e-3;
	
// ------- Declaration of functions to be interpolated -------
double f(double x) {return exp(x);} 

int main() {
	double *poly;
	uint64_t N = 100;
	
	int polyDegree = 4;
	double boundA = 0;
	double boundB = 5; 
	
	//poly = remezAuto(boundA, boundB, f, 0.001, polyDegree);
	poly = remez(boundA, boundB, f, polyDegree);
	
	cout << "\nDegree: " << polyDegree << endl;
	
	cout << "\nPolynome: " << endl;
	for(int i = 0; i < polyDegree + 1; ++i) {
		cout << "a" << i << " = " << poly[i] << endl;
	}

	double* interpolationErrorsAbs = allocateVector(N);
	double* interpolationErrorsRel = allocateVector(N);

	estimateError(interpolationErrorsAbs, poly, f, boundA, boundB, polyDegree, N, false);
	estimateError(interpolationErrorsRel, poly, f, boundA, boundB, polyDegree, N, true);
	
    double maxErrorAbs = getMax(interpolationErrorsAbs, N);
	double maxErrorRel = getMax(interpolationErrorsRel, N);

    cout << "\nMax absolute error : " << maxErrorAbs << endl;
    cout << "Max relative error : " << maxErrorRel << "\n" << endl;
	
	double *interpolationPoints = allocateVector(N);
	createInterpolationPoints(interpolationPoints, boundA, boundB, N);
	
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wwrite-strings"	
	plot(interpolationPoints, interpolationErrorsAbs, N, "Absolute approximation error for exp(x) with a degree-4 polynomial", "testPlotFuncAbsErrorDegree4[0,5]");
	plot(interpolationPoints, interpolationErrorsRel, N, "Relative approximation error for exp(x) with a degree-12 polynomial", "testPlotFuncRelErrorDegree4[0,5]");
	#pragma GCC diagnostic pop
	
	free(poly);
	free(interpolationErrorsAbs);
	free(interpolationErrorsRel);
	
}
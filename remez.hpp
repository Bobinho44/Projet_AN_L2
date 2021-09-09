#include <cassert>
#include <math.h>
#include <stdint.h>
#include "util.hpp"
#include "matrix.hpp"


/*
	Precision marker of the roots found by the Newton method.
*/
extern double epsilon;

/*
	Precision marker of the precision limit of the absolute error of interpolation of a 
	function by the remez method.
*/
extern double lambda;

/*
	Maximum number of iterations of the Newton method.
*/
extern int stop;


typedef double fonction(double); 


/*
	Create an array of nElements equidistant points in the interval [boundA, boundB].
	The array must contain at least two points (boudA and boundB).
	After the procedure, interpolationPoints is an array composed of the nElements points.
	We assume that interpolationPoints has been allocated outside the function.
*/
void createInterpolationPoints(double *interpolationPoints, double boundA, double boundB, uint64_t nElements);


/*
	Search for a root of the derivative of the function interpolation(x) - interpolated function(x) 
	on the [boundA, boundB] interval.
	We assume that poly has been allocated outside the function.
	Return a double:
		*The root of the derivative of the function interpolation(x) - interpolated function(x) 
		* on the [boundA, boundB] interval corresponds to an extremum of the function 
		* interpolation(x) - interpolated function(x) on that interval.
*/
double extremum_local(const double* poly, double boundA, double boundB, fonction f, uint64_t polyDegree);


/*
	At the end of the procedure, extremumPoints is a an array containing the extremes of 
	the function interpolation(x) - interpolated function(x) on the interval 
	[interpolationPoints[0], interpolationPoints[polyDegree +1]].
	We assume that extremumPoints, poly and interpolationPoints has been allocated outside the function.
*/
void extremums(double *extremumPoints, const double *poly, double *interpolationPoints, fonction f, uint64_t polyDegree);


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
void createA(double *A, double *interpolationPoints, uint64_t polyDegree);


/*
	Create the vector B (the right-hand vector of the system).
	At the end of the procedure, B is a vector containing the images of the 
	interpolation points by the interpolation function
	We assume that b and interpolationPoints has been allocated outside the function.
*/
void createB(double *b, double *interpolationPoints, fonction f, uint64_t n);


/*
	Constructs the best polynomial approximation of degree polyDegree of a function f, 
	continuous over a bounded interval [boudA, boundB].
	After the procedure, a good interpolation of the function f was found.
	Return a double* variable:
		* A dynamic array (poly) containing the coefficients of the polynomial found.
		* This array can contain the extremes of the function interpolation(x) - interpolated function(x)
		* if isAuto is true.
*/
double* remez(double boundA, double boundB, fonction f, uint64_t polyDegree, bool isAuto = false);


/*
	Constructs the best polynomial approximation of a function f, 
	continuous over a bounded interval [boudA, boundB].
	The polyDegree is searched by the procedure, so that the relative interpolation 
	error of the found polynomial is bounded by delta.
	After the procedure, a good interpolation of the function f was found.
	Return a double* variable:
		* A dynamic array (poly) containing the coefficients of the polynomial found.
*/
double* remezAuto(double boundA, double boundB, fonction f, double delta, int &polyDegree);
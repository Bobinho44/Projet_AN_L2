#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <iostream>


/* 
	This funciton evaluates a polynomial of degree polyDegree using Horner's method at the value x.
	Attention, polyDegree is the degree of the polynomial, i.e. poly has polyDegree+1 elements.
*/
double horner(const double* poly, double x, uint64_t polyDegree);


/* 
	Given a polynomial poly of degree polyDegree, described as an array of double-precision coefficients, 
	a pointer to a function f, and an interval [a;b], this function fills out the array errors
	with the relative errors |poly(x_i) - f(x_i)/func(x_i)| or the absolute error |poly(x_i) - f(x_i)|
	for each sample point x_i in [a;b], where the number of element is given as a parameter.
	The array errors is assumed to be allocated outside the function. 
	This is a naive and rough estimation of the approximation error.
*/
void estimateError(double *interpolationError, const double* poly, double f(double), double boundA, double boundB, uint64_t polyDegree, uint64_t nElements, bool isRel);


/* 
	Given a n-element array of doubles, this function computes the maximum element.
 */
double getMax(const double* array, uint64_t n);


/* 
	This function plots a function defined by N 2D points using the gnuplot software (http://www.gnuplot.info)
	Inputs : an array of points along x-axis, an array of points along y-axis, the number of points (N<= the length og arrays), the title of the plot and a filename.
*/
void plot(double *xvals, double *yvals, uint64_t N, char* title, char* filename);



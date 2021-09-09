#include "util.hpp"
#include <math.h>
#include <cassert>

using namespace std;


/* 
	This funciton evaluates a polynomial of degree polyDegree using Horner's method at the value x.
	Attention, polyDegree is the degree of the polynomial, i.e. poly has polyDegree+1 elements.
*/
double horner(const double* poly, double x, uint64_t polyDegree) { 
    double result = poly[polyDegree];  // Initialize result 
    //cout << "horner" << endl;
    //cout << result << endl;
    // Evaluate value of polynomial using Horner's method 
    for (int i=polyDegree-1; i>=0; i--) {
        result = result*x + poly[i]; 
        //cout << result << endl;
    }
    //cout << result << endl;
    return result; 
} 


/* 
	Given a polynomial poly of degree polyDegree, described as an array of double-precision coefficients, 
	a pointer to a function f, and an interval [a;b], this function fills out the array errors
	with the relative errors |poly(x_i) - f(x_i)/func(x_i)| or the absolute error |poly(x_i) - f(x_i)|
	for each sample point x_i in [a;b], where the number of element is given as a parameter.
	The array errors is assumed to be allocated outside the function. 
	This is a naive and rough estimation of the approximation error.
*/
void estimateError(double *interpolationError, const double* poly, double f(double), double boundA, double boundB, uint64_t polyDegree, uint64_t nElements, bool isRel) {

	// The array must contain at least two points (boudA and boundB)
    assert(nElements>1);

	// Interval between two points
	double interval = (boundB - boundA) / (nElements - 1);

	// Create the array containing nElements interpolation point errors
    for (int i = 0; i < nElements; i++) {
		double x = boundA + i * interval;
        interpolationError[i] = isRel ? fabs((horner(poly, x, polyDegree) - f(x)) / f(x)) : fabs(horner(poly, x, polyDegree) - f(x));
		if (isinf(interpolationError[i]) || isnan(interpolationError[i])) {
			interpolationError[i] = -1;
		}
    }
}


/* 
	Given a n-element array of doubles, this function computes the maximum element.
 */
double getMax(const double* array, uint64_t n){
    double max =  array[0];
    for(int i = 1; i < n; ++i){
        if(array[i]>max)    max = array[i];
    }
    return max;
}


/* 
	This function plots a function defined by N 2D points using the gnuplot software (http://www.gnuplot.info)
	Inputs : an array of points along x-axis, an array of points along y-axis, the number of points (N<= the length og arrays), the title of the plot and a filename.
*/
void plot(double *xvals, double *yvals, uint64_t N, char* title, char* filename)
{
    char gnuplot_command1[1000] = {0};
    sprintf(gnuplot_command1, "set terminal postscript eps enhanced color solid colortext 9 font \",15\" \n  set output '%s.eps' \n set key left below \n set style data line \n plot \"%s\" title \"%s\" \n", filename, filename, title);

    FILE * temp = fopen(filename, "w");

    

    for (int i=0; i < N; i++)
    {
        fprintf(temp, "%le %le \n", xvals[i], yvals[i]); //Write the data to a temporary file
    }
    fclose(temp);

    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "%s", gnuplot_command1);
  //  fclose(gnuplotPipe);
}



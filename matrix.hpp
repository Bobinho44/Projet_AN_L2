#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


/* 
	Memory allocation for a matrix of size n x m and initilization to 0  
*/
double* allocateMatrix(uint64_t n,uint64_t m);


/* 
	Frees the memory allocated to matrix A
*/
double* allocateVector(uint64_t n);


/* 
	Allocates a n sized vector and initializes all entries to 0 
*/
void freeMatrix(double *A);


/* 
	Trees the memory allocated to a vector
*/
void freeVector (double *v);


/* 
	Sets a n * m matrix A to all zeros 
*/
void setMatrixZero(double *A, uint64_t n, uint64_t m);


/* 
	Sets a n * n matrix A to identity 
*/
void setMatrixIdentity(double *A, uint64_t n);


/* 
	Copies a matrix  
*/
void copyMatrix(double *B, double *A, uint64_t n, uint64_t m);


/*
	Writes a matrix to a stream. For example, writing a matrix to standard output is
	writeMatrix(stdout, A, n, m);
	A sream can also be a file. 
*/
void writeMatrix(FILE *stream, double *A, uint64_t n, uint64_t m);

/*
	The function computes the element-by-element abs of matrix A
*/
void absMatrix(double *Aabs,double *A, uint64_t n, uint64_t m);


/* 
	For a double m x n matrix A the function returns its maximum in absolute value
	element. 
*/
double getMaxInMatrix(double max, double *A, uint64_t n, uint64_t m);


/*
	Performs subtraction of two matrix A (size n x m) and B (size n x m).
	The result S = A - B is a n x m matrix.
	We consider that S is allocated outside the function.
*/
void matrixSub(double *S, double *A, double *B, uint64_t n, uint64_t m);


/*
	Performs addition of two matrix A (size n x m) and B (size n x m).
	The result S = A + B is a n x m matrix.
	We consider that S is allocated outside the function.
*/
void matrixAdd(double *S, double *A, double *B, uint64_t n, uint64_t m);


/*
	Convert the position of an element in a matrix into an indice for storage in an array.
*/
int convert_ij_vers_ind(int i, int j, int k);


/* 
	Performs naive multiplication of matrix A (size p x k) by a matrix B (size k x r).
	The result matrix S = A*B  is of size (k x r).
	We assume that S has already been allocated outside the function.
*/
void matrixMultiplyNaive(double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r);


/* 
	Performs a multiplication of two sqaure matrices A and B (size n x n) by Strassen algorithm.
    We assume that S has already been allocated outside the function.
*/
void matrixMultiplyStrassen(double *S, double *A, double *B, uint64_t n);


/* 
    Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
    Uses iterative ascension algorithm. 
    After the procedure, x contains the solution of Ax=b.
    We assume that x has been allocated outside the function.
*/
void solveTriangularSystemUP(double *x, double *A, double *b, uint64_t n);


/* 
    Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/
bool triangularize(double *A, double *b, uint64_t n);


/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool solveSystemGauss(double *x, double *A, double *b, uint64_t n);

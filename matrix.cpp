#include "matrix.hpp"


/*
    Nous allons stocker des matrices dans des tablaux à 1 dimension:
    Par exemple, matrice de taille n x m est stockué comme:
        A_tab = [A_11, A_12, ... A_1m, A_21, ... A_2m, ... , A_1n, ..., A_nm]

    Par conséquant, un élement A_ij aurait quelle indice dans le tableau A_tab ? 

 */


/* 
	Memory allocation for a matrix of size n x m and initilization to 0  
*/
double *allocateMatrix(uint64_t n,uint64_t m) {
  double *A;
  A = (double *) calloc (n * m, sizeof(double));
  return A;
}


/* 
	Frees the memory allocated to matrix A
*/
void freeMatrix(double *A) {
    free(A);
}


/* 
	Allocates a n sized vector and initializes all entries to 0 
*/
double *allocateVector(uint64_t n) {
  double *v; 
  v = (double *) calloc(n, sizeof(double));
  return v;
}


/* 
	Trees the memory allocated to a vector
*/
void freeVector(double *v) {
  free(v);
}


/* 
	Sets a n * m matrix A to all zeros 
*/
void setMatrixZero(double *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
        /* Note that for a n x m matrix flattened to a 1D array, 
        element A_ij has index i * m + j
        */
      A[i * m + j] = 0.0; 
    }
  }
}


/* 
	Sets a n * n matrix A to identity 
*/
void setMatrixIdentity(double *A, uint64_t n) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
     A[i * n + j] = 0.0;
    }
    A[i * n + i] = 1.0;
  }
}


/* 
	Copies a matrix  
*/
void copyMatrix(double *B, double *A, uint64_t n, uint64_t m) {
  uint64_t i,j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      B[i * m + j] = A[i * m + j]; 
    }
  }
}


/*
	Writes a matrix to a stream. For example, writing a matrix to standard output is
	writeMatrix(stdout, A, n, m);
	A sream can also be a file. 
*/
void writeMatrix(FILE *stream, double *A, uint64_t n, uint64_t m) {
	fprintf(stream, "%d %d \n", (int)n, (int)m);
	int i, j;
	for(i = 0; i < n; ++i)
	{
	      for(j = 0; j < m; ++j)
	      {
		      fprintf(stream, "%f \t", A[i * m + j]);
	      }
	      fprintf(stream, "\n");
	}
}


/*
	The function computes the element-by-element abs of matrix A
*/
void absMatrix(double *Aabs,double *A, uint64_t n, uint64_t m) {
	uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            Aabs[i*m + j] = fabs(A[i*m + j]);
		}
	}

}


/* 
	For a double m x n matrix A the function returns its maximum in absolute value
	element. 
*/
double getMaxInMatrix(double max, double *A, uint64_t n, uint64_t m) {
	double maxA = fabs(A[0]);
	double current = fabs(A[0]);
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			current = fabs(A[i * m + j]);
			if(current > maxA)
				maxA = current;
		}
	}
    return maxA;
}


/*
	Performs addition of two matrix A (size n x m) and B (size n x m).
	The result S = A + B is a n x m matrix.
	We consider that S is allocated outside the function.
*/
void matrixAdd(double *S, double *A, double *B, uint64_t n, uint64_t m) {
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] + B[i*m + j];
		}
	}
}


/*
	Performs subtraction of two matrix A (size n x m) and B (size n x m).
	The result S = A - B is a n x m matrix.
	We consider that S is allocated outside the function.
*/
void matrixSub(double *S, double *A, double *B, uint64_t n, uint64_t m) {
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] - B[i*m + j];
		}
	}
}


/*
	Convert the position of an element in a matrix into an indice for storage in an array.
*/
int convert_ij_vers_ind(int i, int j, int k) {
	return i * k + j;
}


/* 
	Performs naive multiplication of matrix A (size p x k) by a matrix B (size k x r).
	The result matrix S = A*B  is of size (k x r).
	We assume that S has already been allocated outside the function.
*/
void matrixMultiplyNaive(double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r) {
	
	// Browse rows of the matrix A
      for (int h = 0; h < p; h++) {
		  
		// Browse columns of marice B
        for (int i = 0; i < r; i++) {
          int ind = convert_ij_vers_ind(h, i, r);
          S[ind] = 0;
		  
		  // Browse columns of the matrix A
          for (int j = 0; j < k; j++) {
			  
			// Intermediate product (item of column i of A * item of row i of B)
            S[ind] += A[convert_ij_vers_ind(h, j, k)]*B[convert_ij_vers_ind(j, i, r)];
          }
        }
      }
}


/* 
	Performs a multiplication of two sqaure matrices A and B (size n x n) by Strassen algorithm.
    We assume that C has already been allocated outside the function.
*/
void matrixMultiplyStrassen (double *C, double *A, double *B, uint64_t n) {
	
	// If n is not a multiple of 2
	if (!((n&(n-1)) == 0) && (n!=0)) {

		// Give the closest multiple of 2 greater than n
		int newN = pow(2, int(ceil(log2(n))));

		double *newC = allocateMatrix(newN, newN);
		double *newA = allocateMatrix(newN, newN);
		double *newB = allocateMatrix(newN, newN);

		// Rewrite the values of A, B and C in the matrices newA, newB, newC 
		// (square matrix multiple of 2)
		for (int i = 0; i < n; i++) {
      		for (int j = 0; j < n; j++) {
				newA[convert_ij_vers_ind(i, j, n) + i*1] = A[convert_ij_vers_ind(i, j, n)];
        		newB[convert_ij_vers_ind(i, j, n) + i*1] = B[convert_ij_vers_ind(i, j, n)];
    		}
		}

		// Use the algorithm on matrices of good size
		matrixMultiplyStrassen(newC, newA, newB, newN);

		// Rewrite the values of newC in the matrices C
		for (int i = 0; i < n; i++) {
      		for (int j = 0; j < n; j++) {
				C[convert_ij_vers_ind(i, j, n)] = newC[convert_ij_vers_ind(i, j, newN)];
    		}
		}
	}

	// Condition for stopping the recurrence
	else if (n == 1) {
		C[0] = A[0] * B[0];
	}

	// Recurrence
	else {

		// Partition A, B and C
		// Declaration of sum  matrices and product matrices
		double *C11 = allocateMatrix(n/2, n/2);
		double *C12 = allocateMatrix(n/2, n/2);
		double *C21 = allocateMatrix(n/2, n/2);
		double *C22 = allocateMatrix(n/2, n/2);

		double *A11 = allocateMatrix(n/2, n/2);
		double *A12 = allocateMatrix(n/2, n/2);
		double *A21 = allocateMatrix(n/2, n/2);
		double *A22 = allocateMatrix(n/2, n/2);

		double *B11 = allocateMatrix(n/2, n/2);
		double *B12 = allocateMatrix(n/2, n/2);
		double *B21 = allocateMatrix(n/2, n/2);
		double *B22 = allocateMatrix(n/2, n/2);

		// Switch from C to [C11, C12, C21, C22]
		// Switch from A to [A11, A12, A21, A22]
		// Switch from B to [B11, B12, B21, B22]
		int Q1 = 0, Q2 = 0, Q3 = 0, Q4 = 0;
		for (int i = 0; i < n; i++) {
  			for (int j = 0; j < n; j++) {
				// Q1
		 		if (i < n/2 && j < n/2) {
					C11[Q1] = C[convert_ij_vers_ind(i, j, n)];
					A11[Q1] = A[convert_ij_vers_ind(i, j, n)];
					B11[Q1] = B[convert_ij_vers_ind(i, j, n)];
			  		Q1++;
				 }

				// Q2
		 		if (i < n/2 && j >= n/2) {
					C12[Q2] = C[convert_ij_vers_ind(i, j, n)];
					A12[Q2] = A[convert_ij_vers_ind(i, j, n)];
					B12[Q2] = B[convert_ij_vers_ind(i, j, n)];
			  		Q2++;
				 }

				// Q3
		 		if (i >= n/2 && j < n/2) {
					C21[Q3] = C[convert_ij_vers_ind(i, j, n)];
					A21[Q3] = A[convert_ij_vers_ind(i, j, n)];
					B21[Q3] = B[convert_ij_vers_ind(i, j, n)];
			  		Q3++;
				 }

				// Q4
		 		if (i >= n/2 && j >= n/2) {
					C22[Q4] = C[convert_ij_vers_ind(i, j, n)];
					A22[Q4] = A[convert_ij_vers_ind(i, j, n)];
					B22[Q4] = B[convert_ij_vers_ind(i, j, n)];
			  		Q4++;
				 }
			  }
		}

		double *S1 = allocateMatrix(n/2, n/2);
		double *S2 = allocateMatrix(n/2, n/2);
		double *S3 = allocateMatrix(n/2, n/2);
		double *S4 = allocateMatrix(n/2, n/2);
		double *S5 = allocateMatrix(n/2, n/2);
		double *S6 = allocateMatrix(n/2, n/2);
		double *S7 = allocateMatrix(n/2, n/2);
		double *S8 = allocateMatrix(n/2, n/2);
		double *S9 = allocateMatrix(n/2, n/2);
		double *S10 = allocateMatrix(n/2, n/2);

		double *P1 = allocateMatrix(n/2, n/2);
		double *P2 = allocateMatrix(n/2, n/2);
		double *P3 = allocateMatrix(n/2, n/2);
		double *P4 = allocateMatrix(n/2, n/2);
		double *P5 = allocateMatrix(n/2, n/2);
		double *P6 = allocateMatrix(n/2, n/2);
		double *P7 = allocateMatrix(n/2, n/2);

		//  calculate  the  sum  matrices 
		matrixSub(S1, B12, B22, n/2, n/2);
		matrixAdd(S2, A11, A12, n/2, n/2);
		matrixAdd(S3, A21, A22, n/2, n/2);
		matrixSub(S4, B21, B11, n/2, n/2);
		matrixAdd(S5, A11, A22, n/2, n/2);
		matrixAdd(S6, B11, B22, n/2, n/2);
		matrixSub(S7, A12, A22, n/2, n/2);
		matrixAdd(S8, B21, B22, n/2, n/2);
		matrixSub(S9, A11, A21, n/2, n/2);
		matrixAdd(S10, B11, B12, n/2, n/2);

		//  calculate  the  product  matrices 
		matrixMultiplyStrassen(P1, A11, S1, n/2);
		matrixMultiplyStrassen(P2, S2, B22, n/2);
		matrixMultiplyStrassen(P3, S3, B11, n/2);
		matrixMultiplyStrassen(P4, A22, S4, n/2);
		matrixMultiplyStrassen(P5, S5, S6, n/2);
		matrixMultiplyStrassen(P6, S7, S8, n/2);
		matrixMultiplyStrassen(P7, S9, S10, n/2);

		//  calculate  the  final  product  sub  matrices 
		double *P45 = allocateMatrix(n/2, n/2);
		double *P62 = allocateMatrix(n/2, n/2);
		matrixAdd(P45, P4, P5, n/2, n/2);
		matrixSub(P62, P6, P2, n/2, n/2);

		double *P15 = allocateMatrix(n/2, n/2);
		double *P153 = allocateMatrix(n/2, n/2);
		matrixAdd(P15, P1, P5, n/2, n/2);
		matrixSub(P153, P15, P3, n/2, n/2);

		matrixAdd(C11, P45, P62, n/2, n/2);
		matrixAdd(C12, P1, P2, n/2, n/2);
		matrixAdd(C21, P3, P4, n/2, n/2);
		matrixSub(C22, P153, P7, n/2, n/2);
		
		// Switch from [C11, C12, C21, C22] to C
		int nbQ = 0;
		for (int i = 0; i < n/2; i++) {
			for (int j = 0; j < n/2; j++) {
				C[convert_ij_vers_ind(i, j, n)] = C11[nbQ];
				C[convert_ij_vers_ind(i, j + n/2, n)] = C12[nbQ];
				C[convert_ij_vers_ind(i + n/2, j, n)] = C21[nbQ];
				C[convert_ij_vers_ind(i + n/2, j + n/2, n)] = C22[nbQ];
				nbQ++;
			}
		}
	}
}


/* 
    Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
    Uses iterative ascension algorithm. 
    After the procedure, x contains the solution of Ax=b.
    We assume that x has been allocated outside the function.
*/
void solveTriangularSystemUP(double *x, double *A, double *b, uint64_t n) {
	x[n-1] = b[n-1] / A[convert_ij_vers_ind(n-1, n-1, n)];
	
	// Browse rows of the matrix 
	for (int i = n-2; i >= 0; i--) {
		double sum = 0;
		
		// Browse column of the matrix
		for (int j = i+1; j <= n-1; j++) {
			sum += A[convert_ij_vers_ind(i, j, n)] * x[j];
		}
		
		// Find the vector
		x[i] = (b[i] - sum)/A[convert_ij_vers_ind(i, i, n)];
	}
}


/* 
    Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/
bool triangularize(double *A, double *b, uint64_t n) {
    for (int i = 0; i < n; i++) {

		// Line reversal in case of a null pivot
        if (A[convert_ij_vers_ind(i, i, n)] == 0) {
			bool findNewLine = false;
			for (int w = i+1; w < n && !findNewLine; w++) {
				if (A[convert_ij_vers_ind(i, w, n)] != 0) {
					double tempA;
					double tempB;
					for (int x = 0; x < n; x++) {
						tempA = A[convert_ij_vers_ind(i, x, n)];
						A[convert_ij_vers_ind(i, x, n)] = A[convert_ij_vers_ind(w, x, n)];
						A[convert_ij_vers_ind(w, x, n)] = tempA;
					}
					tempB = b[i];
					b[i] = b[w];
					b[w] = tempB;
					findNewLine = true;
				}
			}
			if (!findNewLine) return false;
		}

		// Triangularization
        for (int j = 1+i; j < n; j++) {
          double facteur = (A[convert_ij_vers_ind(j, i, n)]/A[convert_ij_vers_ind(i, i, n)]);
          b[j] -= facteur*b[i];
          for (int k = 0; k < n; k++) {
            A[convert_ij_vers_ind(j, k, n)] -= facteur*A[convert_ij_vers_ind(i, k, n)];
          }
        }
    }
    return true;
}


// Question 7
/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool solveSystemGauss(double *x, double *A, double *b, uint64_t n) {
    if (!triangularize(A, b, n)) return false;
    solveTriangularSystemUP(x, A, b, n);
    return true;
}


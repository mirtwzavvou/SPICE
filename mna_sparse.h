#include "csparse.h"


#ifndef MNASPARSE_H
#define MNASPARSE_H

cs *A_sparse;			// Pinakas A (mna) se morfi triplet
cs *C_sparse;			// Pinakas A (mna) se morfi compressed-column
cs *D_sparse;			// Pinakas C (mna) se morfi triplet
cs *E_sparse;			// Pinakas C (mna) se morfi triplet

css *S;				// Pinakas S gia tin LU kai Cholesky
csn *N;				// Pinakas N gia tin LU kai Cholesky

int sizeA_sparse;		//non zeros stoixeia tou pinaka A
int sizeD_sparse;		//non zeros stoixeia tou pinaka D
double *B_sparse;		//pinakas deksiou melous eksiswshs B ( double, mege8ous : [(n-1)+m2] x 1 ) -> grammi!	

double *x_sparse;

void computeMNASparse();
void solveSparse(double time);
void solve_spdSparse(double time);
void conjugate_gradient_sparse(cs *A, double *b, int n, double *x, double itol);
void bi_conjugate_gradient_sparse(cs *A, double *b, int n, double *x, double itol);
#endif
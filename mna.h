
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

#ifndef MNA_H
#define MNA_H

gsl_matrix *A;			// 2D pinakas aristerou melous eksiswshs A ( double, mege8ous : [(n-1)+m2] x [(n-1)+m2] )
gsl_matrix *C;			// 2D pinakas C aristerou melous eksiswshs  ( double, mege8ous : [(n-1)+m2] x [(n-1)+m2] )
//double *temp;			//ka8e grammi tou pinaka				
int sizeA;				//[(n-1)+m2]x[(n-1)+m2]

gsl_vector *B;			//pinakas deksiou melous eksiswshs B ( double, mege8ous : [(n-1)+m2] x 1 ) -> grammi!	
int sizeB;				//[(n-1)+m2] x 1

gsl_vector *x;

gsl_permutation * p;		//dianusma meta8esewn

void computeMNA();
void solve(double time);
void solve_spd(double time);
void conjugate_gradient(gsl_matrix *a,gsl_vector *b,gsl_vector *X,int n,double tolerance);
void bi_conjugate_gradient(gsl_matrix *a,gsl_vector *b,gsl_vector *X,int n,double tolerance);

void preconditioner_diag(gsl_vector *preconditioner,gsl_matrix *matrix);

void solveEquation(gsl_vector *preconditioner,gsl_vector *right,gsl_vector *left);

double *gslvector2double(gsl_vector *V, int size);
#endif

#include "csparse.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

#ifndef TRANSIENT_H
#define TRANSIENT_H


void transient();

void TRAN_backward_euler_sparse(double *x_sparse, cs *E_sparse,double *e);
void TRAN_backward_euler(gsl_vector *x, gsl_matrix *C,gsl_vector *e);
void TRAN_trapezoidial_sparse(double *x_sparse, cs *E_sparse,double *e,double *e_prev);
void TRAN_trapezoidial(gsl_vector *x, gsl_matrix *C,gsl_vector *e,gsl_vector *e_prev);

double linear_value(double x1, double y1, double x2, double y2, double x);
gsl_vector *compute_E(double t);

#endif

#ifndef EIGENVALUES_H
#define EIGENVALUES_H 1

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "parameters.h"

gsl_vector* E;  /* eigenvalues */
gsl_matrix_complex* EV; /* eigenvectors */
gsl_eigen_hermv_workspace* w; /* workspace for eigen routines */
gsl_vector_complex* An; /* A|n> */

void initialize_eigen ()
{
  E = gsl_vector_alloc (8);
  EV = gsl_matrix_complex_alloc (8, 8);
  w = gsl_eigen_hermv_alloc (8);
  An = gsl_vector_complex_alloc (8);
}

void get_eigen (gsl_matrix_complex* A)
{
  gsl_eigen_hermv (A, E, EV, w);
}

void eigen_free ()
{
  gsl_eigen_hermv_free (w);
  gsl_vector_free (E);
  gsl_matrix_complex_free (EV);
  gsl_vector_complex_free (An);
}

double bracket (gsl_matrix_complex* A, int n)
/* this function return <n|A|n>,
   here A should be hermitian. */
{
  gsl_complex nAn; /* <n|A|n> */
  gsl_vector_complex_view n_ket = gsl_matrix_complex_column (EV, n);
  
  gsl_blas_zgemv (CblasNoTrans, complex_1, 
		  A, &n_ket.vector, complex_0, An); 
  gsl_blas_zdotc (&n_ket.vector, An, &nAn);

  return GSL_REAL (nAn);
}

#endif /* EIGENVALUES_H */

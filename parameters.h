#ifndef PARAMETERS_H
#define PARAMETERS_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>
#include <math.h>

const double J1 = -0.037;
const double J2 = 0.069;
const double c2 = -0.033 * 0;
const double c4 = 0.003 * 0;
const double g = 2;
const double u = 0.6717139; /* this is mu_B/k_B */

/* if T is not initialized in main function, 
   the error of dividing 0 should occur. */
double T = 0; 
double B[3] = {0}; /* Bx, By, Bz */

double sa[3]; /* <sax>, <say>, <saz> */
double sb[3];
double Ha[3]; /* Hax, Hay, Haz */
double Hb[3];

gsl_complex complex_0;
gsl_complex complex_1;
gsl_complex complex_i;

gsl_matrix_complex* Sx; 
gsl_matrix_complex* Sy;
gsl_matrix_complex* Sz;
gsl_matrix_complex* Sx2; /* sx^2 */
gsl_matrix_complex* Sy2; /* sy^2 */
gsl_matrix_complex* Sx4; /* sx^4 */
gsl_matrix_complex* Sy4; /* sy^4 */
gsl_matrix_complex* S2; /* sx^2 + sy^2 */
gsl_matrix_complex* S4; /* sx^4 + sy^4 */


void initialize_matrice ()
{
  complex_0 = gsl_complex_rect (0, 0);
  complex_1 = gsl_complex_rect (1, 0);
  complex_i = gsl_complex_rect (0, 1);

  Sx = gsl_matrix_complex_calloc (8, 8);
  Sy = gsl_matrix_complex_calloc (8, 8);
  Sz = gsl_matrix_complex_calloc (8, 8);
  S2 = gsl_matrix_complex_calloc (8, 8);
  S4 = gsl_matrix_complex_calloc (8, 8);
  Sx2 = gsl_matrix_complex_calloc (8, 8);
  Sy2 = gsl_matrix_complex_calloc (8, 8);
  Sx4 = gsl_matrix_complex_calloc (8, 8);
  Sy4 = gsl_matrix_complex_calloc (8, 8);

  /* set Sx */
  for (int i = 0; i < 7; ++i)
    {
      double sigma = 3.5 - i;
      gsl_matrix_complex_set 
	(Sx, i, i+1, gsl_complex_rect 
	 (0.5*sqrt((3.5+sigma)*(3.5-sigma+1)), 0));
      gsl_matrix_complex_set 
	(Sx, i+1, i, gsl_matrix_complex_get (Sx, i, i+1));
    }
	   
  /* set Sy */
  for (int i = 0; i < 7; ++i)
    {
      double sigma = 3.5 - i;
      gsl_matrix_complex_set 
	(Sy, i, i+1, gsl_complex_rect 
	 (0, -0.5*sqrt((3.5+sigma)*(3.5-sigma+1))));
      gsl_matrix_complex_set
	(Sy, i+1, i, gsl_complex_conjugate 
	 (gsl_matrix_complex_get (Sy, i, i+1)));
    }
      
  /* set Sz */
  for (int i = 0; i < 8; ++i)
    gsl_matrix_complex_set 
      (Sz, i, i, gsl_complex_rect (3.5-i, 0));

  /* set Sx2 = Sx * Sx */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
		  complex_1, Sx, Sx, complex_0, Sx2);
  
  /* set Sx4 = Sx2 * Sx2 */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
		  complex_1, Sx2, Sx2, complex_0, Sx4);

  /* set Sy2 = Sy * Sy */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
		  complex_1, Sy, Sy, complex_0, Sy2);
  
  /* set Sy4 = Sy2 * Sy2 */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
		  complex_1, Sy2, Sy2, complex_0, Sy4);
  
  /* set S2 */
  gsl_matrix_complex_memcpy (S2, Sx2);
  gsl_matrix_complex_add (S2, Sy2);

  /* set S4 */
  gsl_matrix_complex_memcpy (S4, Sx4);
  gsl_matrix_complex_add (S4, Sy4);
}

void matrice_free ()
{
  gsl_matrix_complex_free (Sx);
  gsl_matrix_complex_free (Sy);
  gsl_matrix_complex_free (Sz);
  gsl_matrix_complex_free (Sx2);
  gsl_matrix_complex_free (Sx4);
  gsl_matrix_complex_free (Sy2);
  gsl_matrix_complex_free (Sy4);
  gsl_matrix_complex_free (S2);
  gsl_matrix_complex_free (S4);
}

void print_matrix_complex (gsl_matrix_complex* A)
{
  for (int i = 0; i < A->size1; ++i)
    {
      for (int j = 0; j < A->size2; ++j)
	printf ("%g+%gi\t\t", 
		GSL_REAL (gsl_matrix_complex_get (A, i, j)),
		GSL_IMAG (gsl_matrix_complex_get (A, i, j)));
      printf ("\n");
    }
}

#endif /* PARAMETERS_H */

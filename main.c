#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "parameters.h"
#include "eigenvalues.h"

double exp_aEn[8]; /* exp (-E_n/T) for site a */
double exp_bEn[8]; /* exp (-E_n/T) for site b */

void initialize ()
{
  initialize_matrice ();
  initialize_eigen ();
}

void clean ()
{
  matrice_free ();
  eigen_free ();
}

double a_partition ()
{
  double sum = 0;
  for (int i = 0; i < 8; ++i)
    sum += exp_aEn[i];
  return sum;
}

double b_partition ()
{
  double sum = 0;
  for (int i = 0; i < 8; ++i)
    sum += exp_bEn[i];
  return sum;
}

void iteration ()
{
  gsl_matrix_complex* H =  /* Hamiltonian */
    gsl_matrix_complex_alloc (8, 8);
  gsl_matrix_complex* tmp = gsl_matrix_complex_alloc (8,8);

  sa[0] = sa[1] = -0.2; sa[2] = 0.5; /* <sax>, <say>, <saz> */
  sb[0] = sb[1] = 0.2; sb[2] = -0.54;
  double sa_new[3], sb_new[3];
  double z, diff;

  do 
    {
      for (int i = 0; i < 3; ++i)
	{
	  Ha[i] = 6*J1*sb[i] + 12*J2*sa[i] + g*u*B[i];
	  Hb[i] = 6*J1*sa[i] + 12*J2*sb[i] + g*u*B[i];
	}
      /* Hamiltonian for a type */
      gsl_matrix_complex_memcpy (H, Sx);
      gsl_matrix_complex_scale (H, gsl_complex_rect (Ha[0], 0));
      gsl_matrix_complex_memcpy (tmp, Sy);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (Ha[1], 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_memcpy (tmp, Sz);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (Ha[2], 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_memcpy (tmp, S2);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (c2, 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_memcpy (tmp, S4);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (c4, 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_scale (H, gsl_complex_rect (-1, 0));

      get_eigen (H);
      for (int i = 0; i < 8; ++i)
	exp_aEn[i] = exp (-gsl_vector_get (E,i) / T);
      z = a_partition ();
      
      sa_new[0] = sa_new[1] = sa_new[2] = 0;
      for (int i = 0; i < 8; ++i)
	{
	  sa_new[0] += exp_aEn[i] * bracket (Sx, i);
	  sa_new[1] += exp_aEn[i] * bracket (Sy, i);
	  sa_new[2] += exp_aEn[i] * bracket (Sz, i);
	}
      sa_new[0] /= z; sa_new[1] /= z; sa_new[2] /= z;

      /* Hamiltonian for b type */
      gsl_matrix_complex_memcpy (H, Sx);
      gsl_matrix_complex_scale (H, gsl_complex_rect (Hb[0], 0));
      gsl_matrix_complex_memcpy (tmp, Sy);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (Hb[1], 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_memcpy (tmp, Sz);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (Hb[2], 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_memcpy (tmp, S2);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (c2, 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_memcpy (tmp, S4);
      gsl_matrix_complex_scale (tmp, gsl_complex_rect (c4, 0));
      gsl_matrix_complex_add (H, tmp);
      gsl_matrix_complex_scale (H, gsl_complex_rect (-1, 0));

      get_eigen (H);
      for (int i = 0; i < 8; ++i)
	exp_bEn[i] = exp (-gsl_vector_get (E,i) / T);
      z = b_partition ();
      
      sb_new[0] = sb_new[1] = sb_new[2] = 0;
      for (int i = 0; i < 8; ++i)
	{
	  sb_new[0] += exp_bEn[i] * bracket (Sx, i);
	  sb_new[1] += exp_bEn[i] * bracket (Sy, i);
	  sb_new[2] += exp_bEn[i] * bracket (Sz, i);
	}
      sb_new[0] /= z; sb_new[1] /= z; sb_new[2] /= z;
      
      diff = 0;
      for (int i = 0; i < 3; ++i)
	{
	  diff += fabs (sa[i] - sa_new[i]);
	  diff += fabs (sb[i] - sb_new[i]);
	  sa[i] = sa_new[i];
	  sb[i] = sb_new[i];
	}
    } while (diff > 1e-6);   
  gsl_matrix_complex_free (H);
  gsl_matrix_complex_free (tmp);
}
    
double magnet ()
{
  return sqrt 
    ((sa[0]+sb[0])*(sa[0]+sb[0])+
     (sa[1]+sb[1])*(sa[1]+sb[1])+
     (sa[2]+sb[2])*(sa[2]+sb[2]));
}
      
int main (int argc, char* argv[])
{
  initialize ();
  B[0] = B[1] = 0;
  B[2] = 0;

  T = atof (argv[1]);

  printf ("#B\tM/u_B\n");
  for (B[2] = 0; B[2] <= 7; B[2] += 0.1)
    {
      iteration ();
      printf ("%g\t%g\n", B[2], magnet ());
    }
  clean ();
  return 0;     
}

/*
**      File:   nrutil.h
**      Purpose: Memory allocation routines borrowed from the
**              book "Numerical Recipes" by Press, Flannery, Teukolsky,
**              and Vetterling.
**              state sequence and probablity of observing a sequence
**              given the model.
**      Organization: University of Maryland
**
**      $Id: nrutil.h,v 1.2 1998/02/19 16:32:42 kanungo Exp kanungo $
**
**	Wan-Ping Lee at the Jackson Lab modified the code on 2017-11-02.
*/

float *vector(int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch);
void free_convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
void nrerror(const char *error_text);

#include <stdio.h>
#include <math.h>

double **init_matrix(int,int);
void free_matrix(double **, int);
double ****init_4d_array(int,int,int,int);
void free_4d_array(double ****,int,int,int);

double **Fae_build(int no, int nv, double **t1, double ****t2, double **Fme, double **fock, double ****ints)
{
  int a,e,m,f,n,i,j,b;
  double **Fae;
  double ****tau_tilde;
  double value;

  tau_tilde = init_4d_array(no,no,nv,nv);
  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++)
          tau_tilde[i][j][a][b] = t2[i][j][a][b] + 0.5*(t1[i][a]*t1[j][b] - t1[i][b]*t1[j][a]);

  Fae = init_matrix(nv,nv);

  for(a=0; a < nv; a++)
    for(e=0; e < nv; e++) {
      value = fock[a+no][e+no];
      for(m=0; m < no; m++) {
        value -= 0.5 * fock[m][e+no] * t1[m][a];
        value -= 0.5 * Fme[m][e] * t1[m][a];
        for(f=0; f < nv; f++) {
          value += t1[m][f] * ints[a+no][m][e+no][f+no];
          for(n=0; n < no; n++)
            value -= 0.5 * tau_tilde[m][n][a][f] * ints[m][n][e+no][f+no];
        }
      }

      Fae[a][e] = value;
    }

  free_4d_array(tau_tilde, no, no, nv);

  return Fae;
}

double **Fmi_build(int no, int nv, double **t1, double ****t2, double **fock, double ****ints)
{
  int m,i,e,n,f;
  double **Fmi;
  double value;

  Fmi = init_matrix(no,no);

  for(m=0; m < no; m++) 
    for(i=0; i < no; i++) {

      value = fock[m][i];

      for(e=0; e < nv; e++)
	value += t1[i][e] * fock[m][e+no];

      for(e=0; e < nv; e++)
	for(n=0; n < no; n++)
	  value += t1[n][e] * ints[m][n][i][e+no];

      for(e=0; e < nv; e++)
	for(n=0; n < no; n++)
	  for(f=0; f < nv; f++)
	    value -= 0.5 * t2[i][n][e][f] * ints[n][m][e+no][f+no];

      for(e=0; e < nv; e++)
	for(n=0; n < no; n++)
	  for(f=0; f < nv; f++)
	    value += t1[n][e]* t1[i][f] * ints[n][m][e+no][f+no];

      Fmi[m][i] = value;
    }

  return Fmi;
}

double **Fme_build(int no, int nv, double **t1, double **fock, double ****ints)
{
  int m,e,n,f;
  double **Fme;
  double value;

  Fme = init_matrix(no,nv);

  for(m=0; m < no; m++)
    for(e=0; e < nv; e++) {

      value = fock[m][e+no];

      for(n=0; n < no; n++)
        for(f=0; f < nv; f++)
          value += t1[n][f] * ints[m][n][e+no][f+no];

      Fme[m][e] = value;
    }

  return Fme;
}

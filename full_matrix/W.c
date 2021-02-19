#include <stdio.h>
#include <math.h>

double ****init_4d_array(int,int,int,int);
void free_4d_array(double ****,int,int,int);

double ****Wmnij_build(int no, int nv, double **t1, double ****t2, double **fock, double ****ints)
{
  int m,n,i,j,e,f,a,b;
  double value;
  double ****Wmnij;
  double ****tau;

  tau = init_4d_array(no,no,nv,nv);
  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++)
          tau[i][j][a][b] = t2[i][j][a][b] + t1[i][a]*t1[j][b] - t1[i][b]*t1[j][a];

  Wmnij = init_4d_array(no,no,no,no);

  for(m=0; m < no; m++)
    for(n=0; n < no; n++)
      for(i=0; i < no; i++)
        for(j=0; j < no; j++) {
          value = ints[m][n][i][j];
          for(e=0; e < nv; e++) {
            value += t1[j][e] * ints[m][n][i][e+no] -
              t1[i][e] * ints[m][n][j][e+no];
            for(f=0; f < nv; f++)
              value += 0.5*tau[i][j][e][f]*ints[m][n][e+no][f+no];
          }
          Wmnij[m][n][i][j] = value;
	}

  free_4d_array(tau, no, no, nv);

  return Wmnij;
}


double ****Wabef_build(int no, int nv, double **t1, double ****t2, double **fock, double ****ints)
{
  int a,b,e,f,m,n,i,j;
  double ****Wabef;
  double ****tau;
  double value;

  tau = init_4d_array(no,no,nv,nv);
  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++)
          tau[i][j][a][b] = t2[i][j][a][b] + t1[i][a]*t1[j][b] - t1[i][b]*t1[j][a];

  Wabef = init_4d_array(nv,nv,nv,nv);

  for(a=0; a < nv; a++)
    for(b=0; b < nv; b++)
      for(e=0; e < nv; e++)
        for(f=0; f < nv; f++) {
          value = ints[a+no][b+no][e+no][f+no];
          for(m=0; m < no; m++) {
            value -= t1[m][b] * ints[a+no][m][e+no][f+no] - t1[m][a] * ints[b+no][m][e+no][f+no];
            for(n=0; n < no; n++)
              value += 0.5 * tau[m][n][a][b] * ints[m][n][e+no][f+no];
          }
          Wabef[a][b][e][f] = value;
        }

  free_4d_array(tau,no,no,nv);

  return Wabef;
}


double ****Wmbej_build(int no, int nv, double **t1, double ****t2, double **fock, double ****ints)
{
  int m,b,e,j,f,n;
  double ****Wmbej;
  double value;

  Wmbej = init_4d_array(no,nv,nv,no);

  for(m=0; m < no; m++)
    for(b=0; b < nv; b++) 
      for(e=0; e < nv; e++)
	for(j=0; j < no; j++) {

	  value = ints[m][b+no][e+no][j];

	  for(n=0; n < no; n++)
	    value -= t1[n][b] * ints[m][n][e+no][j];

	  for(f=0; f < nv; f++)
	    value += t1[j][f] * ints[m][b+no][e+no][f+no];

	  for(n=0; n < no; n++)
	    for(f=0; f < nv; f++)
	      value += t2[j][n][b][f] * ints[m][n][e+no][f+no];

	  for(n=0; n < no; n++)
	    for(f=0; f < nv; f++)
              value -= t1[j][f] * t1[n][b] * ints[m][n][e+no][f+no];

	  Wmbej[m][b][e][j] = value;

	}
  return Wmbej;
}

double ****Wmnie_build(int no, int nv, double **t1, double ****t2, double **fock, double ****ints)
{
  int m,n,i,e,f;
  double ****Wmnie;
  double value;

  Wmnie = init_4d_array(no,no,no,nv);

  for(m=0; m < no; m++)
    for(n=0; n < no; n++)
      for(i=0; i < no; i++)
	for(e=0; e < nv; e++) {

	  value = ints[m][n][i][e+no];

	  for(f=0; f < nv; f++)
	    value -= t1[i][f] * ints[m][n][e+no][f+no];

	  Wmnie[m][n][i][e] = value;

	}

  return Wmnie;
}

double ****Wamef_build(int no, int nv, double **t1, double ****t2, double **fock, double ****ints)
{
  int a,m,e,f,n;
  double ****Wamef;
  double value;

  Wamef = init_4d_array(nv,no,nv,nv);

  for(a=0; a < nv; a++)
    for(m=0; m < no; m++)
      for(e=0; e < nv; e++)
	for(f=0; f < nv; f++) {

	  value = ints[a+no][m][e+no][f+no];

	  for(n=0; n < no; n++)
	    value += t1[n][a] * ints[m][n][e+no][f+no];

	  Wamef[a][m][e][f] = value;

	}

  return Wamef;
}

double ****Wmbij_build(int no, int nv, double **t1, double ****t2, double **Fme, double ****Wmnij, double ****ints)
{
  int m,b,i,j,a,e,n,f;
  double ****Wmbij;
  double ****tau;
  double value;

  tau = init_4d_array(no,no,nv,nv);
  for(i=0; i < no; i++)
    for(j=0; j < no; j++) 
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++)
          tau[i][j][a][b] = t2[i][j][a][b] + t1[i][a]*t1[j][b] -  t1[i][b]*t1[j][a];

  Wmbij = init_4d_array(no,nv,no,no);

  for(m=0; m < no; m++)
    for(b=0; b < nv; b++)
      for(i=0; i < no; i++)
        for(j=0; j < no; j++) {
          value = ints[m][b+no][i][j];
          for(e=0; e < nv; e++)
            value -= Fme[m][e] * t2[i][j][b][e];
          for(n=0; n < no; n++)
            value -= t1[n][b] * Wmnij[m][n][i][j];
          for(e=0; e < nv; e++)
            for(f=0; f < nv; f++) 
              value += 0.5 * ints[m][b+no][e+no][f+no] * tau[i][j][e][f];
          for(n=0; n < no; n++)
            for(e=0; e < nv; e++)
              value += ints[m][n][i][e+no] * t2[j][n][b][e];
          for(n=0; n < no; n++)
            for(e=0; e < nv; e++)
              value -= ints[m][n][j][e+no] * t2[i][n][b][e];
          for(e=0; e < nv; e++)
            value += t1[i][e] * ints[m][b+no][e+no][j] -
              t1[j][e] * ints[m][b+no][e+no][i];
          for(e=0; e < nv; e++)
            for(n=0; n < no; n++)
              for(f=0; f < nv; f++)
                value -= t1[i][e] * t2[n][j][b][f] *
                  ints[m][n][e+no][f+no] -
                  t1[j][e] * t2[n][i][b][f] *
                  ints[m][n][e+no][f+no];

          Wmbij[m][b][i][j] = value;

        }

  free_4d_array(tau, no, no, nv);

  return Wmbij;
}

double ****Wabei_build(int no, int nv, double **t1, double ****t2, double **Fme, double ****Wmnie, double ****ints)
{
  int a,b,e,i,j,m,f,n;
  double ****Wabei;
  double ****tau, ****tau2;
  double value;

  tau = init_4d_array(no,no,nv,nv);
  tau2 = init_4d_array(no,no,nv,nv);
  for(i=0; i < no; i++)
    for(j=0; j < no; j++) 
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++) {
          tau[i][j][a][b] = t2[i][j][a][b] + t1[i][a]*t1[j][b] -  t1[i][b]*t1[j][a];
          tau2[i][j][a][b] = t2[i][j][a][b] + t1[i][a]*t1[j][b];
        }

  Wabei = init_4d_array(nv,nv,nv,no);

  for(a=0; a < nv; a++)
    for(b=0; b < nv; b++)
      for(e=0; e < nv; e++)
        for(i=0; i < no; i++) {

          /** Term I **/
          value = ints[a+no][b+no][e+no][i];

          /** Term II **/
          for(m=0; m < no; m++)
            value -= Fme[m][e] * t2[m][i][a][b];

          /** Term IIIa **/
          for(f=0; f < nv; f++)
            value += t1[i][f] * ints[a+no][b+no][e+no][f+no];

          /** Terms IIIc + IIId + IV **/
          for(m=0; m < no; m++)
            for(n=0; n < no; n++)
              value -= 0.5 * Wmnie[m][n][i][e] * tau[m][n][a][b];

          /** Terms IIIb + V **/          
          for(m=0; m < no; m++)
            for(f=0; f < nv; f++)
              value -= tau2[m][i][b][f] * ints[a+no][m][e+no][f+no] - tau2[m][i][a][f] * ints[b+no][m][e+no][f+no];

          /** Term VI **/
          for(m=0; m < no; m++)
            value -= t1[m][a] * ints[m][b+no][e+no][i] - t1[m][b] * ints[m][a+no][e+no][i];

          /** Term VII **/
          for(m=0; m < no; m++)
            for(n=0; n < no; n++)
              for(f=0; f < nv; f++)
                value += ( t1[m][a] * t2[n][i][b][f] - t1[m][b] * t2[n][i][a][f] ) * ints[m][n][e+no][f+no];

          Wabei[a][b][e][i] = value;
        }

   free_4d_array(tau, no, no, nv);
   free_4d_array(tau2, no, no, nv);

  return Wabei;
}

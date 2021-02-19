#include <stdio.h>
#include <math.h>

double ******init_6d_array(int,int,int,int,int,int);
void free_6d_array(double ******,int,int,int,int,int);

double ******Wmabeij_build(int no, int nv, double ****t2, double ****Wamef, double ****Wmnie)
{
  int m,a,b,e,i,j,f,n;
  double value;
  double ******Wmabeij;

  Wmabeij = init_6d_array(no,nv,nv,nv,no,no);

  for(m=0; m < no; m++)
    for(a=0; a < nv; a++)
      for(b=0; b < nv; b++)
        for(e=0; e < nv; e++)
          for(i=0; i < no; i++)
            for(j=0; j < no; j++) {
              value=0.0;
              for(f=0; f < nv; f++) {
                value -= Wamef[a][m][e][f] * t2[i][j][f][b];
                value += Wamef[b][m][e][f] * t2[i][j][f][a];
              }
              for(n=0; n < no; n++) {
                value += Wmnie[m][n][i][e] * t2[n][j][a][b];
                value -= Wmnie[m][n][j][e] * t2[n][i][a][b];
              }
              Wmabeij[m][a][b][e][i][j] = value;
            }

  return Wmabeij;
}

double ******Wmaneij_build(int no, int nv, double ****t2, double ****ints)
{
  int m,a,n,e,i,j,f;
  double value;
  double ******Wmaneij;

  Wmaneij = init_6d_array(no,nv,no,nv,no,no);

  for(m=0; m < no; m++)
    for(a=0; a < nv; a++)
      for(n=0; n < no; n++)
        for(e=0; e < nv; e++)
          for(i=0; i < no; i++)
            for(j=0; j < no; j++) {
              value = 0.0;
              for(f=0; f < nv; f++)
                value += ints[m][n][e+no][f+no] * t2[i][j][a][f];

              Wmaneij[m][a][n][e][i][j] = value;
            }

  return Wmaneij;
}

double ******Wmabeif_build(int no, int nv, double ****t2, double ****ints)
{
  int m,a,b,e,i,f,n;
  double value;
  double ******Wmabeif;

  Wmabeif = init_6d_array(no,nv,nv,nv,no,nv);

  for(m=0; m < no; m++)
    for(a=0; a < nv; a++)
      for(b=0; b < nv; b++)
        for(e=0; e < nv; e++)
          for(i=0; i < no; i++)
            for(f=0; f < nv; f++) {
              value = 0.0;
              for(n=0; n < no; n++)
                value -= ints[m][n][e+no][f+no] * t2[i][n][a][b];

              Wmabeif[m][a][b][e][i][f] = value;
            }

  return Wmabeif;
}


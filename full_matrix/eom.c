/* EOM: Computes EOM-CCSD excitation energy via explicit diagonalization of
** HBAR.
**
** TDC, 3/3/09
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *init_array(int);
double **init_matrix(int,int);
void free_matrix(double **, int);
void print_mat(double **a, int m, int n,FILE *out);
void zero_matrix(double **, int, int);
void zero_array(double *, int);
double ****init_4d_array(int,int,int,int);
void free_4d_array(double ****,int,int,int);
void diag(int nm, int n, double **array, double *e_vals, int matz,
        double **e_vecs, double toler);

double **block_matrix(int, int);
void free_block(double **);
void sort_vector(double *, int);

extern int DGEEV_(char *, char *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);

double **Fae_build(int,int,double **,double ****,double **,double **,double ****);
double **Fmi_build(int,int,double **,double ****,double **,double ****);
double **Fme_build(int,int,double **,double **,double ****);
double ****Wmbej_build(int,int,double **,double ****,double **,double ****);
double ****Wabef_build(int,int,double **,double ****,double **,double ****);
double ****Wmnij_build(int,int,double **,double ****,double **,double ****);
double ****Wmnie_build(int,int,double **,double ****,double **,double ****);
double ****Wamef_build(int,int,double **,double ****,double **,double ****);
double ****Wmbij_build(int,int,double **,double ****,double **,double ****,double ****);
double ****Wabei_build(int,int,double **,double ****,double **,double ****,double ****);
double ******Wmabeij(int no, int nv, double ****t2, double ****Wamef, double ****Wmnie);
double ******Wmaneij(int no, int nv, double ****t2, double ****ints);
double ******Wmabeif(int no, int nv, double ****t2, double ****ints);

int main(int argc, char *argv[])
{
  FILE *input;
  int no, nv, nmo;
  double **fock, ****ints;
  double **t1, ****t2;
  double Escf, Ecc, Emp2;
  int p, q, r, s;
  int i, j, k, l, a, b, c, d;
  double **Fae, **Fmi, **Fme; 
  double ****Wmbej, ****Wabef, ****Wmnij; 
  double ****Wamef, ****Wmnie, ****Wmbij, ****Wabei;
  double ******Wmabeij;
  double ******Wmaneij;
  double ******Wmabeif;
  int dim, M, N, lwork, info;
  double **HBAR, *wr, *wi, **vl, **vr, *work;
  char junk;

  /* Grab the CC Data */
  input = fopen("CC_data.txt", "r");
  fscanf(input, "%d %d", &no, &nv);
  fscanf(input, "%lf", &Escf);
  fscanf(input, "%lf", &Ecc);

  nmo = no + nv;

  fock = init_matrix(nmo,nmo);
  ints = init_4d_array(nmo,nmo,nmo,nmo);
  t1 = init_matrix(no,nv);
  t2 = init_4d_array(no,no,nv,nv);

  for(p=0; p < nmo; p++)
    for(q=0; q < nmo; q++)
      fscanf(input, "%lf", &fock[p][q]);

  for(p=0; p < nmo; p++)
    for(q=0; q < nmo; q++)
      for(r=0; r < nmo; r++)
        for(s=0; s < nmo; s++)
          fscanf(input, "%lf", &ints[p][q][r][s]);

  for(i=0; i < no; i++)
    for(a=0; a < nv; a++)
      fscanf(input, "%lf", &t1[i][a]);

  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++)
          fscanf(input, "%lf", &t2[i][j][a][b]);

  fclose(input);

  /* MP2 Energy Test */
  Emp2=0.0;
  for(i=0; i < no; i++)
    for(j=0; j < i; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < a; b++)
          Emp2 += ints[i][j][a+no][b+no] * ints[i][j][a+no][b+no]/(fock[i][i] + fock[j][j] - fock[a+no][a+no] - fock[b+no][b+no]);

  printf("E(MP2) = %20.12f\n", Emp2);

  /* HBAR Components */
  Fme = Fme_build(no, nv, t1, fock, ints);
  Fae = Fae_build(no, nv, t1, t2, Fme, fock, ints);
  Fmi = Fmi_build(no, nv, t1, t2, fock, ints);
  Wmbej = Wmbej_build(no, nv, t1, t2, fock, ints);
  Wabef = Wabef_build(no, nv, t1, t2, fock, ints);
  Wmnij = Wmnij_build(no, nv, t1, t2, fock, ints);
  Wmnie = Wmnie_build(no, nv, t1, t2, fock, ints);
  Wamef = Wamef_build(no, nv, t1, t2, fock, ints);
  Wmbij = Wmbij_build(no, nv, t1, t2, Fme, Wmnij, ints);
  Wabei = Wabei_build(no, nv, t1, t2, Fme, Wmnie, ints);
  Wmabeij = Wmabeij_build(no, nv, t2, Wamef, Wmnie);
  Wmaneij = Wmaneij_build(no, nv, t2, ints);
  Wmabeif = Wmabeif_build(no, nv, t2, ints);

  /* Build HBAR as one big matrix */
  dim = no*nv + (no*(no-1)/2)*(nv*(nv-1)/2);
  HBAR = block_matrix(dim,dim);

  /* HBAR SS */
  for(i=0,M=0; i < no; i++)
    for(a=0; a < nv; a++,M++) {

      for(j=0,N=0; j < no; j++)
        for(b=0; b < nv; b++,N++) {
          HBAR[M][N] += (i==j) * Fae[a][b];
          HBAR[M][N] -= (a==b) * Fmi[j][i];
          HBAR[M][N] += Wmbej[j][a][b][i];
        }
    }

  /* HBAR SD */
  for(i=0,M=0; i < no; i++)
    for(a=0; a < nv; a++,M++) {

      for(j=0,N=no*nv; j < no; j++)
        for(k=0; k < j; k++)
          for(b=0; b < nv; b++)
            for(c=0; c < b; c++,N++) {
              HBAR[M][N] -= (i==j) * (a==c) * Fme[k][b];
              HBAR[M][N] += (i==j) * (a==b) * Fme[k][c];
              HBAR[M][N] += (i==k) * (a==c) * Fme[j][b];
              HBAR[M][N] -= (i==k) * (a==b) * Fme[j][c];

              HBAR[M][N] += (a==c) * Wmnie[j][k][i][b];
              HBAR[M][N] -= (a==b) * Wmnie[j][k][i][c];

              HBAR[M][N] += (i==j) * Wamef[a][k][b][c];
              HBAR[M][N] -= (i==k) * Wamef[a][j][b][c];
            }
    }

  /* HBAR DS */
  for(i=0,M=no*nv; i < no; i++)
    for(j=0; j < i; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < a; b++,M++) {

          for(k=0,N=0; k < no; k++)
            for(c=0; c < nv; c++,N++) {
             HBAR[M][N] += (
              - (j == k) * Wabei[a][b][c][i]
              + (i == k) * Wabei[a][b][c][j]
              + (b == c) * Wmbij[k][a][i][j] 
              - (a == c) * Wmbij[k][b][i][j]
             );

             HBAR[M][N] += Wmabeij[k][a][b][c][i][j];
            }
        }

  /* HBAR DD */
  for(i=0,M=no*nv; i < no; i++)
    for(j=0; j < i; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < a; b++,M++) {

          for(k=0,N=no*nv; k < no; k++)
            for(l=0; l < k; l++)
            for(c=0; c < nv; c++)
                for(d=0; d < c; d++,N++) {
                  HBAR[M][N] += (
                   - (b == d) * (i == l) * (j == k) * Fae[a][c]
                   + (b == d) * (i == k) * (j == l) * Fae[a][c]
                   + (b == c) * (i == l) * (j == k) * Fae[a][d]
                   - (b == c) * (i == k) * (j == l) * Fae[a][d]
                   + (a == d) * (i == l) * (j == k) * Fae[b][c]
                   - (a == d) * (i == k) * (j == l) * Fae[b][c]
                   - (a == c) * (i == l) * (j == k) * Fae[b][d]
                   + (a == c) * (i == k) * (j == l) * Fae[b][d]
                  );

                  HBAR[M][N] += (
                   + (a == d) * (b == c) * (j == l) * Fmi[k][i]
                   - (a == c) * (b == d) * (j == l) * Fmi[k][i]
                   - (a == d) * (b == c) * (i == l) * Fmi[k][j]
                   + (a == c) * (b == d) * (i == l) * Fmi[k][j]
                   - (a == d) * (b == c) * (j == k) * Fmi[l][i]
                   + (a == c) * (b == d) * (j == k) * Fmi[l][i]
                   + (a == d) * (b == c) * (i == k) * Fmi[l][j]
                   - (a == c) * (b == d) * (i == k) * Fmi[l][j]
                  );

                  HBAR[M][N] += (
                   - (i == l) * (j == k) * Wabef[a][b][c][d]
                   + (i == k) * (j == l) * Wabef[a][b][c][d]

                   - (a == d) * (b == c) * Wmnij[k][l][i][j]
                   + (a == c) * (b == d) * Wmnij[k][l][i][j]

                   + (b == d) * (j == l) * Wmbej[k][a][c][i]
                   - (b == c) * (j == l) * Wmbej[k][a][d][i]
                   - (b == d) * (i == l) * Wmbej[k][a][c][j]
                   + (b == c) * (i == l) * Wmbej[k][a][d][j]
                   - (a == d) * (j == l) * Wmbej[k][b][c][i]
                   + (a == c) * (j == l) * Wmbej[k][b][d][i]
                   + (a == d) * (i == l) * Wmbej[k][b][c][j]
                   - (a == c) * (i == l) * Wmbej[k][b][d][j]
                   - (b == d) * (j == k) * Wmbej[l][a][c][i]
                   + (b == c) * (j == k) * Wmbej[l][a][d][i]
                   + (b == d) * (i == k) * Wmbej[l][a][c][j]
                   - (b == c) * (i == k) * Wmbej[l][a][d][j]
                   + (a == d) * (j == k) * Wmbej[l][b][c][i]
                   - (a == c) * (j == k) * Wmbej[l][b][d][i]
                   - (a == d) * (i == k) * Wmbej[l][b][c][j]
                   + (a == c) * (i == k) * Wmbej[l][b][d][j]
                  );

                  HBAR[M][N] += (
                   + (j == l) * Wmabeif[k][a][b][c][i][d] 
                   - (i == l) * Wmabeif[k][a][b][c][j][d] 
                   - (j == k) * Wmabeif[l][a][b][c][i][d] 
                   + (i == k) * Wmabeif[l][a][b][c][j][d]

                   - (b == d) * Wmaneij[k][a][l][c][i][j] 
                   + (b == c) * Wmaneij[k][a][l][d][i][j] 
                   + (a == d) * Wmaneij[k][b][l][c][i][j] 
                   - (a == c) * Wmaneij[k][b][l][d][i][j] 
                  );
                }
        }

  /* Diagonalize HBAR */
  wr = init_array(dim); wi = init_array(dim); vl = block_matrix(dim,dim);
  vr = block_matrix(dim,dim); work = init_array(20*dim);
  lwork = 20*dim; info = 0; junk = 'V';
  DGEEV_(&junk, &junk, &dim, HBAR[0], &dim, wr, wi, vl[0], &dim, vr[0], &dim, work, &lwork, &info);
  sort_vector(wr, dim);

  printf("EOM-CCSD Excitation Energies:\n");
  printf(" #        Hartree                  eV         \n");
  printf("--  --------------------  --------------------\n");
  for(M=0; M < dim; M++) printf("%2d  %20.10f  %20.10f\n", M, wr[M], wr[M]*27.211);
  free(wr); free(wi); free_block(vl); free_block(vr); free(work);

  free_block(HBAR);

  free_6d_array(Wmaneij, no, nv, no, nv, no);
  free_6d_array(Wmabeif, no, nv, nv, nv, no);
  free_6d_array(Wmabeij, no, nv, nv, nv, no);
  free_4d_array(Wabei, nv, nv, nv);
  free_4d_array(Wmbij, no, nv, no);
  free_4d_array(Wamef, nv, no, nv);
  free_4d_array(Wmnie, no, no, no);
  free_4d_array(Wmnij, no, no, no);
  free_4d_array(Wabef, nv, nv, nv);
  free_4d_array(Wmbej, no, nv, nv);
  free_matrix(Fme, no);
  free_matrix(Fmi, no);
  free_matrix(Fae, nv);

  free_matrix(fock, nmo);
  free_4d_array(ints, nmo, nmo, nmo);
  free_matrix(t1,no);
  free_4d_array(t2,no,no,nv);

  return 0;
}

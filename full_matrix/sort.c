/*!
** sort_vector(): Sort the elements of a vector into increasing order
**
** \param A = array to sort
** \param n = length of array
**
** Returns: none
** \ingroup QT
*/
void sort_vector(double *A, int n)
{
  int i, j, k;
  double val;

  for(i=0; i < n-1; i++) {
    val = A[k=i];

    for(j=i+1; j < n; j++)
      if(A[j] <= val) val = A[k=j];

    if(k != i) {
      A[k] = A[i];
      A[i] = val;
    }
  }
}


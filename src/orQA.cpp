#include <Rcpp.h>

#include <functional>

extern "C" {
void C_isomean(double* y, double* w, int* n, double* ghat)
{
  int c, j, nn; 
  double neu;
  double* gew;
  int* k;
  
  nn = *n; /* nn > 1 */
 
  /* allocate vector - error handling is done by R */
  k = (int *) Calloc((size_t) nn, int);
  gew = (double *) Calloc((size_t) nn, double);

  c = 0;
  k[c] = 0;
  gew[c] = w[0];
  ghat[c] = y[0];

  for (j=1; j < nn; j++)
  {  
    c = c+1;
    k[c] = j;
    gew[c] = w[j];
    ghat[c] = y[j];

    /* c is at least 1 as nn is > 1 */
    while (ghat[c-1] >= ghat[c])
    {
      neu = gew[c]+gew[c-1];
      ghat[c-1] = ghat[c-1]+(gew[c]/neu)*(ghat[c]-ghat[c-1]);
      gew[c-1] = neu;
      c = c-1;

      if (c==0) break;
    }
  }

  while (nn >= 1)
  {
    for (j=k[c]; j < nn; j++)
    {
      ghat[j] = ghat[c];
    }
    nn = k[c];
    c = c-1;
  }

  /* free vector */
  Free(k); 
  Free(gew); 
}
}

extern "C" {
  SEXP misoreg ( SEXP data, SEXP weights );
}

SEXP misoreg ( SEXP data, SEXP weights ) {
SEXP rl=R_NilValue;
RcppMatrix<double> input(data);
RcppVector<double> w(weights);
RcppMatrix<double> output(data);
int n,m;
m= input.getDim2();
n= input.getDim1();
std::vector<double> stlwe(w.stlVector());
std::vector< std::vector<double> > stlin(input.stlMatrix());
std::vector< std::vector<double> > stlo(output.stlMatrix());
for(int i=0; i<n; i++){
C_isomean(&(stlin[i][0]),&stlwe[0],&m,&(stlo[i][0]));
}
RcppResultSet rs;
rs.add("result",stlo);
rs.add("weights",stlwe);
rl = rs.getReturnList();
return(rl);
  Rf_warning("your C program does not return anything!");
  return R_NilValue;
}


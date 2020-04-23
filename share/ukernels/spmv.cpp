



void matvec(Matrix *A, double *x, double *y)
{//
#pragma omp taskloopprivate(j,is,ie,j0,y0) grain_size(500)
  for (i = 0; i < A->n; i++)
    {
      y0 = 0;is = A->ptr[i];
      ie = A->ptr[i + 1];
      for (j = is; j < ie; j++)
	{
	  j0 = index[j];
	  y0 += value[j] * x[j0];
	}

      y[i] = y0;
    }
}


int main(int argc, char** argv)
{
  
#pragma omp parallel
#pragma omp single for
  (iter= 0; iter< sc-  >maxIter; iter++) {precon(A, r, z);vectorDot(r, z, n, &rho);beta = rho / rho_old;xpay(z, beta, n, p);matvec(A, p, q);vectorDot(p, q, n, &dot_pq);alpha = rho / dot_pq;axpy(alpha, p, n, x);axpy(-alpha, q, n, r);sc-  >residual = sqrt(rho) * bnrm2;if (sc-  >residual <= sc->tolerance) break;rho_old= rho;}void matvec(Matrix *A, double *x, double *y) {// ...#pragma omp taskloopprivate(j,is,ie,j0,y0) \grain_size(500)for (i = 0; i < A->n; i++) {y0 = 0;is = A->ptr[i];ie = A->ptr[i + 1];for (j = is; j < ie; j++) {j0 = index[j];y0 += value[j] * x[j0];}y[i] = y0;}// ...}


  #pragma omp parallel forfor (int i=0; i<n; ++i) {y[i] = 0.0;for (int j=0; j<n_diag; ++j) {int start, v_offset, end = n_row();if (diag[j] < 0) {start -= diag[j];v_offset = -start;} else {end -= diag[j];v_offset = diag[j];}ind = j*n_row + i;if ( (i >= start) && (i < end)) {y[i] += val[ind] * x[i+v_offset];}}

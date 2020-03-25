//
// Created by kazem on 12/4/18.
//

#ifndef PROJECT_GMRES_LDLT_H
#define PROJECT_GMRES_LDLT_H

#include <stdio.h>
#include <string.h>
#include <complex.h>

#define conj(A)           (A)

int gmres_ldlt(
  int n, size_t *Ap, int *Ai, double *Ax,
   size_t *Lp, int *Li, double *Lx, size_t NNZ,
  size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
  double *d_val,
  double *x, double *b,
  bool precond, double eps) {
 double *gmHi,*gmH;
 double *gmVi,*gmV;
 double *gmWi,*gmW;
 double *gmsin, *gmcos;
 double *gmG;
 double  tmp;
 //pastix_fixdbl_t     t0, t3;
 double              resid, resid_b;
 double              norm, normb, normx;
 int                 im, im1, itermax;
 int                 i, j,  ldw, iters;
 int                 outflag, inflag;
 int                 savemem = 0;
 double alp[2] = {-1.0, 0};
 double bet[2] = {1.0, 0};
 const int iun = 1;

 /* Get the parameters */
 //GMRES restart parameter
 im      =  25; //pastix_data->iparm[IPARM_GMRES_IM];
 im1     = im + 1;
 itermax = 250; //pastix_data->iparm[IPARM_ITERMAX];
 ldw     = n;

 /*if ( !(pastix_data->steps & STEP_NUMFACT) ) {
  precond = 0;
 }*/

 if ((!precond) || savemem ) {
  ldw = 0;
 }

 gmcos = (double *)malloc(im  * sizeof(double));
 gmsin = (double *)malloc(im  * sizeof(double));
 gmG   = (double *)malloc(im1 * sizeof(double));

 /**
  * H stores the h_{i,j} elements ot the upper hessenberg matrix H (See Alg. 9.5 p 270)
  * V stores the v_{i} vectors
  * W stores the M^{-1} v_{i} vectors to avoid the application of the
  *          preconditioner on the output result (See line 11 of Alg 9.5)
  *
  * If no preconditioner is applied, or the user wants to save memory, W
  * stores only temporarily one vector for the Ax product (ldw is set to 0 to
  * reuse the same vector at each iteration)
  */
 gmH = (double *)malloc(im * im1 * sizeof(double));
 gmV = (double *)malloc(n  * im1 * sizeof(double));
 if (precond && (!savemem) ) {
  gmW = (double *)malloc(n * im  * sizeof(double));
 }
 else {
  gmW = (double *)malloc(n       * sizeof(double));
 }
 memset( gmH, 0, im * im1 * sizeof(double) );

#if defined(PASTIX_DEBUG_GMRES)
 dbg_x = (double *)malloc(n   * sizeof(double));
    dbg_r = (double *)malloc(n   * sizeof(double));
    dbg_G = (double *)malloc(im1 * sizeof(double));
    copy( n, x, dbg_x );
#endif

 normb = norm_dense(1,n,b,2);
 normx = norm_dense(1,n,x,2);
 //normb = solver.norm( pastix_data, n, b );
 //normx = solver.norm( pastix_data, n, x );


 /**
  * Algorithm from Iterative Methods for Sparse Linear systems, Y. Saad, Second Ed. p267-273
  *
  * The version implemented is the Right preconditioned algorithm.
  */
 outflag = 1;
 iters = 0;
 while (outflag) {
  /* Initialize v_{0} and w_{0} */
  gmVi = gmV;
  /* Compute r0 = b - A * x */
  for (int i = 0; i < n; ++i) {
   gmVi[i] = b[i];
  }
  if ( normx > 0. ) {
   spmv_csc_sym(n, Ap, Ai, Ax, -1, alp, bet, 1, x, gmVi);
  }

  //solver.copy( pastix_data, n, b, gmVi );
  //if ( normx > 0. ) {
  // solver.spmv( pastix_data, PastixNoTrans, -1., x, 1., gmVi );
  //}

  /* Compute resid = ||r0||_f */
  resid = norm_dense(1,n,gmVi,2);
  resid_b = resid / normb;
  //resid = solver.norm( pastix_data, n, gmVi );

  /* If residual is small enough, exit */
  if ( resid_b <= eps ) {
   outflag = 0;
   break;
  }

  /* Compute v0 = r0 / resid */
  tmp = (double)( 1.0 / resid );
  VEC_SCAL( n, tmp, gmVi, iun );

  gmG[0] = (double)resid;
  inflag = 1;
  i = -1;
  gmHi = gmH - im1;
  gmWi = gmW - ldw;
/* Starts from here */
  while( inflag ){
   i++;
   /* Set H and W pointers to the beginning of columns i */
   gmHi = gmHi + im1;
   gmWi = gmWi + ldw;

   /* Backup v_{i} into w_{i} for the end */
   for (int k = 0; k < n; ++k) {
    gmWi[k] = gmVi[k];
   }
   //solver.copy( pastix_data, n, gmVi, gmWi );

   /* Compute w_{i} = M^{-1} v_{i} */
   if ( precond ) {
    //solver.spsv( pastix_data, gmWi );
    solve_phase_ldl(n, d_val, gmWi, col2sup, sup2col,
                    Lp, Li, Lx, Li_ptr, supNo, NNZ);
   }
   /* v_{i+1} = A (M^{-1} v_{i}) = A w_{i} */
   gmVi += n;
   //solver.spmv( pastix_data, PastixNoTrans, 1.0, gmWi, 0., gmVi );

   /* Classical Gram-Schmidt */
   for (j=0; j<=i; j++) {
    /* Compute h_{j,i} = < v_{i+1}, v_{j} > */
    //gmHi[j] = solver.dot( pastix_data, n, gmVi, gmV + j * n );
    gmHi[j] = dot(n, gmVi, gmV+(j*n) );
    /* Compute v_{i+1} = v_{i+1} - h_{j,i} v_{j} */
    //solver.axpy( pastix_data, n, -1. * gmHi[j],  gmV + j * n, gmVi );
    double coef = -1. * gmHi[j];
    daxpy(&n,&coef, gmV + j * n, &iun, gmV, &iun);
   }

   /* Compute || v_{i+1} ||_f */
   //norm = solver.norm( pastix_data, n, gmVi );
   norm = norm_dense(1,n,gmVi,2);
   gmHi[i+1] = norm;

   /* Compute v_{i+1} = v_{i+1} / h_{i+1,i} iff h_{i+1,i} is not too small */
   if ( norm > 1e-50 ){
    tmp = (double)(1.0 / norm);
    //solver.scal( pastix_data, n, tmp, gmVi );
    VEC_SCAL( n, tmp, gmVi, iun );
   }

   /* Apply the previous Givens rotation to the new column (should call LAPACKE_zrot_work())*/
   for (j=0; j<i;j++){
    /*
     * h_{j,  i} = cos_j * h_{j,  i} +      sin_{j}  * h_{j+1, i}
     * h_{j+1,i} = cos_j * h_{j+1,i} - conj(sin_{j}) * h_{j,   i}
     */
    tmp = gmHi[j];
    gmHi[j]   = gmcos[j] * tmp       +      gmsin[j]  * gmHi[j+1];
    gmHi[j+1] = gmcos[j] * gmHi[j+1] - conj(gmsin[j]) * tmp;
   }

   /*
    * Compute the new Givens rotation (zrotg)
    *
    * t = sqrt( h_{i,i}^2 + h_{i+1,i}^2 )
    * cos = h_{i,i}   / t
    * sin = h_{i+1,i} / t
    */
   {
    tmp = csqrt( gmHi[i]   * gmHi[i] +
                 gmHi[i+1] * gmHi[i+1] );

    if ( cabs(tmp) <= eps ) {
     tmp = (double)eps;
    }
    gmcos[i] = gmHi[i]   / tmp;
    gmsin[i] = gmHi[i+1] / tmp;
   }
   /* Update the residuals (See p. 168, eq 6.35) */
   gmG[i+1] = -gmsin[i] * gmG[i];
   gmG[i]   =  gmcos[i] * gmG[i];
   /* Apply the last Givens rotation */
   gmHi[i] = gmcos[i] * gmHi[i] + gmsin[i] * gmHi[i+1];
   /* (See p. 169, eq 6.42) */
   resid = cabs( gmG[i+1] );
   resid_b = resid / normb;
   iters++;
   if ( (i+1 >= im) ||
        (resid_b <= eps) ||
        (iters >= itermax) ) {
    inflag = 0;
   }
  /* if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
    solver.output_oneiter( t0, t3, resid_b, iters );*/

#ifdef DEBUG_GMRES
    {
                    double normr2;
                    /* Compute y_m = H_m^{-1} g_m (See p. 169) */
                    memcpy( dbg_G, gmG, im1 * sizeof(pastix_complex64_t) );
                    cblas_ztrsv( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                 i+1, gmH, im1, dbg_G, 1 );

                    solver.copy( pastix_data, n, b, dbg_r );
                    solver.copy( pastix_data, n, x, dbg_x );

                    /* Accumulate the current v_m */
                    solver.gemv( pastix_data, n, i+1, 1.0, (precond ? gmW : gmV), n, dbg_G, 1.0, dbg_x );

                    /* Compute b - Ax */
                    solver.spmv( pastix_data, PastixNoTrans, -1., dbg_x, 1., dbg_r );

                    normr2 = solver.norm( pastix_data, n, dbg_r );
                    fprintf(stdout, OUT_ITERREFINE_ERR, normr2 / normb );
                }
#endif
//   }
  }

  /* Compute y_m = H_m^{-1} g_m (See p. 169) */
  /*cblas_ztrsv( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
               i+1, gmH, im1, gmG, 1 );*/
  dtrsv("U", "N", "N", i+1, gmH,
        im1,gmG,&iun);


  /**
   * Compute x_m = x_0 + M^{-1} V_m y_m
   *             = x_0 +        W_m y_m
   */
  if (precond && savemem) {
   /**
    * Since we saved memory, we do not have (M^{-1} V_m) stored,
    * thus we compute:
    *     w = V_m y_m
    *     w = M^{-1} (V_m y_m)
    *     x = x0 + (M^{-1} (V_m y_m))
    */
   solver.gemv( pastix_data, n, i+1, 1.0, gmV, n, gmG, 0., gmW );
   solver.spsv( pastix_data, gmW );
   solver.axpy( pastix_data, n, 1.,  gmW, x );
  }
  else {
   /**
    * Since we did not saved memory, we do have (M^{-1} V_m) stored in
    * W_m if precond is true, thus we compute:
    *     x = x0 + W_m y_m, if precond
    *     x = x0 + V_m y_m, if not precond
    */
   gmWi = precond ? gmW : gmV;
   //solver.gemv( pastix_data, n, i+1, 1.0, gmWi, n, gmG, 1.0, x );
   dgemv("N", n, i+1, bet, gmWi, n, gmG, iun, bet, y, iun);
  }

  if ((resid_b <= eps) || (iters >= itermax)){
   outflag = 0;
  }
 }

#if defined(PASTIX_DEBUG_GMRES)
 solver.free(dbg_x);
    solver.free(dbg_r);
    solver.free(dbg_G);
#endif

 return iters;
}

#endif //PROJECT_GMRES_LDLT_H

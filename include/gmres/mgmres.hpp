//  mgmres.hpp
//
//  Thanks to Philip Sakievich for correcting mismatches between this HPP
//  file and the corresponding CPP file, 23 March 2016!
//
# include <cstddef>

namespace nasoq {
//  Thanks to Philip Sakievich for correcting mismatches between this HPP
//  file and the corresponding CPP file, 23 March 2016!
//

 void ax_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[],
            double w[]);

 void diagonal_pointer_cr(int n, int nz_num, int ia[], int ja[], int ua[]);

 void mult_givens(double c, double s, int k, double g[]);

 double r8vec_dot(int n, double a1[], double a2[]);

 void rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[]);

 void timestamp();

/******************************************************************************/

 double **dmatrix(int nrl, int nrh, int ncl, int nch);

 double **dmatrix_prealloc(int nrl, int nrh, int ncl, int nch, double *ws);

 void free_dmatrix_prealloc(double **m, int nrl, int nrh, int ncl, int nch);

/******************************************************************************/

 void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);


 void lus_cr(int n, int nz_num, int ia[], int ja[], double l[], int ua[],
             double r[], double z[]);




/******************************************************************************/

 int pmgmres_ldlt_cr(int n, int nz_num, int ia[], int ja[], double a[],
                     size_t *Lp, int *Li, double *Lx, size_t NNZ,
                     size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
                     double *d_val,
                     double x[], double rhs[], int itr_max, int mr,
                     double tol_abs, double tol_rel, int sorted = 1, int blocked = 0,
                     int levels = 0, int *levelPtr = NULL, int *levelSet = NULL,
                     int parts = 0, int *parPtr = NULL, int *partition = NULL,
                     int chunk = 0);


 void free_ws_ir(double *ws, double *c, double *g, double **h, int mr,
                 double *r, double *s, double **v, int n, double *y);

 /******************************************************************************/

 int pmgmres_ldlt_auto(int n, int nz_num, int ia[], int ja[], double a[],
                       size_t *Lp, int *Li, double *Lx, size_t NNZ,
                       size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
                       double *d_val,
                       double x[], double rhs[], int itr_max, int mr,
                       double tol_abs, double tol_rel, int sorted = 1, int blocked = 0,
                       int levels = 0, int *levelPtr = NULL, int *levelSet = NULL,
                       int parts = 0, int *parPtr = NULL, int *partition = NULL,
                       int chunk = 0, double *ws = NULL);


/******************************************************************************/

 int pmgmres_ldlt_auto_update(int n, int nz_num, int ia[], int ja[], double a[],
                              size_t *Lp, int *Li, double *Lx, size_t NNZ,
                              size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
                              double *d_val,
                              double x[], double rhs[], int itr_max, int mr,
                              double tol_abs, double tol_rel, bool *mask, int *mask_col,
                              int sorted = 1, int blocked = 0,
                              int levels = 0, int *levelPtr = NULL, int *levelSet = NULL,
                              int parts = 0, int *parPtr = NULL, int *partition = NULL,
                              int s_level_no = -1, int *s_level_ptr = NULL,
                              int *s_level_set = NULL,
                              int chunk = 0, double *ws = NULL, double *ws_zerod = NULL);

}

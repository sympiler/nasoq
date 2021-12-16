//
// Created by kazem on 7/27/17.
//

#ifndef CHOLOPENMP_SPARSEUTILS_H
#define CHOLOPENMP_SPARSEUTILS_H

#include <cstdint>
#include <vector>

#include "nasoq/common/def.h"

namespace nasoq {
/* Return the number of entries in a sparse matrix.
 *
 * workspace: none
 * integer overflow cannot occur, since the matrix is already allocated.
 */

 long int getNNZ
   (
     /* ---innerPartsSize[l]- input ---- */
     int ncol,
     int *Ap,
     int *Anz,
     int packed,
     /* --------------- */
     int &status
   );


/* ========================================================================== */
/* === cholmod_add_size_t =================================================== */
/* ========================================================================== */

/* Safely compute a+b, and check for integer overflow.  If overflow occurs,
 * return 0 and set OK to FALSE.  Also return 0 if OK is FALSE on input. */

 size_t add_size_t(size_t a, size_t b, int *ok);


/* ========================================================================== */
/* === cholmod_mult_size_t ================================================== */
/* ========================================================================== */

/* Safely compute a*k, where k should be small, and check for integer overflow.
 * If overflow occurs, return 0 and set OK to FALSE.  Also return 0 if OK is
 * FALSE on input. */

 size_t mult_size_t(size_t a, size_t k, int *ok);

 int makeUnique(int *node2Par,
                std::vector<int> &list,
                int n, bool *ws);

/*
 * b = alpha*a+b
 */
 void add_vec(int n, double *a, double alpha, double *b);

/*
 * converting sparse csc to dense matrix col-wise
 */
 void sparse2dense(CSC *A, double *D);

/*
 * transpose a dense matrix
 */
 void dense_traspose(int r, int c, double *a_in, double *a_out);

/// Used in QP
/// \param A
/// \param a
/// \param rec_norm
/// \return
 int compute_recieporical_length(CSC *A, double *a,
                                 double *rec_norm);
}
#endif //CHOLOPENMP_SPARSEUTILS_H

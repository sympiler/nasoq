//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_ETREE_H
#define CHOLOPENMP_ETREE_H

#include <cstdint>
#include <cstdio>

#include "nasoq/common/def.h"

namespace nasoq {
/* ========================================================================== */
/* === update_etree ========================================================= */
/* ========================================================================== */

 void update_etree
   (
     /* inputs, not modified */
     int k,  /* process the edge (k,i) in the input graph */
     int i,
     /* inputs, modified on output */
     int Parent[], /* Parent [t] = p if p is the parent of t */
     int Ancestor[] /* Ancestor [t] is the ancestor of node t in the
			   partially-constructed etree */
   );


/* ========================================================================== */
/* === cholmod_etree ======================================================== */
/* ========================================================================== */

/* Find the elimination tree of A or A'*A */

 int etreeC
   (
     /* ---- input ---- */
     CSC *A,
     /* ---- output --- */
     int *Parent, /* size ncol.  Parent [j] = p if p is the parent of j */
     /* --------------- */
     //cholmod_common *Common
     int &status
   );

/* ========================================================================== */
/* === cholmod_etree ======================================================== */
/* ========================================================================== */

/* Find the elimination tree of A or A'*A */

 int etree1
   (
     /* ---- input ---- */
     // cholmod_sparse *A,
     int ncol,
     int nrow,
     int *Ap,
     int *Ai,
     int *Anz,
     int stype,
     int packed,
     /* ---- output --- */
     int *Parent /* size ncol.  Parent [j] = p if p is the parent of j */
     /* --------------- */
     //cholmod_common *Common
   );


 int *etree(int n, int *Ap, int *Ai, int ata);

}

#endif //CHOLOPENMP_ETREE_H

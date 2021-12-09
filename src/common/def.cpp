//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/def.h"

namespace nasoq {

 void allocateAC(CSC *A, int nrow, int nnz, int sytpe, int sw) {
  if (sw) {
   A->nrow = A->ncol = nrow;
   A->nzmax = nnz;
   A->stype = sytpe;
   A->xtype = CHOLMOD_REAL;//TODO removed later
   A->packed = TRUE; // Always
   A->p = new int[nrow + 1]();
   A->i = new int[nnz]();
   A->x = new double[nnz]();
   A->nz = NULL;

  } else {
   delete[]A->p;
   delete[]A->i;
   delete[]A->x;
   delete A;
  }


 }

 void allocateLC(BCSC *L, int sw) {
  int sNo = L->nsuper;
  if (sw) {
   L->super = new int[sNo + 1]();
//  L->col2Sup = new int[L->n](); //TODO HACK
   L->p = new size_t[L->ssize + 1]();
   L->pi = new size_t[sNo + 1]();
   L->i_ptr = new size_t[L->xsize + 1](); // index pointers
   L->i = new int[sNo + 1]();//Nothing for now
   // L->px = new int[L->xsize]();
   L->s = new int[L->xsize]();//Index values
   //L->sParent = new int[sNo](); //TODO HACK
   // L->x = new double[L->xsize]();
   L->is_ll = TRUE;
   L->xtype = CHOLMOD_REAL;
   L->is_super = TRUE;

  } else {
   delete[]L->super;
   delete[]L->col2Sup;
   delete[]L->p;
   delete[]L->pi;
   delete[]L->i;
   delete[]L->i_ptr;
   // delete []L->px;
   delete[]L->s;
   delete[]L->sParent;
   delete[]L->Parent;
   delete[]L->IPerm;
   delete[]L->Perm;
//  delete []L->x;
  }


 }
}
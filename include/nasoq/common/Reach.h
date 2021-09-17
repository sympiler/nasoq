//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_REACH_H
#define CHOLOPENMP_REACH_H

namespace nasoq {
/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
 int reach(int n, int *Gp, int *Gi, int *Bp, int *Bi, int k, int *xi, const int *pinv);

 int reach_sn(int n, int *Gp, int *Gi, int *Bp, int *Bi,
              int k, int *PBset, const int *pinv, int sn,
              int *col2sup, double &analysis);

 int reach_col(int n, int *Gp, int *Gi, int *Bp, int *Bi, int k, int *xi, const int *pinv);

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
 int ereach(int n, int *Ap, int *Ai, int k, const int *parent,
            int *s, int *w);

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
 int ereach_sn(int n, int *Ap, int *Ai, int col1, int col2, int *col2sup,
               const int *parent, int *s, int *w);
}
#endif //CHOLOPENMP_REACH_H

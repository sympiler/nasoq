//
// The modified version of CHOLMOD code
//

#ifndef CHOLOPENMP_POSTORDER_H
#define CHOLOPENMP_POSTORDER_H

#include <cstddef>
#include <cstdint>

namespace nasoq {
 int postOrderC /* return # of nodes postordered */
   (
     /* ---- input ---- */
     int *Parent, /* size n. Parent [j] = p if p is the parent of j */
     size_t n,
     int *Weight, /* size n, optional. Weight [j] is weight of node j */
     /* ---- output --- */
     int *Post,  /* size n. Post [k] = j is kth in postordered tree */
     /* --------------- */
     // cholmod_common *Common
     int status
   );


//CSPARSE
 int *postOrder(const int *parent, int n);

}
#endif //CHOLOPENMP_POSTORDER_H

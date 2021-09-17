//
// Created by kazem on 1/27/18.
//

#ifndef CHOLOPENMP_PERFORMANCEMODEL_H
#define CHOLOPENMP_PERFORMANCEMODEL_H

#include <cstddef>

namespace nasoq {
#define OPS_SCAL(i)     (i)
#define OPS_COPY(i)     (i)
/*
#define OPS_AXPY(i)     (i)
#define OPS_SWAP(i)     (i)
#define OPS_GEMV(i,j)   ((i)*(j))
#define OPS_TRSV(i,j)   (((i)*(i-1)*(j))/2.)
#define OPS_GER(i,j)    ((i)*(j))
*/
#define OPS_GEMM(i, j, k) (2.0*(i)*(j)*(k))
#define OPS_TRSM(i, j)   ((i)*(i-1)*(j))
#define OPS_PPF(i)      (0.3333*(i)*(i)*(i)+2.*(i)*(i))
#define OPS_GEAM(i, j)   ((i)*(j))


/**GEMM**/
#define GEMM_A  2.429169e-10
#define GEMM_B  2.724804e-10
#define GEMM_C  1.328900e-09
#define GEMM_D  1.148989e-07
#define GEMM_E  -2.704179e-10
#define GEMM_F  1.216278e-06
#define PERF_GEMM(i, j, k) (GEMM_A*(double)(i)*(double)(j)*(double)(k)+GEMM_B*(double)(i)*(double)(j)+GEMM_C*(double)(j)*(double)(k)+GEMM_D*(double)(i)+GEMM_E*(double)(j)+GEMM_F)


/**GEAM**/
#define GEAM_A  1.358111e-09
#define GEAM_B -4.416379e-09
#define GEAM_C 2.270780e-08
#define GEAM_D  -3.335563e-07
#define PERF_GEAM(i, j)   (GEAM_A*(double)(i)*(double)(j)+GEAM_B*(double)(i)+GEAM_C*(double)(j)+GEAM_D)

/**TRSM**/
#define TRSM_A 2.626177e-10
#define TRSM_B 3.976198e-08
#define TRSM_C 3.255168e-06
#define PERF_TRSM(i, j)   (TRSM_A*(double)(i)*(double)(i)*(double)(j)+TRSM_B*(double)(i)+TRSM_C)

/**POF**/
#define POF_A 2.439599e-11
#define POF_B 1.707006e-08
#define POF_C -1.469893e-07
#define POF_D 4.071507e-07
#define PERF_POF(i) (POF_A*(double)(i)*(double)(i)*(double)(i)+POF_B*(double)(i)*(double)(i)+POF_C*(double)(i)+POF_D)

/**PPF**/
#define PPF_A 2.439599e-11
#define PPF_B 1.707006e-08
#define PPF_C -1.469893e-07
#define PPF_D 4.071507e-07
#define PERF_PPF(i) (PPF_A*(double)(i)*(double)(i)*(double)(i)+PPF_B*(double)(i)*(double)(i)+PPF_C*(double)(i)+PPF_D)

/**SCAL**/
#define SCAL_A 4.371793e-10
#define SCAL_B 2.052399e-07
#define PERF_SCAL(i) (SCAL_A*(double)(i)+SCAL_B)

/**COPY**/
#define COPY_A 9.177969e-10
#define COPY_B 2.266129e-07
#define PERF_COPY(i) (COPY_A*(double)(i)+COPY_B)

/**AXPY**/
#define AXPY_A 4.620143e-10
#define AXPY_B 2.101008e-07
#define PERF_AXPY(i) (AXPY_A*(double)(i)+AXPY_B)

/**GEMV**/
#define GEMV_A  6.192657e-10
#define GEMV_B -2.884799e-09
#define GEMV_C 7.594831e-10
#define GEMV_D  3.575035e-07
#define PERF_GEMV(i, j)   (GEMV_A*(double)(i)*(double)(j)+GEMV_B*(double)(i)+GEMV_C*(double)(j)+GEMV_D)

/**TRSV**/
#define TRSV_A 3.224536e-10
#define TRSV_B 1.709178e-08
#define TRSV_C 1.947268e-07
#define PERF_TRSV(i) (TRSV_A*(double)(i)*(double)(i)+TRSV_B*(double)(i)+TRSV_C)


 double computeCostColFact(double colNo, double len);

/*
 * number of FLOPs
 */
 double computeCostColFact2(double colNo, double len);

/*
 * The number of col as for the cost
 */
 double computeCostColFact3(double colNo, double len);

/*
 * The number of NNZ as for the cost
 */
 double computeCostColFact4(double colNo, double len);

 double computeContribCost(double dim1, double midDim, double dim3);

 double computeContribCost2(double dim1, double midDim, double dim3);

/*
 * using the number of col only
 */
 double computeContribCost3(double dim1, double midDim, double dim3);


/*
 * using the number of NNZ only
 */
 double computeContribCost4(double dim1, double midDim, double dim3);

/*
 * For CSC Matrix
 */
 double computeCostperCol(int n, int colNo, int *eTree,
                          int *cT, int *rT, int *xi, int &top);


 double computeCostperBlock(int n, int curCol, int nxtCol, int *eTree,
                            int *cT, int *rT, int *col2Sup, int *blockSet,
                            int *lR, size_t *Li_ptr, int *xi);

 double computeComCost(int n, int *c, int *r, double *values,
                       size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                       int *blockSet, int supNo, double *timing,
                       int *aTree, int *cT, int *rT, int *col2Sup,
                       int nLevels, int *levelPtr, int *levelSet,
                       int nPar, int *parPtr, int *partition,
   //Outputs
                       long &totalIntraContribNNZs, long &totalInterCotribNNZs,
                       int &numberOfIntraCore, int &numberOfInterCore,
                       long &numberOfEdgeCuts, double &avgNNZperThread
 );
}
#endif //CHOLOPENMP_PERFORMANCEMODEL_H

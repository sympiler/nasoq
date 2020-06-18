//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_UTIL_H
#define CHOLOPENMP_UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
//#include <getopt.h>
#include <map>
#include "mkl.h"
#include "def.h"
#include "transpose_unsym.h"

namespace nasoq {
/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
 bool readMatrix_old(std::string fName, int &n, int &NNZ, int *&col,
                     int *&row, double *&val) {
  /*This function reads the input matrix from "fName" file and
   * allocate memory for matrix A, L and U.
   * - The input file is a coordinate version and e
   * ach row of the file shows (col, row, nnz)
   * - The matrices are zero-indexed
   */
  std::ifstream inFile;
  inFile.open(fName);
  inFile >> n;
  inFile >> n;
  inFile >> NNZ;
  int factorSize = (n * n) / 2;//Worst case assumption
  if (n <= 0 || NNZ <= 0)
   return false;
  col = new int[n + 1]();
  // colL = new int[n + 1]; colU = new int[n + 1];
  row = new int[NNZ];
  // rowL = new int[factorSize]; rowU = new int[factorSize];
  val = new double[NNZ];
  // valL = new double[factorSize]; valU = new double[factorSize];
  if (!val || !col || !row)
   return false;
  //Initializing the result vector
  int y, x, colCnt = 0, nnzCnt = 0;
  double value;

  col[0] = 0;
  for (int i = 0; nnzCnt < NNZ;) {//Reading from file row by row
   inFile >> x;
   x--;
   inFile >> y;
   y--;//zero indexing
   inFile >> value;
   if (y > n)
    return false;
   if (y == i) {
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    colCnt++;
    nnzCnt++;
   } else {//New col
    col[i + 1] = col[i] + colCnt;
    i++;//next iteration
    colCnt = 1;
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    nnzCnt++;
   }

  }
  col[n] = col[n - 1] + colCnt;//last col

  return true;
 }

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
 bool readMatrix(std::string fName, size_t &n, size_t &NNZ, int *&col,//FIXME change col type to size_t
                 int *&row, double *&val) {
  /*This function reads the input matrix from "fName" file and
   * allocate memory for matrix A, L and U.
   * - The input file is a coordinate version and e
   * ach row of the file shows (col, row, nnz)
   * - The matrices are zero-indexed
   */

  std::ifstream inFile;
  inFile.open(fName);
  std::string line, banner, mtx, crd, arith, sym;
  /*  File format:
   *    %%MatrixMarket matrix coordinate real general/symmetric/...
   *    % ...
   *    % (optional comments)
   *    % ...
   *    #rows    #non-zero
   *    Triplet in the rest of lines: row    col    value
   */
  std::getline(inFile, line);
  for (unsigned i = 0; i < line.length(); line[i] = tolower(line[i]), i++);
  std::istringstream iss(line);
  if (!(iss >> banner >> mtx >> crd >> arith >> sym)) {
   std::cout << "Invalid header (first line does not contain 5 tokens)\n";
   return false;
  }

  if (banner.compare("%%matrixmarket")) {
   std::cout << "Invalid header (first token is not \"%%%%MatrixMarket\")\n";
   return false;
  }
  if (mtx.compare("matrix")) {
   std::cout << "Not a matrix; this driver cannot handle that.\"\n";
   return false;
  }
  if (crd.compare("coordinate")) {
   std::cout << "Not in coordinate format; this driver cannot handle that.\"\n";
   return false;
  }
  if (arith.compare("real") && arith.compare("integer")) {
   if (!arith.compare("complex")) {
    std::cout << "Complex matrix; use zreadMM instead!\n";
    return false;
   } else if (!arith.compare("pattern")) {
    std::cout << "Pattern matrix; values are needed!\n";
    return false;
   } else {
    std::cout << "Unknown arithmetic\n";
    return false;
   }
  }
  while (!line.compare(0, 1, "%")) {
   std::getline(inFile, line);
  }
  std::istringstream issDim(line);
  if (!(issDim >> n >> n >> NNZ)) {
   std::cout << "The matrix dimension is missing\n";
   return false;
  }
  if (n <= 0 || NNZ <= 0)
   return false;
  col = new int[n + 1]();
  // colL = new int[n + 1]; colU = new int[n + 1];
  row = new int[NNZ];
  // rowL = new int[factorSize]; rowU = new int[factorSize];
  val = new double[NNZ];
  // valL = new double[factorSize]; valU = new double[factorSize];
  if (!val || !col || !row)
   return false;
  //Initializing the result vector
  int y, x, colCnt = 0, nnzCnt = 0;
  double value;

  col[0] = 0;
  int sw = 1;
  for (int i = 0; nnzCnt < NNZ;) {//Reading from file row by row
   inFile >> x;
   x--;
   inFile >> y;
   y--;//zero indexing
   inFile >> value;
   if (i == 0 && sw) {//matrix starts with empty cols
    if (y != 0) {
     i = y;
    }
    sw = 0;
   }
   if (y > n)
    return false;
   if (y == i) {
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    colCnt++;
    nnzCnt++;
   } else {//New col
    int nnz_c = col[i] + colCnt;
    col[i + 1] = nnz_c;
    i++;//next iteration
    for (int j = i; j < y; ++j) {//empty cols
     col[j + 1] = nnz_c;
     i++;
    }
    colCnt = 1;
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    nnzCnt++;
   }
  }
  //col[y+1] = col[y] + colCnt;
  if (y == n - 1) {//each column has something in it
   col[n] = col[n - 1] + colCnt;//last col
  } else { // we have some empty columns
   int nnzLast = col[y] + colCnt;
   col[y + 1] = nnzLast;
   y++;
   assert(nnzLast == NNZ);
   for (int i = y + 1; i < n + 1; ++i) {
    col[i] = nnzLast;
   }
  }

  return true;
 }


/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
 bool readMatrix_rect(std::string fName, size_t &n_row, size_t &n_col,
                      size_t &NNZ, int *&col,//FIXME change col type to size_t
                      int *&row, double *&val) {
  /*This function reads the input matrix from "fName" file and
   * allocate memory for matrix A, L and U.
   * - The input file is a coordinate version and e
   * ach row of the file shows (col, row, nnz)
   * - The matrices are zero-indexed
   */

  std::ifstream inFile;
  inFile.open(fName);
  std::string line, banner, mtx, crd, arith, sym;
  /*  File format:
   *    %%MatrixMarket matrix coordinate real general/symmetric/...
   *    % ...
   *    % (optional comments)
   *    % ...
   *    #rows    #non-zero
   *    Triplet in the rest of lines: row    col    value
   */
  std::getline(inFile, line);
  for (unsigned i = 0; i < line.length(); line[i] = tolower(line[i]), i++);
  std::istringstream iss(line);
  if (!(iss >> banner >> mtx >> crd >> arith >> sym)) {
   std::cout << "Invalid header (first line does not contain 5 tokens)\n";
   return false;
  }

  if (banner.compare("%%matrixmarket")) {
   std::cout << "Invalid header (first token is not \"%%%%MatrixMarket\")\n";
   return false;
  }
  if (mtx.compare("matrix")) {
   std::cout << "Not a matrix; this driver cannot handle that.\"\n";
   return false;
  }
  if (crd.compare("coordinate")) {
   std::cout << "Not in coordinate format; this driver cannot handle that.\"\n";
   return false;
  }
  if (arith.compare("real") && arith.compare("integer")) {
   if (!arith.compare("complex")) {
    std::cout << "Complex matrix; use zreadMM instead!\n";
    return false;
   } else if (!arith.compare("pattern")) {
    std::cout << "Pattern matrix; values are needed!\n";
    return false;
   } else {
    std::cout << "Unknown arithmetic\n";
    return false;
   }
  }
  while (!line.compare(0, 1, "%")) {
   std::getline(inFile, line);
  }
  std::istringstream issDim(line);
  if (!(issDim >> n_row >> n_col >> NNZ)) {
   std::cout << "The matrix dimension is missing\n";
   return false;
  }
  if (n_col < 0 || NNZ < 0)
   return false;
  if (n_col == 0 || NNZ == 0 || n_row == 0) {
   col = NULL;
   row = NULL;
   val = NULL;
   return true;
  }
  col = new int[n_col + 1]();
  // colL = new int[n + 1]; colU = new int[n + 1];
  row = new int[NNZ];
  // rowL = new int[factorSize]; rowU = new int[factorSize];
  val = new double[NNZ];
  // valL = new double[factorSize]; valU = new double[factorSize];
  if (!val || !col || !row)
   return false;
  //Initializing the result vector
  int y, x, colCnt = 0, nnzCnt = 0;
  double value;

  col[0] = 0;
  for (int i = 0; nnzCnt < NNZ;) {//Reading from file row by row
   inFile >> x;
   x--;
   inFile >> y;
   y--;//zero indexing
   inFile >> value;
   if (y > n_col)
    return false;
   if (y == i) {
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    colCnt++;
    nnzCnt++;
   } else if (i + 1 == y) {//New col
    col[i + 1] = col[i] + colCnt;
    i++;//next iteration
    colCnt = 1;
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    nnzCnt++;
   } else {  // new non-consecutive col, we have gap y > i+1
    col[i + 1] = col[i] + colCnt;
    i++;
    for (int j = i + 1; j <= y; ++j, ++i) { //fill up the col gap
     col[j] = col[i];
    }
    val[nnzCnt] = value;
    row[nnzCnt] = x;
    nnzCnt++;
    colCnt = 1;
    // now y == i

   }
  }
  //col[y+1] = col[y] + colCnt;
  if (y == n_col - 1) {//each column has something in it
   col[n_col] = col[n_col - 1] + colCnt;//last col
  } else { // we have some empty columns
   int nnzLast = col[y] + colCnt;
   col[y + 1] = nnzLast;
   y++;
   assert(nnzLast == NNZ);
   for (int i = y + 1; i < n_col + 1; ++i) {
    col[i] = nnzLast;
   }
  }

  return true;
 }

/*
 *
 */
 void read_dense(std::string fname, int &n_row, int &n_col, double *&mat) {
  std::ifstream in_file;
  in_file.open(fname);
  in_file >> n_row;
  in_file >> n_col;
  mat = new double[n_col * n_row];
  for (int i = 0; i < n_col * n_row; i++) {//writing from file row by row
   in_file >> mat[i];
  }
 }

/*
 *
 */
 int write_dense(std::string fname, int n_row, int n_col, double *mat) {
  if (n_row <= 0 || n_col <= 0)
   return 0;
  double value;
  std::ofstream out_file;
  out_file.precision(30);
  out_file.open(fname);
  out_file << n_row << " ";
  out_file << n_col << " ";
  out_file << "\n";
  for (int i = 0; i < n_col * n_row; i++) {//writing from file row by row
   value = mat[i];
   out_file << value << " ";
  }
  return 1;
 }

/*
 * expand the empty columns by one zero
 * e.g., in indefinite matrices where empty columns exist
 * returns false if the matrix is not expanded.
 */
 bool expandMatrix(size_t n, size_t NNZ, const int *col,
                   const int *row, const double *val,
                   size_t &newNNZ, int *&newCol,//FIXME change col type to size_t
                   int *&newRow, double *&newVal,
                   double insDiag = 0) {
  int emptycols = 0;
  //int cur_c=0;//case for empty col at the begining
  for (int i = 0; i < n; ++i) {
   if (col[i + 1] - col[i] == 0)
    emptycols++;
/*   if(row[col[cur_c]] == i){
    cur_c++;
   }*/
  }
  if (emptycols == 0)
   return false;
  newNNZ = NNZ + emptycols;
  newCol = new int[n + 1]();
  newRow = new int[newNNZ];
  newVal = new double[newNNZ];
  int curNNZ = 0;
  int insertedNNZ = 0;
  newCol[0] = 0;
  for (int i = 0; i < n; ++i) {
   if (col[i + 1] - col[i] == 0) {
    newVal[col[i] + insertedNNZ] = insDiag;//inserted val
    newRow[col[i] + insertedNNZ] = i; // inseted diagonal row
    insertedNNZ++;
   } else {
    for (int j = col[i]; j < col[i + 1]; ++j) {
     newVal[j + insertedNNZ] = val[j];
     newRow[j + insertedNNZ] = row[j];
    }
   }
   newCol[i + 1] = col[i + 1] + insertedNNZ;
  }
  assert(newCol[n] == newNNZ);

  return true;
 }

/*
 * reading a ordering from PasTiX stored format.
 */
 int read_vector(std::string fName, size_t n, double *perm) {
  std::ifstream inFile;
  inFile.open(fName);
  if (n <= 0)
   return 0;

  std::string line, banner, mtx, crd, arith, sym;
  /*  File format:
   *    %%MatrixMarket matrix coordinate real general/symmetric/...
   *    % ...
   *    % (optional comments)
   *    % ...
   *    #rows    #non-zero
   *    Triplet in the rest of lines: row    col    value
   */
  std::getline(inFile, line);
  for (unsigned i = 0; i < line.length(); line[i] = tolower(line[i]), i++);
  std::istringstream iss(line);
  if (!(iss >> banner >> mtx >> crd >> arith >> sym)) {
   //std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
   //return false;

   //FIXME: will fix this part once all vectors are mtx
   //Initializing the result vector
   inFile.clear();
   inFile.seekg(0, std::ios::beg);
   double value;
   for (int i = 0; i < n; i++) {//Reading from file row by row
    inFile >> value;
    perm[i] = value;
   }
   return 1;
  }

  if (banner.compare("%%matrixmarket")) {
   std::cout << "Invalid header (first token is not \"%%%%MatrixMarket\")\n";
   return false;
  }
  if (mtx.compare("matrix")) {
   std::cout << "Not a matrix; this driver cannot handle that.\"\n";
   return false;
  }
  if (crd.compare("array")) {
   std::cout << "Not in coordinate format; this driver cannot handle that.\"\n";
   return false;
  }
  if (arith.compare("real") && arith.compare("integer")) {
   if (!arith.compare("complex")) {
    std::cout << "Complex matrix; use zreadMM instead!\n";
    return false;
   } else if (!arith.compare("pattern")) {
    std::cout << "Pattern matrix; values are needed!\n";
    return false;
   } else {
    std::cout << "Unknown arithmetic\n";
    return false;
   }
  }
  while (!line.compare(0, 1, "%")) {
   std::getline(inFile, line);
  }
  int n_row, n_col;
  std::istringstream issDim(line);
  if (!(issDim >> n_row >> n_col)) {
   std::cout << "The matrix dimension is missing\n";
   return false;
  }
  assert(n_row == n);
  //Initializing the result vector
  double value;

  for (int i = 0; i < n; i++) {//Reading from file row by row
   inFile >> value;
   perm[i] = value;
  }
  return 1;
 }

/*
 * writing a vector to a file.
 */
 int write_vector(std::string fName, size_t n,
                  double *vec_vals, std::string header = "",
                  int prec = 30) {
  if (n <= 0)
   return 0;
  double value;
  std::ofstream out_file;
  out_file.precision(prec);
  out_file.open(fName);
  out_file << header;
  for (int i = 0; i < n; i++) {//writing from file row by row
   value = vec_vals[i];
   out_file << value;
   if (i != n - 1)
    out_file << "\n";
  }
  out_file.close();
  return 1;
 }

/*
 * reading a ordering from PasTiX stored format.
 */
 bool readOrdering(std::string fName, size_t n, size_t *perm) {
  std::ifstream inFile;
  inFile.open(fName);
  std::string line, banner, mtx, crd, arith, sym;
  size_t dontKnow;
  /*  File format:
   *    # 0
   *    #??    #rows
   *    values values ...
   */
  std::getline(inFile, line);
  std::istringstream iss(line);

  while (!line.compare(0, 1, "%")) {
   std::getline(inFile, line);
  }


  if (!(iss >> n)) {
   std::cout << "The matrix dimension is missing\n";
   return false;
  }

  if (n <= 0)
   return false;
  //Initializing the result vector
  int y, x, colCnt = 0, nnzCnt = 0;
  double value;

  for (int i = 0; colCnt < n; i++, colCnt++) {//Reading from file row by row
   inFile >> x;
   perm[i] = x;
  }
  assert(colCnt == n);
  return true;
 }

 bool enableColdCache(int n, std::ifstream &f) {
  /*
   * n specifies the size of data for double computation. It depends
   * on the cache size
   */
  //TODO check file during read
  assert(!f.fail());
  double curVal;
  double **waste = new double *[n];
  for (int i = 0; i < n; ++i) {
   waste[i] = new double[n];
  }
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
    f >> curVal;
    waste[i][j] = curVal;
   }
  }
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
     waste[i][j] += waste[i][k] * waste[k][j];
    }
   }
  }
  for (int i = 0; i < n; ++i) {
   delete waste[i];
  }
  delete[] waste;
  return true;
 }

/*
 * For linear solver using A and A'
 */

 void rhsInit_linearSolver(int n, int *Ap, int *Ai, double *Ax, //A
                           int *Bp, int *Bi, double *Bx, //AT
                           double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (MKL_INT cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
    b[Ai[cc]] += Ax[cc];
   }
   for (int i = Bp[c] + 1; i < Bp[c + 1]; ++i) {
    b[Bi[i]] += Bx[i];
   }
  }
 }

/*
 * Assuming we have only lower part of a symmetric matrix
 * return false if rhs is zero
 */
 bool rhsInit_linearSolver_onlylower(int n, int *Ap, int *Ai, double *Ax, //A
                                     double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
    b[Ai[cc]] += Ax[cc];
    if (Ai[cc] != c) { // upper part
     b[c] += Ax[cc];
    }
   }
  }
  return true;
 }

/*
 * For triangular solve
 */

 void rhsInit(int n, MKL_INT *Ap, MKL_INT *Ai, double *Ax, double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (MKL_INT cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
    b[Ai[cc]] += Ax[cc];
   }
  }
 }

/*
 * RHS initilization for blocked
 */

 void rhsInitBlocked(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, double *Ax, double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (int cc = Ap[c], j = 0; cc < Ap[c + 1]; ++cc, ++j) {
    size_t curRow = Ai[AiP[c] + j];
    b[curRow] += Ax[cc];
   }
  }
 }

/*
 * RHS initilization for blocked L'
 */

 void rhsInitBlockedLT(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, double *Ax, double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (int cc = Ap[c], j = 0; cc < Ap[c + 1]; ++cc, ++j) {
    b[c] += Ax[cc]; // will be the row of L'
   }
  }
 }

/*
 * Making full symmetric matrix from lower and upper
 */

 void make_full(int ncol, int nnz, int *Ap, int *Ai, double *Ax, //A
                int *ATransp, int *ATransi, double *ATransx,
                int &nnzFull, int *&AFullp, int *&AFulli, double *&AFullx) {
  //Making a full symmetric matrix with
  //both upper and lower parts
  nnzFull = nnz * 2 - ncol;
  AFullp = new int[ncol + 1]();
  AFulli = new int[nnzFull]();
  AFullx = new double[nnzFull]();
  int ncolIDXT = ncol;
  AFullp[0] = 0;
  for (int i = 0; i < ncol; ++i) {
   int nnzOfCurCol = ATransp[i + 1] - ATransp[i] - 1;
   nnzOfCurCol += Ap[i + 1] - Ap[i];
   AFullp[i + 1] = (long int) AFullp[i] + nnzOfCurCol;
   //copying Upper part, ignoring diagonal since it is in L already
   int base = AFullp[i];
   for (int j = ATransp[i], k = 0; j < ATransp[i + 1] - 1; ++j, ++k) {
    AFulli[base + k] = (long int) ATransi[j];
    AFullx[base + k] = ATransx[j];
   }
   //copying L part
   base += ATransp[i + 1] - ATransp[i] - 1;
   for (int j = Ap[i], k = 0; j < Ap[i + 1]; ++j, ++k) {
    AFulli[base + k] = (long int) Ai[j];
    AFullx[base + k] = Ax[j];
   }
  }
 }

/*
 * Testing lower triangular solve to make sure all x elements are ONE.
 * Warning: This works only for when RHS does not have zero
 */

 int testTriangular(size_t n, const double *x, double epsilon = 1e-3) {//Testing
  int test = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(1 - x[i]) < epsilon) {
    test++;
   } /*else{
   std::cout<<i<<" : "<<1-x[i]<<";";
  }*/
   //else
   // cout<<i<<";";
  }
  if (n - test > 0) {
   return false;
  }
  return true;
 }



/*
 * converting CSC to CSR
 */  //TODO
 void csc2csr(int n, int *Ap, int *Ar, double *Ax,
              int *row_ptr, int *col_idx, double *Bx) {
  int *row_len = new int[n]();
  int *finger = new int[n]();
  for (int i = 0; i < Ap[i + 1]; ++i) {
   row_len[Ar[i]]++;
  }
  col_idx[0] = 0;
  for (int l = 1; l < n + 1; ++l) {
   col_idx[l] = col_idx[l - 1] + row_len[l];
  }
  for (int k = 0; k < n; ++k) {
   for (int i = Ap[i]; i < Ap[i + 1]; ++i) {
    //Bx[col_idx[]]
   }
  }
 }


/*
 * Converting BCSC to CSC
 */
 int bcsc2csc(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col, double *Ax,
   //out
   MKL_INT *Cp, MKL_INT *Ci, double *Cx
 ) {
  size_t actualNNZ = 0;
  Cp[0] = 0;
  for (int i = 0; i < nBlocks; ++i) {
   int curCol = sup2col[i];
   int nxtCol = sup2col[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);

   for (int j = curCol; j < nxtCol; ++j) {
    for (MKL_INT k = Ap[j] + (j - curCol), kk = AiP[curCol] + (j - curCol);
         k < Ap[j + 1]; ++k, ++kk) {
     Cx[actualNNZ] = Ax[k];
     assert(Ai[kk] < n);
     Ci[actualNNZ] = Ai[kk];
     actualNNZ++;
    }
    Cp[j + 1] = actualNNZ;
   }
  }
  return 1;
 }


/*
 * Converting BCSC to CSC, removing all zeros, assuming
 * fill-ins are initially zero.
 * Warning: could be problematic in some exceptional cases where,
 * exact numerical cancellation happens.
 */
 int bcsc2csc_aggressive(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col, double *Ax,
   //out
   //MKL_INT *Cp, MKL_INT *Ci, double *Cx
   size_t *Cp, int *Ci, double *Cx
 ) {
  size_t actualNNZ = 0;
  Cp[0] = 0;
  for (int i = 0; i < nBlocks; ++i) {
   int curCol = sup2col[i];
   int nxtCol = sup2col[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);
/*  if(i==660)
   printf("here");*/
   for (int j = curCol; j < nxtCol; ++j) {
    for (MKL_INT k = Ap[j] + (j - curCol), kk = AiP[curCol] + (j - curCol);
         k < Ap[j + 1]; ++k, ++kk) {
     if (Ax[k] != 0) {
      Cx[actualNNZ] = Ax[k];
      assert(Ai[kk] < n);
      Ci[actualNNZ] = Ai[kk];
      actualNNZ++;
     }
    }
    Cp[j + 1] = actualNNZ;
   }
  }
  return 1;
 }


 int check_row_idx_l(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col
 ) {
  size_t actualNNZ = 0;
  for (int i = 0; i < nBlocks; ++i) {
   int curCol = sup2col[i];
   int nxtCol = sup2col[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);
   int prev_idx = curCol - 1;
   for (int l = AiP[curCol]; l < AiP[curCol + 1]; ++l) {
    if (Ai[l] <= prev_idx) {
     std::cout << "Error in " << curCol << "\n";
     assert(false);
    } else {
     prev_idx = Ai[l];
    }
   }
  }
  return 1;
 }

 int bcsc2csc_aggressive_int(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col, double *Ax,
   //out
   //MKL_INT *Cp, MKL_INT *Ci, double *Cx
   int *Cp, int *Ci, double *Cx
 ) {
  size_t actualNNZ = 0;
  Cp[0] = 0;
  for (int i = 0; i < nBlocks; ++i) {
   int curCol = sup2col[i];
   int nxtCol = sup2col[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);

   for (int j = curCol; j < nxtCol; ++j) {
    for (MKL_INT k = Ap[j] + (j - curCol), kk = AiP[curCol] + (j - curCol);
         k < Ap[j + 1]; ++k, ++kk) {
     if (Ax[k] != 0) {
      Cx[actualNNZ] = Ax[k];
      assert(Ai[kk] < n);
      Ci[actualNNZ] = Ai[kk];
      actualNNZ++;
     }
    }
    Cp[j + 1] = actualNNZ;
   }
  }
  return 1;
 }

 void swapWSet(std::vector<std::vector<int>> &List, int i, int j) {
  std::vector<int> tmp;
  if (List[i].size() == 0 || List[j].size() == 0) {
   std::cout << "wrong inputs\n";
   return;
  }
  tmp.insert(tmp.begin(), List[i].begin(), List[i].end());
  assert(tmp.size() == List[i].size());
  List[i].erase(List[i].begin(), List[i].end());
  assert(List[i].size() == 0);
  assert(tmp.size() > 0);
  List[i].insert(List[i].begin(), List[j].begin(), List[j].end());

  List[j].erase(List[j].begin(), List[j].end());
  List[j].insert(List[j].begin(), tmp.begin(), tmp.end());
 }


/*
 * Check whether the matrix is stored in full
 * return the lower half if it is full
 * Makes sense only for symmetric matrices.
 */
 CSC *computeLowerTriangular(CSC *A) {
  CSC *lowerHalf;
  size_t lNNZcnt = 0;
  for (int i = 0; i < A->ncol; ++i) {
   for (int j = A->p[i]; j < A->p[i + 1]; ++j) {
    if (A->i[j] >= i)
     lNNZcnt++;
   }
  }
  //lower is already stored
  if (lNNZcnt == A->nzmax || lNNZcnt == 0)
   return NULL;
  lowerHalf = new CSC;
  allocateAC(lowerHalf, A->ncol, lNNZcnt, 0, true);
  lowerHalf->ncol = lowerHalf->ncol = A->ncol;
  size_t curNNZ = 0;
  lowerHalf->p[0] = 0;
  //copying into new space
  for (int i = 0; i < A->ncol; ++i) {
   for (int j = A->p[i]; j < A->p[i + 1]; ++j) {
    if (A->i[j] >= i) {
     lowerHalf->i[curNNZ] = i;
     lowerHalf->x[curNNZ] = A->x[j];
     curNNZ++;
    }
   }
   lowerHalf->p[i + 1] = curNNZ;
  }
  lowerHalf->nzmax = lNNZcnt;
  return lowerHalf;
 }


 void print_csc(std::string beg, size_t n, int *Ap, int *Ai, double *Ax) {

  std::cout << beg;
  int nnz = n > 0 ? Ap[n] : 0;
  std::cout << n << " " << n << " " << nnz << "\n";
  for (int i = 0; i < n; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    std::cout << Ai[j] + 1 << " " << i + 1 << " " << std::setprecision(12) << (Ax ? Ax[j] : 0);
    if (j + 1 != Ap[n])
     std::cout << "\n";
   }
  }
 }

 void print_csc_l(std::string beg, size_t n,
                  long *Ap, int *Ai, double *Ax) {
  std::setprecision(48);
  std::cout << beg;
  std::cout << n << " " << n << " " << Ap[n] << "\n";
  for (int i = 0; i < n; ++i) {
   for (long int j = Ap[i]; j < Ap[i + 1]; ++j) {
    std::cout << Ai[j] + 1 << " " << i + 1 << " " << Ax[j];
    if (j + 1 != Ap[n])
     std::cout << "\n";
   }
  }
 }

 template<class type>
 void print_vec(std::string header, int beg,
                int end, type *vec) {

  std::cout << header;
  for (int i = beg; i < end; ++i) {
   std::cout << std::setprecision(15) << vec[i] << ", ";
  }
  std::cout << "\n";
 }

/*
 * Print dense matrix, stored col-wise
 *
 */
 template<class type>
 void print_mat(std::string header, int row_no,
                int col_no, type *mat) {
  std::cout << header;
  for (int i = 0; i < col_no; ++i) {
   for (int j = 0; j < row_no; ++j) {
    std::cout << std::setprecision(15) << mat[i * row_no + j] << ", ";
   }
   std::cout << "\n";
  }
  std::cout << "\n";
 }

 int *compute_inv_perm(int n, int *perm, int *pinv) {
  if (n <= 0)
   return NULL;
  for (int i = 0; i < n; ++i) {
   assert(perm[i] >= 0 && perm[i] < n);
   pinv[perm[i]] = i;
  }
  return pinv;
 }


 void combine_perms(int n, int *perm1, int *perm2, int *result) {
  for (int i = 0; i < n; ++i) {
   result[i] = perm1[perm2[i]];
  }
 }

 int zero_diag_removal(size_t n, size_t NNZ, const int *col,
                       const int *row, double *val,
                       double threshold, double reg_diag) {
  int num_zeros = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(val[col[i]]) < threshold) {
    val[col[i]] = reg_diag;
    num_zeros++;
   }
  }
  return num_zeros;
 }

/*
 * return 0 if equal
 */
 int is_vector_equal(int n, double *v1, double *v2,
                     double eps = std::numeric_limits<double>::min()) {
  int num_ineq = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(v1[i] - v2[i]) > eps)
    num_ineq++;
  }
  return num_ineq;
 }

/*
 * if all vector values are zero
 */
 int is_vector_zero(int n, double *v, double threshold) {
  int num_nonzero = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(v[i]) > threshold)
    num_nonzero++;
  }
  return num_nonzero;
 }


/*
 * returns True:1 if equal
 */
 int is_value_equal(double a, double b, double eps) {
  return std::abs(a - b) < eps;
 }


/*
 * removes duplicates of a vector
 */
 void make_unique(std::vector<int> &sn_arrays) {
  std::sort(sn_arrays.begin(), sn_arrays.end());
  std::vector<int>::iterator it;
  it = std::unique(sn_arrays.begin(), sn_arrays.end());
  sn_arrays.resize(std::distance(sn_arrays.begin(), it));
 }

/*
 * Gathering rows from a CSC matrix
 */
 void gather_rows(std::vector<int> row_ids,
                  size_t An, size_t Am,
                  int *Ap,
                  int *Ai, double *Ax,
                  size_t &Bn, size_t &Bm, size_t &Bnz, int *&Bp,
                  int *&Bi, double *&Bx) {
  CSC *AT = new CSC;
  //transpose to access rows easily.
  transpose_unsym(An, Am, Ap, Ai, Ax,
                  AT->nrow, AT->ncol, AT->p, AT->i, AT->x);
  //buidling transpose of B
  int *BTp, *BTi;
  double *BTx;
  Bn = row_ids.size();//reduced row
  Bm = Am;
  BTp = new int[Bn + 1];
  BTp[0] = 0;
  for (int i = 0; i < Bn; ++i) {
   int j = row_ids[i];
   BTp[i + 1] = BTp[i] + (AT->p[j + 1] - AT->p[j]);
  }
  Bnz = BTp[Bn];
  BTi = new int[Bnz];
  BTx = new double[Bnz];
  int nnz_cnt = 0;
  for (int k = 0; k < Bn; ++k) {
   int j = row_ids[k];
   for (int i = AT->p[j]; i < AT->p[j + 1]; ++i) {
    BTi[nnz_cnt] = AT->i[i];
    BTx[nnz_cnt] = AT->x[i];
    nnz_cnt++;
   }
  }
  //Creating the original B
  transpose_unsym(Bm, Bn, BTp, BTi, BTx,
                  Bn, Bm, Bp, Bi, Bx);
  delete[]BTi;
  delete[]BTp;
  delete[]BTx;
  allocateAC(AT, 0, 0, 0, FALSE);
 }


/*
* Build super matrix from A and B,  All are stored in CSC.
*/
 int build_super_matrix_test(CSC *A, CSC *B, double reg_diag,
                             size_t &SM_size, size_t &SM_nz,
                             int *&SMp, int *&SMi, double *&SMx,
                             double *&sm_rhs) {
  int status;
  SM_size = A->ncol + B->nrow;
  SMp = new int[SM_size + 1];
  SMp[0] = 0;
  for (int i = 1; i < A->ncol + 1; ++i) {
   SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
            (B->p[i] - B->p[i - 1]);
  }
  //Adding diagonal for columns with zero values.
  for (int k = A->ncol + 1; k < SM_size + 1; ++k) {
   SMp[k] = SMp[k - 1] + 1;
  }
  SM_nz = SMp[SM_size];
  SMi = new int[SM_nz];
  SMx = new double[SM_nz]();

  int base1 = A->ncol;
  size_t stp = 0;
  for (int j = 0; j < A->ncol; ++j) {
   stp = SMp[j];
   //Adding Hessian
   for (int i = A->p[j]; i < A->p[j + 1]; ++i) {
    SMi[stp] = A->i[i];
    SMx[stp] = A->x[i];
    stp++;
   }
   //Adding inequalities
   if (B->nrow > 0) {
    for (int i = B->p[j]; i < B->p[j + 1]; ++i) {
     SMi[stp] = base1 + B->i[i];
     SMx[stp] = B->x[i];
     stp++;
    }
   }
   assert(stp == SMp[j + 1]);
  }
  //Putting a small value in diagonals
  for (int l = SMp[A->ncol], j = 0; l < SM_nz; ++l, ++j) {
   //SMx[l] = reg_diag;
   SMi[l] = base1 + j;
  }


  sm_rhs = new double[SM_size]();
  //CSC *sKKTt = ptranspose(SM,2,NULL,NULL,0,status);
  //print_csc("skkt\n",SM_size,SMp,SMi,SMx);
  //print_csc("\nskkt Trans\n",sKKTt->ncol,sKKTt->p,sKKTt->i,sKKTt->x);
  //std::cout<<"\n";
  return status;
 }

 int make_lower(size_t An, size_t Anz, int *Ap, int *Ai, double *Ax,
                size_t &Bn, size_t &Bnz, int *&Bp, int *&Bi, double *&Bx) {
  Bn = An;
  Bnz = (Anz - An) / 2 + An;
  Bp = new int[Bn + 1];
  Bp[0] = 0;
  Bi = new int[Bnz];
  Bx = new double[Bnz];
  int nnz_cnt = 0;
  for (int i = 0; i < An; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    if (Ai[j] >= i) {
     Bi[nnz_cnt] = Ai[j];
     Bx[nnz_cnt] = Ax[j];
     nnz_cnt++;
    }
   }
   Bp[i + 1] = nnz_cnt;
  }
  return 1;
 }

 void compressed_set_to_vector(int set_size, int *set_p, int *set_n,
                               std::vector<std::vector<int>> &res_set) {
  res_set.resize(set_size);
  for (int i = 0; i < set_size; ++i) {
   for (int j = set_p[i]; j < set_p[i + 1]; ++j) {
    res_set[i].push_back(set_n[j]);
   }
  }
 }

 void max_min_spmat(size_t An, int *Ap, int *Ai, double *Ax,
                    double &max_val, double &min_val) {
  max_val = -std::numeric_limits<double>::max();
  min_val = std::numeric_limits<double>::max();
  for (int i = 0; i < An; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    if (Ax[j] > max_val) {
     max_val = Ax[j];
    }
    if (Ax[j] < min_val) {
     min_val = Ax[j];
    }
   }
  }
 }

// bool parse_args(int argc, char **argv, std::map<std::string, std::string> &qp_args) {
//  const char *const short_opts = "m:o:l:a:b:c:d:n:p:r:e:t:h";
//  const option long_opts[] = {
//    {"mode",       required_argument, nullptr, 'm'},
//    {"objective",  required_argument, nullptr, 'o'},
//    {"linear",     required_argument, nullptr, 'l'},
//    {"eq",         required_argument, nullptr, 'a'},
//    {"eqb",        required_argument, nullptr, 'b'},
//    {"ineq",       required_argument, nullptr, 'c'},
//    {"ineqb",      required_argument, nullptr, 'd'},
//    {"nasoq",      required_argument, nullptr, 'n'},
//    {"perturb",    required_argument, nullptr, 'p'},
//    {"refinement", required_argument, nullptr, 'r'},
//    {"epsilon",    required_argument, nullptr, 'e'},
//    {"toli",       required_argument, nullptr, 't'},
//    {"help",       no_argument,       nullptr, 'h'},
//    {nullptr,      no_argument,       nullptr, 0}
//  };
//
//  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
//  while (true) {
//   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
//   if (-1 == opt)
//    break;
//
//   switch (opt) {
//    case 'm':
//     qp_args.insert(std::pair<std::string, std::string>("mode", optarg));
//     break;
//    case 'o':
//     qp_args.insert(std::pair<std::string, std::string>("H", optarg));
//     break;
//    case 'l':
//     qp_args.insert(std::pair<std::string, std::string>("q", optarg));
//     break;
//    case 'a':
//     qp_args.insert(std::pair<std::string, std::string>("A", optarg));
//     break;
//    case 'b':
//     qp_args.insert(std::pair<std::string, std::string>("b", optarg));
//     break;
//    case 'c':
//     qp_args.insert(std::pair<std::string, std::string>("C", optarg));
//     break;
//    case 'd':
//     qp_args.insert(std::pair<std::string, std::string>("d", optarg));
//     break;
//    case 'n':
//     qp_args.insert(std::pair<std::string, std::string>("nasoq", optarg));
//     break;
//    case 'p':
//     qp_args.insert(std::pair<std::string, std::string>("perturbation", optarg));
//     break;
//    case 'r':
//     qp_args.insert(std::pair<std::string, std::string>("iterations", optarg));
//     break;
//    case 'e':
//     qp_args.insert(std::pair<std::string, std::string>("epsilon", optarg));
//     break;
//    case 't':
//     qp_args.insert(std::pair<std::string, std::string>("tolerance", optarg));
//     break;
//    case 'h': // -h or --help
//    case '?': // Unrecognized option
//    default:
//     print_help();
//     break;
//   }
//  }
//  return true;
// }
//
// bool parse_nasoq_args(int argc, char **argv, std::map<std::string, std::string> &qp_args) {
//  const char *const short_opts = "i:o:l:d:v:p:r:e:t:h";
//  const option long_opts[] = {
//    {"input",      required_argument, nullptr, 'i'},
//    {"output",      required_argument, nullptr, 'o'},
//    {"log",      required_argument, nullptr, 'l'},
//    {"header",      required_argument, nullptr, 'd'},
//    {"variant",      required_argument, nullptr, 'v'},
//    {"perturb",    required_argument, nullptr, 'p'},
//    {"refinement", required_argument, nullptr, 'r'},
//    {"epsilon",    required_argument, nullptr, 'e'},
//    {"toli",       required_argument, nullptr, 't'},
//    {"help",       no_argument,       nullptr, 'h'},
//    {nullptr,      no_argument,       nullptr, 0}
//  };
//  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
//  while (true) {
//   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
//
//   if (-1 == opt)
//    break;
//   switch (opt) {
//    case 'i':
//     qp_args.insert(std::pair<std::string, std::string>("input", optarg));
//     break;
//    case 'o':
//     qp_args.insert(std::pair<std::string, std::string>("output", optarg));
//     break;
//     case 'l':
//     qp_args.insert(std::pair<std::string, std::string>("log", optarg));
//     break;
//    case 'd':
//     qp_args.insert(std::pair<std::string, std::string>("header", optarg));
//     break;
//     case 'v':
//     qp_args.insert(std::pair<std::string, std::string>("variant", optarg));
//     break;
//    case 'p':
//     qp_args.insert(std::pair<std::string, std::string>("perturbation", optarg));
//     break;
//    case 'r':
//     qp_args.insert(std::pair<std::string, std::string>("iterations", optarg));
//     break;
//    case 'e':
//     qp_args.insert(std::pair<std::string, std::string>("epsilon", optarg));
//     break;
//    case 't':
//     qp_args.insert(std::pair<std::string, std::string>("tolerance", optarg));
//     break;
//    case 'h': // -h or --help
//    case '?': // Unrecognized option
//    default:
//     print_help();
//     break;
//   }
//  }
//  if(qp_args.find("input") == qp_args.end()){
//   print_help();
//   return false;
//  }
//  return true;
// }


}
#endif //CHOLOPENMP_UTIL_H

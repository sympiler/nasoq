//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/Util.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "nasoq/common/cxxopts.hpp"

namespace nasoq {

 bool readMatrix_old(std::string fName, int &n, int &NNZ, int *&col, int *&row, double *&val) {
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

 bool readMatrix(std::string fName, size_t &n, size_t &NNZ, int *&col, int *&row, double *&val) {
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

 bool readMatrix_rect(std::string fName, size_t &n_row, size_t &n_col, size_t &NNZ, int *&col, int *&row, double *&val) {
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

 bool
 expandMatrix(size_t n, size_t NNZ, const int *col, const int *row, const double *val, size_t &newNNZ, int *&newCol,
              int *&newRow, double *&newVal, double insDiag) {
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

 int write_vector(std::string fName, size_t n, double *vec_vals, std::string header, int prec) {
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

 void rhsInit_linearSolver(int n, int *Ap, int *Ai, double *Ax, int *Bp, int *Bi, double *Bx, double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
    b[Ai[cc]] += Ax[cc];
   }
   for (int i = Bp[c] + 1; i < Bp[c + 1]; ++i) {
    b[Bi[i]] += Bx[i];
   }
  }
 }

 bool rhsInit_linearSolver_onlylower(int n, int *Ap, int *Ai, double *Ax, double *b) {
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

 void rhsInit(int n, int *Ap, int *Ai, double *Ax, double *b) {
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j] = 0;
  }
  for (int c = 0; c < n; ++c) {
   for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
    b[Ai[cc]] += Ax[cc];
   }
  }
 }

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

 void
 make_full(int ncol, int nnz, int *Ap, int *Ai, double *Ax, int *ATransp, int *ATransi, double *ATransx, int &nnzFull,
           int *&AFullp, int *&AFulli, double *&AFullx) {
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

 int testTriangular(size_t n, const double *x, double epsilon) {//Testing
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

 int bcsc2csc(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, int *sup2col, double *Ax, int *Cp, int *Ci,
              double *Cx) {
  size_t actualNNZ = 0;
  Cp[0] = 0;
  for (int i = 0; i < nBlocks; ++i) {
   int curCol = sup2col[i];
   int nxtCol = sup2col[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);

   for (int j = curCol; j < nxtCol; ++j) {
    for (int k = Ap[j] + (j - curCol), kk = AiP[curCol] + (j - curCol);
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

 int
 bcsc2csc_aggressive(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, int *sup2col, double *Ax, size_t *Cp,
                     int *Ci, double *Cx) {
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
    for (int k = Ap[j] + (j - curCol), kk = AiP[curCol] + (j - curCol);
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

 int check_row_idx_l(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, int *sup2col) {
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

 int
 bcsc2csc_aggressive_int(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, int *sup2col, double *Ax, int *Cp,
                         int *Ci, double *Cx) {
  size_t actualNNZ = 0;
  Cp[0] = 0;
  for (int i = 0; i < nBlocks; ++i) {
   int curCol = sup2col[i];
   int nxtCol = sup2col[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);

   for (int j = curCol; j < nxtCol; ++j) {
    for (int k = Ap[j] + (j - curCol), kk = AiP[curCol] + (j - curCol);
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

 void print_csc_l(std::string beg, size_t n, long *Ap, int *Ai, double *Ax) {
  std::cout <<  std::setprecision(48) << beg << n << " " << n << " " << Ap[n] << "\n";
  for (int i = 0; i < n; ++i) {
   for (long int j = Ap[i]; j < Ap[i + 1]; ++j) {
    std::cout << Ai[j] + 1 << " " << i + 1 << " " << Ax[j];
    if (j + 1 != Ap[n])
     std::cout << "\n";
   }
  }
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

 int
 zero_diag_removal(size_t n, size_t NNZ, const int *col, const int *row, double *val, double threshold, double reg_diag) {
  int num_zeros = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(val[col[i]]) < threshold) {
    val[col[i]] = reg_diag;
    num_zeros++;
   }
  }
  return num_zeros;
 }

 int is_vector_equal(int n, double *v1, double *v2, double eps) {
  int num_ineq = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(v1[i] - v2[i]) > eps)
    num_ineq++;
  }
  return num_ineq;
 }

 int is_vector_zero(int n, double *v, double threshold) {
  int num_nonzero = 0;
  for (int i = 0; i < n; ++i) {
   if (std::abs(v[i]) > threshold)
    num_nonzero++;
  }
  return num_nonzero;
 }

 int is_value_equal(double a, double b, double eps) {
  return std::abs(a - b) < eps;
 }

 void make_unique(std::vector<int> &sn_arrays) {
  std::sort(sn_arrays.begin(), sn_arrays.end());
  std::vector<int>::iterator it;
  it = std::unique(sn_arrays.begin(), sn_arrays.end());
  sn_arrays.resize(std::distance(sn_arrays.begin(), it));
 }

 void gather_rows(std::vector<int> row_ids, size_t An, size_t Am, int *Ap, int *Ai, double *Ax, size_t &Bn, size_t &Bm,
                  size_t &Bnz, int *&Bp, int *&Bi, double *&Bx) {
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

 int build_super_matrix_test(CSC *A, CSC *B, double reg_diag, size_t &SM_size, size_t &SM_nz, int *&SMp, int *&SMi,
                             double *&SMx, double *&sm_rhs) {
  int status=0;
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
   //SMx[l] = diag_perturb;
   SMi[l] = base1 + j;
  }


  sm_rhs = new double[SM_size]();
  //CSC *sKKTt = ptranspose(SM,2,NULL,NULL,0,status);
  //print_csc("skkt\n",SM_size,SMp,SMi,SMx);
  //print_csc("\nskkt Trans\n",sKKTt->ncol,sKKTt->p,sKKTt->i,sKKTt->x);
  //std::cout<<"\n";
  return status;
 }

 int make_lower(size_t An, size_t Anz, int *Ap, int *Ai, double *Ax, size_t &Bn, size_t &Bnz, int *&Bp, int *&Bi,
                double *&Bx) {
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

 void compressed_set_to_vector(int set_size, int *set_p, int *set_n, std::vector<std::vector<int>> &res_set) {
  res_set.resize(set_size);
  for (int i = 0; i < set_size; ++i) {
   for (int j = set_p[i]; j < set_p[i + 1]; ++j) {
    res_set[i].push_back(set_n[j]);
   }
  }
 }

 void max_min_spmat(size_t An, int *Ap, int *Ai, double *Ax, double &max_val, double &min_val) {
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

 bool parse_args(int argc, const char **argv, std::map<std::string, std::string> &qp_args) {
  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  try
  {
   cxxopts::Options options(argv[0], " - example command line options");
   options
     .positional_help("[optional args]")
     .show_positional_help();

   options
     .allow_unrecognised_options()
     .add_options()
       ("v,variant", "NASOQ Variant", cxxopts::value<std::string>())
       ("i,input", "Input SMP format", cxxopts::value<std::string>()) //nasoq driver
       ("u,output", "Output SMP format", cxxopts::value<std::string>()) // not imp
       ("d,header", "Print CSV header", cxxopts::value<std::string>()) // only for nasoq driver
       ("o,objective", "Quadratic matrix",  cxxopts::value<std::string>())
       ("l,linear", "linear vector",  cxxopts::value<std::string>())
       ("a,eq", "Equality matrix",  cxxopts::value<std::string>() )
       ("b,eqb", "Equality vector", cxxopts::value<std::string>())
       ("c,ineq", "Inequality matrix", cxxopts::value<std::string>())
       ("g,ineqb", "Inequality vector", cxxopts::value<std::string>())
       ("p,perturb", "Pertubation power value",
        cxxopts::value<std::string>())
       ("r,refinement", "Refinement iterations",
        cxxopts::value<std::string>())
       ("e,epsilon", "Accuracy threshold power",
        cxxopts::value<std::string>())
       ("t,toli", "Iterative refinement threshold",
        cxxopts::value<std::string>())
       ("help", "Print help")
#ifdef CXXOPTS_USE_UNICODE
    ("unicode", u8"A help option with non-ascii: à. Here the size of the"
        " string should be correct")
#endif
     ;

   auto result = options.parse(argc, argv);
   if(result.count("v"))
    qp_args.insert(std::pair<std::string, std::string>("variant",
                                                       result["v"].as<std::string>()));

   if(result.count("i"))
    qp_args.insert(std::pair<std::string, std::string>("input",
                                                       result["i"].as<std::string>()));

   if(result.count("o"))
    qp_args.insert(std::pair<std::string, std::string>("H",
                                                       result["o"].as<std::string>()));

   if(result.count("l"))
    qp_args.insert(std::pair<std::string, std::string>("q",
                                                       result["l"].as<std::string>()));

   if(result.count("a"))
    qp_args.insert(std::pair<std::string, std::string>("A",
                                                       result["a"].as<std::string>()));

   if(result.count("b"))
    qp_args.insert(std::pair<std::string, std::string>("b",
                                                       result["b"].as<std::string>()));

   if(result.count("c"))
    qp_args.insert(std::pair<std::string, std::string>("C",
                                                       result["c"].as<std::string>()));

   if(result.count("g"))
    qp_args.insert(std::pair<std::string, std::string>("d",
                                                       result["g"].as<std::string>()));

   if(result.count("p"))
    qp_args.insert(std::pair<std::string, std::string>("perturbation",
                                                       result["p"].as<std::string>()));

   if(result.count("e"))
    qp_args.insert(std::pair<std::string, std::string>("epsilon",
                                                       result["e"].as<std::string>()));

   if(result.count("r"))
    qp_args.insert(std::pair<std::string, std::string>("iterations",
                                                       result["r"].as<std::string>()));

   if(result.count("t"))
    qp_args.insert(std::pair<std::string, std::string>("tolerance",
                                                       result["t"].as<std::string>()));

   if(result.count("d"))
    qp_args.insert(std::pair<std::string, std::string>("header",
                                                       result["d"].as<std::string>()));

   if(result.count("u"))
    qp_args.insert(std::pair<std::string, std::string>("output",
                                                       result["u"].as<std::string>()));
   if (result.count("h"))
   {
    std::cout << "needs to write it! :|" << std::endl;
    exit(0);
   }
  }
  catch (const cxxopts::OptionException& e)
  {
   std::cout << "error parsing options: " << e.what() << std::endl;
   exit(1);
  }
  return true;
 }
}
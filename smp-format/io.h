//
// Created by kazem on 2020-05-15.
//

#ifndef SCO_CONVERTOR_IO_H
#define SCO_CONVERTOR_IO_H

#include <sstream>
#include <vector>
#include <fstream>

#include "def.h"
#include "utils.h"

namespace format{
int precision = 48;
 enum TYPE{
  REAL,INT,COMPLEX,PATTERN
 };

 enum SHAPE{// LOWER and UPPER both are symmetric matrices.
  LOWER,UPPER,GENERAL
 };
 enum FORMAT{
  COORDINATE,ARRAY
 };

 struct triplet{
  int row{}; int col{}; double val{};
 };

 /// Write calls
 /*
  * Printing a compressed matrix, CSC to output
  * output can be stdout or a file (The out will be file.rdbuf())
  */
 template<class type> void print_csc( size_t m, size_t n, int *Ap,
    int *Ai, type *Ax,
    std::streambuf* out = std::cout.rdbuf(),
    const std::string indent="  ",
    const std::string &beg="%%MatrixMarket matrix coordinate real symmetric"){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout<<indent<<beg << "\n";
  size_t nnz = n>0?Ap[n]:0;
  std::cout<<indent<<m<<" "<<n<<" "<<nnz<<"\n";
  for (auto i = 0; i < n; ++i) {
   for (auto j = Ap[i]; j < Ap[i+1]; ++j) {
    assert(j<nnz);
    std::cout<<indent<<Ai[j]+1<<" "<<i+1<<" "<<std::setprecision(precision)<<Ax[j];
    std::cout<<"\n";
   }
  }
  std::cout.rdbuf(sb_cout_backup);
 }

 /*
 * Print dense matrix, stored col-wise
 * The out can be file.rdbuf()
 */
 template<class type> void print_dense( int row_no,
       int col_no, int lda,
       type *mat,
       std::streambuf* out = std::cout.rdbuf(),
       const std::string indent="  ",
       const std::string &header="%%MatrixMarket matrix array real general"){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout<<indent<<header<<"\n";
  std::cout<<indent<<row_no<<" "<<col_no<<"\n";
  for (int i = 0; i < col_no*row_no; i+=lda) {
   for (int j = 0; j < lda; ++j) {
    std::cout<<indent<<std::setprecision(precision)<<mat[i+j];
   }
   std::cout<<"\n";
  }
  std::cout.rdbuf(sb_cout_backup);
 }

 /// prints a constant value with indention
 /// \tparam type
 /// \param val
 /// \param out
 /// \param indent
 template<class type> void print_constant( type val,
                                        std::streambuf* out = std::cout.rdbuf(),
                                        const std::string indent="  "){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout<<indent<<std::setprecision(precision)<<val<<"\n";
  std::cout.rdbuf(sb_cout_backup);
 }


 /// print a string to output with an indent
 /// \param context
 /// \param out
 /// \param indent
 void print_string(std::string context, std::streambuf* out = std::cout.rdbuf(),
      std::string indent = "  "){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::stringstream ss(context);
  std::string line;
  while (std::getline(ss,line)){
   std::cout<<indent<<line<<"\n";
  }
  std::cout.rdbuf(sb_cout_backup);
 }

 /// Read calls

/// Reads the header information of a mtx format.
/// \param inFile input stream
/// \param n_row
/// \param n_col
/// \param n_nnz
/// \param type
/// \param shape
/// \param mtx_format
/// \return
 int read_header(std::ifstream &inFile, int &n_row, int &n_col,
   size_t &n_nnz, int &type, int &shape, int &mtx_format){
  std::string line,banner, mtx, crd, arith, sym;
  std::getline(inFile,line);
  trim(line);
  for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
  std::istringstream iss(line);
  if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
   std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
   return false;
  }
  if(!(banner =="%%matrixmarket")) {
   std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")\n";
   return false;
  }
  if(!(mtx =="matrix")) {
   std::cout<<"Not a matrix; this driver cannot handle that.\"\n";
   return false;
  }
  if(crd == "coordinate") {
   mtx_format = COORDINATE;
  } else if(crd == "array") {
   mtx_format = ARRAY;
  } else{
   std::cout<<"Not in coordinate format; this driver cannot handle that.\"\n";
   return false;
  }
  if(arith == "real")
   type = REAL;
  else if(arith == "integer")
   type = INT;
  else if (arith == "complex")
   type = COMPLEX;
  else if(arith == "pattern")
   type = PATTERN;
  else{
   std::cout<<"Unknown arithmetic\n";
   return false;
  }
  if(sym == "symmetric")
   shape = LOWER;
  else if(sym == "general")
   shape = GENERAL;
  else{
   std::cout<<"Unknown shape\n";
   return false;
  }
  while (!line.compare(0,1,"%"))
  {
   std::getline(inFile, line);
   trim(line);
  }
  std::istringstream issDim(line);
  if(mtx_format != ARRAY){
   if (!(issDim >> n_row >> n_col >> n_nnz)){
    std::cout<<"The matrix dimension is missing\n";
    return false;
   }
  } else{
   if (!(issDim >> n_row >> n_col)){
    std::cout<<"The matrix dimension is missing\n";
    return false;
   }
   n_nnz = n_row*n_col;
  }
  return true;
 }


/// Builds triplet from coordinate mtx file
/// \param inFile
/// \param nnz
/// \param triplet_vec
/// \param zero_indexing
 void read_triplets_real(std::ifstream &inFile, int nnz,
   std::vector<triplet>& triplet_vec, bool zero_indexing=false){
  for (int i = 0; i < nnz; ++i) {
   triplet tmp;
   inFile >> tmp.row;
   inFile >> tmp.col;
   inFile >> tmp.val;
   if(!zero_indexing){
    tmp.col--; tmp.row--;
   }
   triplet_vec.push_back(tmp);
  }
 }

 void compress_triplets_to_csc(std::vector<triplet>& triplet_vec, CSC *A,
   bool add_diags= true){
  assert(A->nnz == triplet_vec.size());
  std::sort(triplet_vec.begin(), triplet_vec.end(),
    [](const triplet& a, const triplet& b){return (a.col<b.col) || (a.col==b.col && a.row<b.row);});
  auto *count = new int[A->n]();
  for (auto i = 0; i < A->nnz; ++i) {
   count[triplet_vec[i].col]++;
  }
  A->p[0] = 0;
  for (auto j = 0; j < A->n; ++j) {
   if(count[j] == 0 && add_diags){ // insert zero diag for empty cols
    triplet tmp; tmp.col = tmp.row = j; tmp.val=0;
    triplet_vec.insert(triplet_vec.begin()+A->p[j], tmp);
    A->p[j+1] = A->p[j] + 1;
   }else{
    A->p[j+1] = A->p[j] + count[j];
   }
  }
  delete []count;
  for (auto k = 0; k < A->nnz; ++k) {
   A->i[k] = triplet_vec[k].row;
   A->x[k] = triplet_vec[k].val;
  }
 }

 bool read_mtx_csc_real(std::ifstream &in_file, CSC *&A, bool insert_diag=false){
  int n, m;
  int shape, arith, mtx_format;
  size_t nnz;
  std::vector<triplet> triplet_vec;

  bool ret = read_header(in_file, m, n, nnz, arith, shape, mtx_format);
  if(arith != REAL || !ret || mtx_format != COORDINATE)
   return false;
  A = new CSC(m,n,nnz,false, shape == LOWER);
  read_triplets_real(in_file, nnz, triplet_vec);
  compress_triplets_to_csc(triplet_vec, A, insert_diag);
  A->nnz = A->p[n]; // if insert diag is true, it will be different.
  //print_csc(A->n, A->p, A->i, A->x);
  return true;
 }

 /*
 * Reads an array stored in matrix market format.
 */
 bool read_mtx_array_real(std::ifstream &in_file, Dense *&A) {
  int n, m;
  int shape, arith, mtx_format;
  size_t nnz;
  bool ret = read_header(in_file, m, n, nnz, arith, shape, mtx_format);
  if(arith != REAL || !ret || mtx_format != ARRAY)
   return false;
  A = new Dense(m, n, 1);//
  for (int i = 0; i < m * n; i++) {//writing from file row by row
   in_file >> A->a[i];
  }
  std::ofstream file;
  //print_dense(A->row, A->col, A->lda, A->a);
  return true;
 }

/// Reads a real constant from input
/// \param in_file
/// \param val
 void read_real_constant(std::ifstream &in_file, double &val){
  in_file >> val;
 }

 /// Reading sting from a file with indention
 /// \param in_file
 /// \param context
 /// \param indent
 void read_string(std::ifstream &in_file, std::string &context, std::string indent = "  "){
  std::string line;
  while (!in_file.eof()){
   std::streampos oldpos = in_file.tellg();  // stores the position
   std::getline(in_file,line);
   if(line.compare(0, indent.size(), indent)){
    if(line[0] == '-') //TODO not sure if this is always true
     continue;
    in_file.seekg (oldpos);   // get back to the position
    break;
   } else{
    ltrim(line);
    line+="\n";
    context.append(line);
   }
  }
 }


}

#endif //SCO_CONVERTOR_IO_H

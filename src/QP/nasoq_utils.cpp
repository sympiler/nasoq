//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/QP/nasoq_utils.h"

#include <limits>
#include <iostream>
#include <cmath>
namespace nasoq {

 int max_min_sparse_matrix(const CSC *M, double &max_v, double &min_v, double &avg, double &var, double m_val) {
  size_t NNZ = M->nzmax;
  min_v = std::numeric_limits<double>::max();
  max_v = -std::numeric_limits<double>::max();
  avg = 0;
  for (int i = 0; i < NNZ; ++i) {
   double tmp = std::abs(M->x[i]);
   if (tmp > max_v)
    max_v = tmp;
   else if (tmp < min_v && tmp > m_val)
    min_v = tmp;
   avg += (tmp / NNZ);
  }
  /*if(max_v == -std::numeric_limits<double >::max())
   max_v=0;
  if(min_v == std::numeric_limits<double >::max())
   min_v=0;*/
  var = 0;
  for (int i = 0; i < NNZ; ++i) {
   double tmp = std::abs(M->x[i]) - avg;
   var += ((tmp * tmp) / NNZ);
  }
  return 1;
 }

 int max_min_vector(int n, double *v, double &max_v, double &min_v, double &avg, double &var, double m_val) {
  var = 0;
  avg = 0;
  max_v = -std::numeric_limits<double>::max();
  min_v = std::numeric_limits<double>::max();
  for (int i = 0; i < n; ++i) {
   double tmp = std::abs(v[i]);
   if (tmp > max_v)
    max_v = tmp;
   else if (tmp < min_v && tmp > m_val)
    min_v = tmp;
   avg += (tmp / n);
  }
/* if(max_v == -std::numeric_limits<double >::max())
  max_v=0;
 if(min_v == std::numeric_limits<double >::max())
  min_v=0;*/
  for (int i = 0; i < n; ++i) {
   double tmp = std::abs(v[i]) - avg;
   var += ((tmp * tmp) / n);
  }
  return 0;

 }

 int scale_sparse_matrix(CSC *M, double scalar) {
  size_t NNZ = M->nzmax;
  for (int i = 0; i < NNZ; ++i) {
   M->x[i] = scalar * M->x[i];
  }
  return 1;
 }

 void scale_vector(int n, double *vec, double scalar) {
  for (int i = 0; i < n; ++i) {
   vec[i] *= scalar;
  }
 }
}

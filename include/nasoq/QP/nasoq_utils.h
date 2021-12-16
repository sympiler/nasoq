//
// Created by kazem on 9/22/19.
//

#ifndef PARS_NASOQ_UTILS_H
#define PARS_NASOQ_UTILS_H

#include "nasoq/common/def.h"

namespace nasoq {
 int max_min_sparse_matrix(const CSC *M, double &max_v,
                           double &min_v, double &avg, double &var, double m_val = 1e-15);

 int max_min_vector(int n, double *v, double &max_v,
                    double &min_v, double &avg, double &var, double m_val = 1e-15);

 int scale_sparse_matrix(CSC *M, double scalar);

 void scale_vector(int n, double *vec, double scalar);
}
#endif //PARS_NASOQ_UTILS_H

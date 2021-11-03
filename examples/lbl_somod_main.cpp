//
// Created by kazem on 2/21/21.
//

#include <iostream>

#include <nasoq/nasoq.h>
#include <nasoq/QP/linear_solver_wrapper.h>
#include <cmath>


int main(int argc, char** argv) {
 /// Declaring input matrices
 size_t sizeH;
 size_t nnzH;
 sizeH = 2;
 nnzH = 2;
 auto *q = new double[sizeH];
 auto *Hp = new int[sizeH+1];
 auto *Hi = new int[nnzH];
 auto *Hx = new double[nnzH];

 q[0] = -4; q[1] = -4;

 Hp[0]=0;Hp[1]=1;Hp[2]=2;
 Hi[0]=0;Hi[1]=1;
 Hx[0]=2;Hx[1]=2;

 auto *H = new nasoq::CSC; H->nzmax = nnzH; H->ncol= H->nrow =sizeH;
 H->p = Hp; H->i = Hi; H->x = Hx; H->stype=-1; H->packed=1;

 size_t rowsC;
 size_t nnzC;
 rowsC = 2;
 nnzC = 2;
 auto* d = new double[sizeH];
 auto *Cp = new int[sizeH+1];
 auto *Ci = new int[nnzC];
 auto *Cx = new double[nnzC];

 d[0] = 0; d[1] = 0;

 Cp[0]=0;Cp[1]=1;Cp[2]=2;
 Ci[0]=0;Ci[1]=1;
 Cx[0]=1;Cx[1]=1;

 auto *C = new nasoq::CSC; C->nzmax = nnzC; C->ncol = sizeH; C->nrow = rowsC;
 C->p = Cp; C->i = Ci; C->x = Cx; C->stype=0; C->packed=1;

 /// Solving the linear system
 auto *lbl = new nasoq::SolverSettings(H,q,C,d);
 lbl->ldl_variant = 2; // set it 4 if your C++ compiler supports OpenMP
 lbl->req_ref_iter = 2;
 lbl->solver_mode = 1;
 lbl->reg_diag = pow(10,-9);
 lbl->symbolic_analysis();
 lbl->numerical_factorization();
 double *x = lbl->solve_only();


 /// Printing results with no modification
 // expected x={-2,-2};
 std::cout<<"Solution: ";
 for (int i = 0; i < sizeH; ++i) {
  std::cout<<x[i]<<",";
 }

 /// Updating with constraint 0
 std::vector<int> constraint_inds;
 constraint_inds.push_back(0);
 d[0] = 2; d[1] = 1;
 lbl->update_somod(constraint_inds, 1, d);
 lbl->update_factorization();
 x = lbl->solve_only();
 /// Printing results
 // expected x={0,-2, -4};
 std::cout<<"Solution: ";
 for (int i = 0; i < sizeH + 2; ++i) {
  std::cout<<x[i]<<",";
 }

 /// Updating with constraint 1
 constraint_inds.clear(); constraint_inds.push_back(1);
 lbl->update_somod(constraint_inds, 1, d);
 lbl->update_factorization();
 x = lbl->solve_only();
 /// Printing results
 // expected x={0,0,-4,-4};
 std::cout<<"Solution: ";
 for (int i = 0; i < sizeH + 2; ++i) {
  std::cout<<x[i]<<",";
 }

 /// Downdating constraint 0
 constraint_inds.clear(); constraint_inds.push_back(0);
 lbl->downdate_somod(constraint_inds, 1);
 lbl->update_factorization();
 x = lbl->solve_only();
 /// Printing results
 // expected x={-2, 0, 0, -4};
 std::cout<<"Solution: ";
 for (int i = 0; i < sizeH + 2; ++i) {
  std::cout<<x[i]<<",";
 }

 delete lbl;
 delete []Hp;
 delete []Hi;
 delete []Hx;
 delete []q;
 //delete []x;
 delete H;
 delete []Cp;
 delete []Ci;
 delete []Cx;
 delete []d;
 delete C;
 return 0;
}

//
// Created by kazem on 11/18/20.
//


#include <nasoq/nasoq_step.h>
#include <iostream>
#include <cmath>

/*
 * Minimizing 1/2 x^THx + q^Tx + C; Cx <= d
 * H and C are sparse CSC matrices
 * q and d are dense arrays
 */

int main(int argc, char *argv[]){
 /// Declaring inputs
 size_t sizeH;
 size_t nnzH;
 size_t CRows; size_t CCols; size_t nnzA;
 sizeH = 2; nnzH = 2;
 CRows = 4; CCols = 2; nnzA = 8;
 auto *q = new double[sizeH];
 auto *Hp = new int[sizeH+1];
 auto *Hi = new int[nnzH];
 auto *Hx = new double[nnzH];
 auto *Cp = new int[CCols+1];
 auto *Ci = new int[nnzA];
 auto *Cx = new double[nnzA];
 auto *d = new double[CRows];

 q[0] = -4; q[1] = -4;

 Hp[0]=0;Hp[1]=1;Hp[2]=2;
 Hi[0]=0;Hi[1]=1;
 Hx[0]=2;Hx[1]=2;

 Cp[0]=0;Cp[1]=4;Cp[2]=8;
 Ci[0]=0;Ci[1]=1;Ci[2]=2;Ci[3]=3;
 Ci[4]=0;Ci[5]=1;Ci[6]=2;Ci[7]=3;
 Cx[0]=2;Cx[1]=1;Cx[2]=-1;Cx[3]=-2;
 Cx[4]=1;Cx[5]=-1;Cx[6]=-1;Cx[7]=1;

 d[0]=2;d[1]=1;d[2]=1;d[3]=2;

 /// Solving the QP pronlem
 nasoq::NasoqStep *qm;
 qm = new nasoq::NasoqStep(sizeH,Hp,Hi,Hx,
                       q,CRows,CCols,Cp, Ci, Cx, d);
 qm->diag_perturb=pow(10,-9);
 qm->eps_abs=pow(10,-3);
 qm->max_iter = 0;
 qm->variant = nasoq::PREDET;
 if(qm->solve_init()){
  std::cout<<"Initialization: \n";
  std::cout<<qm->primal_vars[0] <<","<<qm->primal_vars[1]<<"\n";
  std::cout<<qm->dual_vars[0] <<","<<qm->dual_vars[1]<<"\n";
 }

  while( true ){
   auto ret_val = qm->solve_step();
   if(ret_val != nasoq::NotFinished)
    break;
   std::cout<<"Iteration :"<< qm->num_iter<<"\n";
   std::cout<<qm->primal_vars[0] <<","<<qm->primal_vars[1]<<"\n";
   std::cout<<qm->dual_vars[0] <<","<<qm->dual_vars[1]<<"\n";
  }


 /// Printing results
 if(qm->ret_val == nasoq::Optimal)
  std::cout<<"The problem is converged\n";

 // expected x={0.4,1.2};
 auto *x = qm->primal_vars;
 std::cout<<"Primal variables: ";
 for (int i = 0; i < sizeH; ++i) {
  std::cout<<x[i]<<",";
 }

 // expected z = {1.6,0,0,0}
 std::cout<<"\nDual variables: ";
 auto *z = qm->dual_vars;
 for (int i = 0; i < CRows; ++i) {
  std::cout<<z[i]<<",";
 }

 delete qm;
 delete []Hp; delete []Hi; delete []Hx; delete []q;
 delete []Cp;  delete []Ci; delete []Cx; delete []d;
 return 0;
}


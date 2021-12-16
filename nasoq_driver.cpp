//
// Created by kazem on 3/6/19.
//

#include <string>
#include <iostream>
#include <cassert>

#ifdef MKL_BLAS
#include "mkl.h"
#endif

#include "nasoq/nasoq.h"

#include "utils.h"
#include <mp_format_converter.h>

using namespace nasoq;

int nasoq_demo(int argc, const char **argv);

int main(int argc, const char *argv[]) {

 nasoq_demo(argc, argv);

}


void get_num_iter(double acc_thresh, int& inner_iter_ref,
                  int& outer_iter_ref){
 if(acc_thresh >= 1e-3){
  inner_iter_ref = outer_iter_ref = 1;
 } else if (acc_thresh <= 1e-4 && acc_thresh > 1e-10){
  inner_iter_ref = outer_iter_ref = 2;
 } else if (acc_thresh <= 1e-11 && acc_thresh > 1e-13){
  inner_iter_ref = outer_iter_ref = 3;
 } else{
  inner_iter_ref = outer_iter_ref = 9;
 }
}


/*
 *
 */
int nasoq_demo(int argc, const char **argv) {
 std::map<std::string, std::string> qp_args;
 if (!format::parse_args(argc, argv, qp_args))
  return -1;

 /// New settings if provided
 std::string input_qp_path = qp_args["input"];
 double reg_diag = pow(10,-9);
 double zero_threshold = reg_diag;
 double eps = 1e-6;
 int inner_iter = 2;
 int outer_iter = 2;
 double stop_tol = 1e-15;
 int solver_mode = Fixed;
 bool print_header = false;
 if (qp_args.find("variant") != qp_args.end()) {
  std::string nasoq_mode = qp_args["variant"];
  if (nasoq_mode == "tuned") {
   solver_mode = Tuned;
  } else if (nasoq_mode == "auto"){
   solver_mode = AUTO;
  } else if (nasoq_mode == "predet") {
   solver_mode = PREDET;
  } else {
   solver_mode = Fixed;
  }
 }
 if(qp_args.find("perturbation") != qp_args.end())
  reg_diag = pow(10, std::stoi(qp_args["perturbation"]) );
 if(qp_args.find("iterations") != qp_args.end())
  inner_iter = std::stoi(qp_args["iterations"]);
 outer_iter = inner_iter;
 if(qp_args.find("epsilon") != qp_args.end())
  eps = pow(10, std::stoi(qp_args["epsilon"]) );
 if(qp_args.find("tolerance") != qp_args.end())
  stop_tol = pow(10, std::stoi(qp_args["tolerance"]) );
 if(qp_args.find("header") != qp_args.end())
  print_header = true;


 auto *QPFC = new format::QPFormatConverter();
 if(!QPFC->load_smp(input_qp_path))
  return -1;
 QPFC->smp_to_ie();
 int num_ineq = QPFC->num_ineq_constraints();
#if defined(MKL_BLAS)
 int num_thread = mkl_get_max_threads ();
 MKL_Set_Num_Threads(num_thread);
#elif defined(_OPENMP)
 int num_thread = omp_get_max_threads();
 omp_set_num_threads(num_thread);
#else
  int num_thread = 1;
#endif
//QPFC->ief_->print();
 Nasoq *qm;
 qm = new nasoq::Nasoq(QPFC->ief_->H->n,QPFC->ief_->H->p,QPFC->ief_->H->i,QPFC->ief_->H->x,
                       QPFC->ief_->q->a,QPFC->ief_->A->m,QPFC->ief_->A->n,QPFC->ief_->A->p,QPFC->ief_->A->i,
                       QPFC->ief_->A->x,QPFC->ief_->b ? QPFC->ief_->b->a : NULLPNTR,
                       QPFC->ief_->C->m,QPFC->ief_->C->n,QPFC->ief_->C->p,
                       QPFC->ief_->C->i,QPFC->ief_->C->x,
                       QPFC->ief_->d ? QPFC->ief_->d->a : NULLPNTR);
 qm->diag_perturb=reg_diag;
 qm->eps_abs=eps;
 qm->max_iter = inner_iter;
 qm->stop_tol = stop_tol;
 qm->max_iter_nas=200000;
 int converged;
 //solver_mode = 0;
 if(solver_mode==0 || solver_mode==1 || solver_mode==3){
  if(solver_mode == 0)
   qm->variant = Fixed;
  else if(solver_mode == 1)
   qm->variant = AUTO;
  else if(solver_mode == 3)
   qm->variant = PREDET;
  converged = qm->solve();
 }else{
  std::chrono::time_point<std::chrono::system_clock> tot_st, tot_end;
  std::chrono::duration<double> elapsed_seconds;
  tot_st = std::chrono::system_clock::now();
  double etime = qm->qi->tot;
  qm->variant = PREDET;
  qm->max_iter_nas=200000;
  if(num_ineq >= 10000){
   qm->inner_iter_ref = 9;
   qm->outer_iter_ref = 9;
   qm->diag_perturb=qm->zero_thresh=pow(10,-11);
   qm->stop_tol = pow(10,-20);
  }else{
   get_num_iter(eps,qm->inner_iter_ref,qm->outer_iter_ref);
   qm->diag_perturb=qm->zero_thresh=pow(10,-9);
   qm->stop_tol = pow(10,-15);
  }

  converged = qm->solve();
  if(converged != nasoq_status::Optimal){
   delete qm;

   std::vector<nasoq_config> configs;
   int itr_in, itr_out;
   get_num_iter(eps,itr_in,itr_out);
   /*
    * NASOQ	/usr/UFDB//MarosMezaros//QPCBOEI1_q
2	2	1	0.000000001
3	3	1	1E-10
3	3	1	1E-11
4	4	1E-18	1E-10
3	3	1E-17	1E-10

/usr/UFDB//MarosMezaros//QPCBOEI2_q
4	4	1e-18	1e-12
4	4	1E-18	1E-10
3	3	1E-17	1E-10
    4	4	1.00E-18	1.00E-10


/DUALC1
2	2	1E-17	1E-11
2	2	1E-18	1E-11




    YAO
    9	9	1E-18	0.000000001
9	9	1E-17	0.000000001
10	10	1E-18	0.000000001
10	10	1E-17	0.000000001
4	4	1E-18	1E-10
4	4	1E-17	1E-10
9	9	1E-18	1E-10
10	10	1E-18	1E-10
3	3	1E-18	1E-11
3	3	1E-17	1E-11
3	3	1E-16	1E-11

    lamb_lamb1
    4	4	1.00E-17	1.00E-11
    4	4	1.00E-18	1.00E-11
    3	3	1.00E-17	1.00E-11

    */
   if(num_ineq >= 10000)
    configs.emplace_back(nasoq_config(itr_in,itr_out,pow(10,-9),pow(10,-15)));
   configs.emplace_back(nasoq_config(2,2,pow(10,-13),pow(10,-15)));
   configs.emplace_back(nasoq_config(itr_in,itr_out,pow(10,-7),pow(10,-15)));
   configs.emplace_back(nasoq_config(3,3,pow(10,-10),pow(10,-15)));
   configs.emplace_back(nasoq_config(3,3,pow(10,-10),pow(10,-17)));
   configs.emplace_back(nasoq_config(3,3,pow(10,-11),pow(10,-17)));
   configs.emplace_back(nasoq_config(2,2,pow(10,-11),pow(10,-17)));
   for (int i = 0; i < configs.size(); ++i) {
    nasoq_config nc = configs[i];
    qm = new nasoq::Nasoq(QPFC->ief_->H->n,QPFC->ief_->H->p,QPFC->ief_->H->i,
      QPFC->ief_->H->x,QPFC->ief_->q->a,QPFC->ief_->A->m,QPFC->ief_->A->n,
      QPFC->ief_->A->p,QPFC->ief_->A->i,QPFC->ief_->A->x,
      QPFC->ief_->b ? QPFC->ief_->b->a : NULLPNTR,QPFC->ief_->C->m,
      QPFC->ief_->C->n,QPFC->ief_->C->p,QPFC->ief_->C->i,QPFC->ief_->C->x,
      QPFC->ief_->d ? QPFC->ief_->d->a : NULLPNTR);
    qm->variant = PREDET;
    qm->diag_perturb=nc.pert_diag;
    qm->batch_size = 1;
    qm->eps_abs=eps;
    qm->eps_rel=eps;
    qm->scaling=0;
    qm->zero_thresh = nc.pert_diag;
    qm->inner_iter_ref = nc.inner_iter;
    qm->outer_iter_ref = nc.outer_iter;
    qm->stop_tol = nc.stop_tol;
    qm->max_iter_nas=200000;
    converged = qm->solve();
    if(i == configs.size()-1){
     tot_end = std::chrono::system_clock::now();
     qm->variant = Tuned;
     break;
    }
    if(converged != nasoq_status::Optimal ){
     delete qm;
    } else {
     tot_end = std::chrono::system_clock::now();
     qm->variant = Tuned;
     break;
    }
   }
   elapsed_seconds = tot_end-tot_st;
   etime += elapsed_seconds.count();
   qm->qi->tot=etime;

  } else{
   qm->variant = Tuned;
  }

 }
 qm->detect_solver_name();
 //qm->print(1);
 std::string heade_vec = "%%MatrixMarket matrix array real general \n";
 heade_vec += std::to_string(qm->H->ncol);
 heade_vec += " 1\n";
// qm->export_to_file(QPFC->problem_name,heade_vec);
 //std::cout<<"dual res norm: "<<qm->lagrangian_residual_norm()<<"\n";
 //std::cout<<"constraint sat norm: "<<qm->constraint_sat_norm()<<"\n";
 //std::cout<<"primal FB norm: "<<qm->primal_FB_norm()<<"\n";
 //std::cout<<"dual FB norm: "<<qm->dual_FB_norm()<<"\n";
 //std::cout<<"obj: "<<qm->objective<<"\n";
 //std::cout<<QPFC->A_eq->nrow<<"\n";

if(print_header){
 std::cout<<"Tool Name,Problem Name,Hessian dim,Hessian NNZ,# of Eq Constraints,"
            "Eq Constraint NNZ,# of Ineq Const,Ineq Constraint NNZ,# of Threads,"
            "eps_abs,Outer GMRES Iter,Inner GMRES Iter,GMRES Tol,Diagonal Pert,"
            "Status,# of Iterations,Time (s),Active-set Size,Constraint Satisfaction Inf,"
            "Residual Lagrangian inf,Primal Obj,Dual Obj,Obj Value,Non-negativity Inf,Complementarity Inf,"
            "Problem Type,Problem Class,\n";
}
 std::cout<<qm->sol_name<<",";
 QPFC->print_log();
 std::cout<<num_thread<<",";
 qm->print_log();
 std::cout<<QPFC->smp_->desc_struct_.application_<<",";
 std::cout<<QPFC->smp_->desc_struct_.category_<<",";
 //qm->print_qp_range(1);
/* std::cout<<qm->eps_abs<<";"
          <<qm->num_iter<<";"<<qm->qi->tot<<";"<<num_thread<<";"
          <<qm->active_set.size()<<";;";*/
 //std::cout<<qm->qi->factt<<";";
 //qm->ss->psi->print_profiling();


 delete qm;
 delete QPFC;
 return 0;
}

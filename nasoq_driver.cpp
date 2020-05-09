//
// Created by kazem on 3/6/19.
//

#include <qp_utils.h>
#include <nasoq.h>
#include <Norm.h>
#include "qp_format_converter.h"


int QP_demo01(int argc, char **argv);

int main(int argc, char *argv[]) {

 QP_demo01(argc, argv);

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
int QP_demo01(int argc, char **argv){
 if (argc < 1) {
  std::cout << "missing args!\n";
  return -1;
 }

 double reg_diag=1e-9;
 double zero_threshold = reg_diag;
 double eps = 1e-6;
 int inner_iter = 2;
 int outer_iter = 2;
 double stop_tol = 1e-15;
 int solver_mode = 1;
 std::string hessian_file = "";
 std::string linear_file = "";
 std::string ineq_file_l = "";
 std::string ineq_file = "";
 std::string ineq_file_u = "";
 std::string eq_file = "";
 std::string eq_file_u = "";

 QPFormatConverter *QPFC = new QPFormatConverter();

 int qp_type = atoi(argv[1]);// 1=IE format, 2 general ineq
 if(qp_type == 1){
  if (argc < 8) {
   std::cout << "missing args!\n";
   return -1;
  }
  hessian_file = argv[2];
  linear_file = argv[3];
  eq_file = argv[4];
  eq_file_u = argv[5];
  ineq_file = argv[6];
  ineq_file_u = argv[7];
  if(argc>8){
   int reg_diag_in = atoi(argv[8]);
   int outer_iter_in = atoi(argv[9])-1;
   int inner_iter_in = atoi(argv[10])-1;
   int eps_in = atoi(argv[11]);
   int tol_in = atoi(argv[12]);
   int sol_mode = 0;
   if(argc > 13)
    sol_mode = atoi(argv[13]);
   reg_diag = pow(10,-reg_diag_in);
   zero_threshold = reg_diag;
   eps = pow(10,-eps_in);
   stop_tol = pow(10,-tol_in);
   inner_iter = inner_iter_in;
   outer_iter = outer_iter_in;
   solver_mode = sol_mode;
  }
  QPFC->read_IE_format(hessian_file,linear_file,
                       eq_file,eq_file_u,ineq_file,
                       ineq_file_u);

 } else if(qp_type==2){
  if (argc < 7) {
   std::cout << "missing args!\n";
   return -1;
  }
  hessian_file = argv[2];
  linear_file = argv[3];
  ineq_file_l = argv[4];
  ineq_file = argv[5];
  ineq_file_u = argv[6];

  if(argc>7){
   int reg_diag_in = atoi(argv[7]);
   int outer_iter_in = atoi(argv[8])-1;
   int inner_iter_in = atoi(argv[9])-1;
   int eps_in = atoi(argv[10]);
   int tol_in = atoi(argv[11]);
   int sol_mode = 0;
   if(argc > 12)
    sol_mode = atoi(argv[12]);
   reg_diag = pow(10,-reg_diag_in);
   zero_threshold = reg_diag;
   eps = pow(10,-eps_in);
   stop_tol = pow(10,-tol_in);
   inner_iter = inner_iter_in;
   outer_iter = outer_iter_in;
   solver_mode = sol_mode;
  }
  QPFC->read_bounded_format(hessian_file,linear_file,
                            ineq_file_l,ineq_file,ineq_file_u);
  //QPFC->print_bounded_format();

  QPFC->B2IE();
  //QPFC->print_IE_format();
  //QPFC->IE_export_to_dense(linear_file);

 }else{
  return -1;
 }

 int num_thread = mkl_get_max_threads ();
 omp_set_num_threads(num_thread);
 MKL_Set_Num_Threads(num_thread);

 Nasoq *qm;
 qm = new Nasoq(QPFC->H->ncol,QPFC->H->p,QPFC->H->i,QPFC->H->x,QPFC->q,
                    QPFC->A_eq->nrow,QPFC->A_eq->ncol,QPFC->A_eq->p,QPFC->A_eq->i,
                    QPFC->A_eq->x,QPFC->a_eq,
                    QPFC->A_ineq->nrow,QPFC->A_ineq->ncol,QPFC->A_ineq->p,
                    QPFC->A_ineq->i,
                    QPFC->A_ineq->x,QPFC->a_ineq);
 qm->reg_diag=reg_diag;
 qm->batch_size = 1;
 qm->eps_abs=eps;
 qm->eps_rel=eps;
 qm->scaling=0;
 qm->zero_thresh = zero_threshold;
 qm->inner_iter_ref = inner_iter;
 qm->outer_iter_ref = outer_iter;
 qm->tol_ref = stop_tol;
 qm->max_iter=200000;
 int converged;
 //solver_mode = 0;
 if(solver_mode==0 || solver_mode==1 || solver_mode==3){
  if(solver_mode == 0)
   qm->nas_mode = Fixed;
  else if(solver_mode == 1)
   qm->nas_mode = AUTO;
  else if(solver_mode == 3)
   qm->nas_mode = PREDET;
  converged = qm->solve();
 }else{
  std::chrono::time_point<std::chrono::system_clock> tot_st, tot_end;
  std::chrono::duration<double> elapsed_seconds;
  tot_st = std::chrono::system_clock::now();
  double etime = qm->qi->tot;
  qm->nas_mode = PREDET;
  qm->max_iter=200000;
  if(QPFC->A_ineq->nrow >= 10000){
   qm->inner_iter_ref = 9;
   qm->outer_iter_ref = 9;
   qm->reg_diag=qm->zero_thresh=pow(10,-11);
   qm->tol_ref = pow(10,-20);
  }else{
   get_num_iter(eps,qm->inner_iter_ref,qm->outer_iter_ref);
   qm->reg_diag=qm->zero_thresh=pow(10,-9);
   qm->tol_ref = pow(10,-15);
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
   if(QPFC->A_ineq->nrow >= 10000)
    configs.emplace_back(nasoq_config(itr_in,itr_out,pow(10,-9),pow(10,-15)));
   configs.emplace_back(nasoq_config(2,2,pow(10,-13),pow(10,-15)));
   configs.emplace_back(nasoq_config(itr_in,itr_out,pow(10,-7),pow(10,-15)));
   configs.emplace_back(nasoq_config(3,3,pow(10,-10),pow(10,-15)));
   configs.emplace_back(nasoq_config(3,3,pow(10,-10),pow(10,-17)));
   configs.emplace_back(nasoq_config(3,3,pow(10,-11),pow(10,-17)));
   configs.emplace_back(nasoq_config(2,2,pow(10,-11),pow(10,-17)));
   for (int i = 0; i < configs.size(); ++i) {
    nasoq_config nc = configs[i];
    qm = new Nasoq(QPFC->H->ncol,QPFC->H->p,QPFC->H->i,QPFC->H->x,QPFC->q,
                       QPFC->A_eq->nrow,QPFC->A_eq->ncol,QPFC->A_eq->p,QPFC->A_eq->i,
                       QPFC->A_eq->x,QPFC->a_eq,
                       QPFC->A_ineq->nrow,QPFC->A_ineq->ncol,QPFC->A_ineq->p,
                       QPFC->A_ineq->i,
                       QPFC->A_ineq->x,QPFC->a_ineq);
    qm->nas_mode = PREDET;
    qm->reg_diag=nc.pert_diag;
    qm->batch_size = 1;
    qm->eps_abs=eps;
    qm->eps_rel=eps;
    qm->scaling=0;
    qm->zero_thresh = nc.pert_diag;
    qm->inner_iter_ref = nc.inner_iter;
    qm->outer_iter_ref = nc.outer_iter;
    qm->tol_ref = nc.stop_tol;
    qm->max_iter=200000;
    converged = qm->solve();
    if(i == configs.size()-1){
     tot_end = std::chrono::system_clock::now();
     qm->nas_mode = Tuned;
     break;
    }
    if(converged != nasoq_status::Optimal ){
     delete qm;
    } else {
     tot_end = std::chrono::system_clock::now();
     qm->nas_mode = Tuned;
     break;
    }
   }
   elapsed_seconds = tot_end-tot_st;
   etime += elapsed_seconds.count();
   qm->qi->tot=etime;

  } else{
   qm->nas_mode = Tuned;
  }

 }
 qm->detect_solver_name();
 //qm->print(1);
 std::string heade_vec = "%%MatrixMarket matrix array real general \n";
 heade_vec += std::to_string(qm->H->ncol);
 heade_vec += " 1\n";
 qm->export_to_file(QPFC->problem_name,heade_vec);
 //std::cout<<"dual res norm: "<<qm->lagrangian_residual_norm()<<"\n";
 //std::cout<<"constraint sat norm: "<<qm->constraint_sat_norm()<<"\n";
 //std::cout<<"primal FB norm: "<<qm->primal_FB_norm()<<"\n";
 //std::cout<<"dual FB norm: "<<qm->dual_FB_norm()<<"\n";
 //std::cout<<"obj: "<<qm->objective<<"\n";
 //std::cout<<QPFC->A_eq->nrow<<"\n";


 std::cout<<qm->sol_name<<",";
 QPFC->print_log();
 std::cout<<num_thread<<",";
 qm->print_log();
 //qm->print_qp_range(1);
/* std::cout<<qm->eps_abs<<";"
          <<qm->num_iter<<";"<<qm->qi->tot<<";"<<num_thread<<";"
          <<qm->active_set.size()<<";;";*/
 //std::cout<<qm->qi->factt<<";";
 //qm->ss->psi->print_profiling();


 delete qm;
 delete QPFC;

}

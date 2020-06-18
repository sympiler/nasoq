//
// Created by kazem on 2020-05-17.
//

#ifndef SCO_CONVERTOR_SMP_FORMAT_H
#define SCO_CONVERTOR_SMP_FORMAT_H

#include <utility>

#include <string>
#include "io.h"
#include "sparse_utilities.h"

namespace format{

 struct  Description{
  std::string ID_; //unique integer
  // belongs to Computer Graphics, Robotics, Finance, MM, ML
  std::string category_;
  // QP name, the group name, the source tool/paper, exporter author, export time.
  std::string name_, group_,source_, author_,date_;
  //Contact simulation, Shape deformation, Model Reconstruction, MPC, ...
  std::string application_;
  Description():ID_("00000"),category_("computer graphics"),name_(""),
  group_(""),source_(""),author_("Kazem Cheshmi"),date_("06/20"),
  application_("contact simulation"){}

  std::string get_desc(){
   std::string out="Sparse Mathematical Programming Repository\n";
   out += ("ID = "+ ID_ + "\n"); out += ("category = "+ category_+ "\n");
   out += ("application = "+ application_ + "\n");
   out += ("name = "+ name_ + "\n"); out += ("group = "+ group_ + "\n");
   out += ("source = "+ source_ + "\n"); out += ("author = "+ author_ + "\n");
   out += ("date = "+ date_ + "\n");
   return out;
  }

  void set_desc_from_string(std::string d){
   for (unsigned i=0; i<d.length(); d[i]=tolower(d[i]),i++);
   std::stringstream ss(d);
   std::string line;
   while (std::getline(ss, line, '\n')){
    auto p1 = line.find('=');
    if(p1 != std::string::npos){
     std::string key = line.substr(0, p1);
     std::string value = line.substr(p1+1);
     trim(key);
     trim(value);
     if(key == "id")
      ID_ = value;
     else if(key == "category")
      category_ = value;
     else if(key == "application")
      application_ = value;
     else if(key =="name")
      name_ = value;
     else if(key == "group")
      group_= value;
     else if(key == "source")
      group_ = value;
     else if(key == "author")
      author_ = value;
     else if(key == "date")
      date_ = value;
    }
   }
  }


 };

 /*
  * Storage for mat for Sparse Constraint Optimization
  * Min 1/2 x^T H x + q^T x + r
  *      A x = b
  *      l <= C x <= u
  */
 struct SMP{
  int num_vars_;
  std::string in_path_;
  std::string desc_;
  Description desc_struct_;
  // Input attributes
  CSC *H_, *A_, *C_;
  CSC *AT_, *CT_;
  Dense  *q_, *b_, *l_, *u_;
  double r_;

  // Output attributes
  Dense *primals_, *duals_;
  double optimal_obj_;

  /*
   * Constructor for reading from file
   */
  explicit  SMP(std::string path):in_path_(std::move(path)), H_(NULLPNTR),
  A_(NULLPNTR), C_(NULLPNTR), q_(NULLPNTR),
  b_(NULLPNTR), l_(NULLPNTR), u_(NULLPNTR), r_(0),AT_(NULLPNTR),CT_(NULLPNTR),
  primals_(NULLPNTR),duals_(NULLPNTR), optimal_obj_(0),num_vars_(0),desc_("\n"){
  }


  /*
   * constructor with only inputs, no file
   */
  SMP(CSC *o, Dense *q, CSC *A, CSC *C, Dense *b, Dense *l, Dense *u,
      double r, std::string  desc)
    :H_(o), A_(A), C_(C), q_(q), b_(b), l_(l), u_(u), r_(r),
    AT_(NULLPNTR),CT_(NULLPNTR),
     primals_(NULLPNTR),duals_(NULLPNTR), optimal_obj_(0),desc_(std::move(desc)),
     in_path_(""),num_vars_(0){
   set_num_vars();
   desc_struct_.set_desc_from_string(desc_);
  };


  SMP(const SMP *smp):SMP(""){
   if(smp){
    num_vars_ = smp->num_vars_;
    in_path_ = smp->in_path_;
    desc_ = smp->desc_;
    desc_struct_.set_desc_from_string(desc_);
    H_ = sym_lib::copy_sparse(smp->H_);
    A_ = sym_lib::copy_sparse(smp->A_);
    C_ = sym_lib::copy_sparse(smp->C_);
    AT_ = sym_lib::copy_sparse(smp->AT_);
    CT_ = sym_lib::copy_sparse(smp->CT_);
    b_ = sym_lib::copy_dense(smp->b_);
    q_ = sym_lib::copy_dense(smp->q_);
    l_ = sym_lib::copy_dense(smp->l_);
    u_ = sym_lib::copy_dense(smp->u_);
    primals_ = sym_lib::copy_dense(smp->primals_);
    duals_ = sym_lib::copy_dense(smp->duals_);
    r_ = smp->r_;
   }
  }


  ~SMP(){
   delete H_;
   delete q_;
   delete A_;
   delete b_;
   delete C_;
   delete l_;
   delete u_;
   delete AT_;
   delete CT_;
   delete primals_;
   delete duals_;
  }

  int N_EQ(){
   return A_ ? A_->m : 0;
  }
  int N_INEQ(){
   return C_ ? C_->m : 0;
  }
  int set_num_vars(){
   num_vars_ = H_ ? H_->n : 0;
   num_vars_ = (q_ && num_vars_ == 0) ? q_->row : num_vars_;
   num_vars_ = (A_ && num_vars_ == 0) ? A_->n : num_vars_;
   num_vars_ = (C_ && num_vars_ == 0) ? C_->n : num_vars_;
   num_vars_ = (primals_ && num_vars_ == 0) ? primals_->row : num_vars_;
   num_vars_ = (duals_ && num_vars_ == 0) ? duals_->row : num_vars_;
   return num_vars_;
  }

  bool check(){
   return !(H_->m != H_->n || H_->n != A_->n || H_->n != A_->n || H_->n != C_->n
            || l_->row != u_->row || l_->row != C_->m || A_->m != b_->row);
  }

  bool load(){
   std::ifstream fin(in_path_);
   if(fin.is_open()){
    while (!fin.eof()){
     std::string line;
     std::getline(fin, line);
     for (auto i=0; i<line.length(); line[i]=tolower(line[i]),i++);
     trim(line);
     if(line == "\"quadratic\": |" && !H_) {
      read_mtx_csc_real(fin, H_, true);
     } else if(line == "\"linear\": |" && !q_) {
      read_mtx_array_real(fin, q_);
     } else if(line == "\"equality\": |" && !A_) {
      read_mtx_csc_real(fin, A_, false);
      //smp->A_=A;
     } else if(line == "\"equality bounds\": |" && !b_) {
      read_mtx_array_real(fin, b_);
      //smp->b_=b;
     } else if(line == "\"inequality l-bounds\": |" && !l_) {
      read_mtx_array_real(fin, l_);
      //smp->l_=l;
     } else if(line == "\"inequality u-bounds\": |" && !u_) {
      read_mtx_array_real(fin, u_);
      //smp->u_=u;
     } else if(line == "\"inequality\": |" && !C_) {
      read_mtx_csc_real(fin, C_, false);
      //smp->C_=C;
     }else if(line == "\"fixed\": |") {
      read_real_constant(fin, r_);
     } else if(line == "\"primals\": |" && !primals_) {
      read_mtx_array_real(fin, primals_);
      //smp->primals_ = primals;
     } else if(line == "\"duals\": |" && !duals_) {
      read_mtx_array_real(fin, duals_);
      //smp->duals_=duals;
     } else if(line == "\"optimal objective\": |") {
      read_real_constant(fin, optimal_obj_);
     }else if(line == "\"description\": |" || line == "\"description\":") {
      read_string(fin, desc_);
      desc_struct_.set_desc_from_string(desc_);
     }else if(line.find(':') != std::string::npos){
      std::cout<<"Key in "<<line <<" is unknown\n"; fin.close();
      return false;
     }
    }
    if(set_num_vars() == 0){
     std::cout<<"Warning: Input QP is all zeros \n";
    }
   } else{
    std::cout<<"input path:"<< in_path_ <<"does not exist. \n";
    return false;
   }
   fin.close();
   if(desc_struct_.name_ == ""){//SMP with empty descriptor
    desc_struct_.name_ = strip_name(in_path_);
   }
   return true;
  }


  bool write(const std::string& out_path){
   std::ofstream fout(out_path, std::ios::out);
   if(fout.is_open()){

    fout << "\"Description\": |\n";
    print_string(desc_,fout.rdbuf());

    if(duals_){
     fout << "\"Duals\": |\n";
     print_dense(duals_->row, duals_->col, duals_->lda, duals_->a, fout.rdbuf());
    }

    if(A_){
     fout << "\"Equality\": |\n";
     print_csc(A_->m, A_->n, A_->p, A_->i, A_->x, fout.rdbuf());
    }

    if(b_){
     fout << "\"Equality bounds\": |\n";
     print_dense(b_->row, b_->col, b_->lda, b_->a, fout.rdbuf());
    }

    fout << "\"Fixed\": |\n";
    print_constant(r_,fout.rdbuf());

    if(C_){
     fout << "\"Inequality\": |\n";
     print_csc(C_->m, C_->n, C_->p, C_->i, C_->x, fout.rdbuf());
    }

    if(l_){
     fout << "\"Inequality l-bounds\": |\n";
     print_dense(l_->row, l_->col, l_->lda, l_->a, fout.rdbuf());
    }

    if(u_){
     fout << "\"Inequality u-bounds\": |\n";
     print_dense(u_->row, u_->col, u_->lda, u_->a, fout.rdbuf());
    }

    if(q_){
     fout << "\"Linear\": |\n";
     print_dense(q_->row, q_->col, q_->lda, q_->a, fout.rdbuf());
    }

    if(primals_ || duals_){
     fout << "\"Optimal objective\": |\n";
     print_constant(optimal_obj_,fout.rdbuf());
    }

    if(primals_){
     fout << "\"Primals\": |\n";
     print_dense(primals_->row, primals_->col, primals_->lda, primals_->a,
                 fout.rdbuf());
    }

    if(H_){
     fout << "\"Quadratic\": |\n";
     print_csc(H_->m, H_->n, H_->p, H_->i, H_->x, fout.rdbuf());
    }

   } else{
    std::cout<<"The path: "<< out_path <<" is not available for writing.\n";
    return false;
   }
   fout.close();
   return true;
  }

  bool equality_check(const SMP* smp, bool is_out = false){
   auto infinity_two_vector = [](Dense *v1, Dense *v2){
    if(!v1 && v2)
     return !v2->is_finite();
    else if(v1 && !v2)
     return !v1->is_finite();
    return sym_lib::are_equal(v1, v2);
   };

   bool h_c = sym_lib::are_equal(H_, smp->H_);
   auto a_c = sym_lib::are_equal(A_, smp->A_);
   auto c_c = sym_lib::are_equal(C_, smp->C_);
   auto l_c = infinity_two_vector(l_, smp->l_);
   auto u_c = infinity_two_vector(u_, smp->u_);
   auto q_c = infinity_two_vector(q_, smp->q_);
   auto b_c = infinity_two_vector(b_, smp->b_);
   bool p_c = true, d_c = true, o_c = true;
   if(is_out){
    p_c = sym_lib::are_equal(primals_, smp->primals_);
    d_c = sym_lib::are_equal(duals_, smp->duals_);
    o_c = is_equal(optimal_obj_, smp->optimal_obj_);
   }
   return h_c && a_c && c_c  && l_c && u_c && q_c &&
   b_c && b_c && p_c && d_c && o_c;
  }

  void set_description(std::string d){
   desc_ = d;
  desc_struct_.set_desc_from_string(desc_);}

 };
}
#endif //SCO_CONVERTOR_SMP_FORMAT_H

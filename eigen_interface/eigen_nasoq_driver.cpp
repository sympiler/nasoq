//
// Created by kazem on 2020-05-09.
//

#include <unsupported/Eigen/SparseExtra>
#include <Util.h>
#include "nasoq_eigen.h"


int main(int argc, char *argv[]){

 std::map<std::string,std::string> qp_args;
 parse_args(argc, argv, qp_args);

 /// Reading input matrices.
 Eigen::SparseMatrix<double,Eigen::ColMajor,int> H, A, C;
 Eigen::VectorXd q, b, d;
 std::string message = "Could not load ";
 if( !Eigen::loadMarket( H, qp_args["H"] ) ){ std::cout<<message<<"H"; return 1; }
 if( !Eigen::loadMarketVector( q, qp_args["q"] ) ){ std::cout<<message<<"q"; return 1; }
 if( !Eigen::loadMarket( A, qp_args["A"] ) ){ std::cout<<message<<"A"; return 1; }
 if( !Eigen::loadMarketVector( b, qp_args["b"] ) ){ std::cout<<message<<"b"; return 1; }
 if( !Eigen::loadMarket( C, qp_args["C"] ) ){ std::cout<<message<<"C"; return 1; }
 if( !Eigen::loadMarketVector( d, qp_args["d"] ) ){ std::cout<<message<<"d"; return 1; }

 /// New settings if provided
 int mode = std::stoi(qp_args["mode"]);
 std::string nasoq_mode = qp_args["nasoq"];
 double pert = pow(10, -std::stoi(qp_args["perturbation"]) );
 int iter = std::stoi(qp_args["iterations"]);
 double eps = pow(10, -std::stoi(qp_args["epsilon"]) );
 double tol = pow(10, -std::stoi(qp_args["tolerance"]) );
 auto *qs = new QPSettings();
 qs->eps=eps; qs->reg_diag=pert; qs->inner_iter_ref=qs->outer_iter_ref=iter;
 qs->tol_ref=tol; qs->nasoq_mode=nasoq_mode;

 /// output vectors
 Eigen::VectorXd x, y, z;

 /// call the wrapper.
 int ret = eigen::nasoq::quadprog(H,q,A,b,C,d,x,y,z,qs);

 //std::cout<<" -->"<< ret<<"\n";
 delete qs;
 return 0;
}
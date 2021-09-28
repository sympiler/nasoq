//
// Created by kazem on 2020-05-18.
//

#include <unsupported/Eigen/SparseExtra>
#include <nasoq/common/Util.h>
#include <nasoq/lbl_eigen.h>

using namespace nasoq;

int main(int argc, const char *argv[]){

 std::map<std::string,std::string> qp_args;
 parse_args(argc, argv, qp_args);

 /// Reading input matrices.
 Eigen::SparseMatrix<double,Eigen::ColMajor,int> H;
 Eigen::VectorXd q;
 std::string message = "Could not load ";
 if( !Eigen::loadMarket( H, qp_args["H"] ) ){ std::cout<<message<<"H"; return 1; }
 if( !Eigen::loadMarketVector( q, qp_args["q"] ) ){ std::cout<<message<<"q"; return 1; }
 // WARNINGL for now, if there is a zero diagonal in H,
 // the row/col indices of zero diagonal should be there.
 //TODO: getting more parameters here and replace qs

 /// output vectors
 Eigen::VectorXd x;

 /// call the wrapper.
 int ret = nasoq::linear_solve(H,q,x);

 //print_vec("sol:\n", 0, H.rows(), x.data());
 return ret;
}
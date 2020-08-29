//
// Created by kazem on 2020-05-09.
//

#include <iostream>
#include <map>
#include <cmath>

//#include "qp_format_converter.h"
#include "utils.h"
#include "mp_format_converter.h"


using namespace format;

int main(int argc,const char *argv[]){

 std::map<std::string,std::string> qp_args;
 parse_args(argc, argv, qp_args);
  if(qp_args.find("input") == qp_args.end())
   return 0;
 /*SMP *smp = new SMP(qp_args["H"]);
 smp->load();
 std::cout<<"\n";
 smp->write("tmp.yml");
 delete smp;*/
 Description test;
 std::string dtes=test.get_desc();
 auto *qfc = new QPFormatConverter();
 qfc->load_smp(qp_args["input"]);
//qfc->smp_->write("tmp2.yml");
qfc->smp_to_bounded();

 auto *qfc2 = new QPFormatConverter(qfc->bf_);
//qfc->smp_->write("tmp2.yml");
 qfc2->bounded_to_smp();

 if(!sym_lib::are_equal(qfc->smp_, qfc2->smp_) )
  std::cout<<"WRONG conversion\n";

 if(!sym_lib::are_equal(qfc->bf_, qfc2->bf_) )
  std::cout<<"WRONG conversion\n";
 delete qfc;
 delete qfc2;


 auto *qfc3 = new QPFormatConverter();
 qfc3->load_smp(qp_args["input"]);
//qfc->smp_->write("tmp2.yml");
 qfc3->smp_to_ie();

 auto *qfc4 = new QPFormatConverter(qfc3->ief_);
//qfc->smp_->write("tmp2.yml");
 qfc4->ie_to_smp();

// if(!sym_lib::are_equal(qfc3->smp_, qfc4->smp_) )
//  std::cout<<"WRONG IE conversion\n";
//
// if(!sym_lib::are_equal(qfc3->ief_, qfc4->ief_) )
//  std::cout<<"WRONG conversion\n";

 qfc3->smp_->set_description(dtes);
 qfc3->smp_->write(qp_args["output"]);

 delete qfc3;
 delete qfc4;
 /// Reading input matrices.
// Eigen::SparseMatrix<double,Eigen::ColMajor,int> H, A, C;
// Eigen::VectorXd q, b, d;
// std::string message = "Could not load ";
// if( !Eigen::loadMarket( H, qp_args["H"] ) ){ std::cout<<message<<"H"; return 1; }
// if( !Eigen::loadMarketVector( q, qp_args["q"] ) ){ std::cout<<message<<"q"; return 1; }
// if( !Eigen::loadMarket( A, qp_args["A"] ) ){ std::cout<<message<<"A"; return 1; }
// if( !Eigen::loadMarketVector( b, qp_args["b"] ) ){ std::cout<<message<<"b"; return 1; }
// if( !Eigen::loadMarket( C, qp_args["C"] ) ){ std::cout<<message<<"C"; return 1; }
// if( !Eigen::loadMarketVector( d, qp_args["d"] ) ){ std::cout<<message<<"d"; return 1; }

 return 0;
}
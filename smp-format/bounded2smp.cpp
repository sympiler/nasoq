//
// Created by kazem on 6/12/20.
//

#include "utils.h"
#include "mp_format_converter.h"


using namespace format;

int main(int argc, const char *argv[]){

 std::map<std::string,std::string> qp_args;
 parse_args_bounded(argc, argv, qp_args);

 std::string p_name, q_name, l_name, a_name, u_name;
 std::string output = "noname.yml";

 if(qp_args.find("quadratic") != qp_args.end())
  p_name = qp_args["quadratic"];
 if(qp_args.find("linear") != qp_args.end())
  q_name = qp_args["linear"];
 if(qp_args.find("l-bounds") != qp_args.end())
  l_name = qp_args["l-bounds"];
 if(qp_args.find("constraints") != qp_args.end())
  a_name = qp_args["constraints"];
 if(qp_args.find("u-bounds") != qp_args.end())
  u_name = qp_args["u-bounds"];
 if(qp_args.find("output") != qp_args.end())
  output = qp_args["output"];

 auto *bf = load_bounded(p_name, q_name, l_name, a_name, u_name);
 Description test;
 std::string dtes=test.get_desc();
 auto *qfc = new QPFormatConverter(bf);
//qfc->smp_->write("tmp2.yml");
 qfc->bounded_to_smp();
 qfc->smp_->write(output);

 delete bf;
 delete qfc;

 return 0;
}

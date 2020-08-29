//
// Created by kazem on 6/12/20.
//


#include "utils.h"
#include "mp_format_converter.h"


using namespace format;

int main(int argc, const char *argv[]){

 std::map<std::string,std::string> qp_args;
 parse_args_ie(argc, argv, qp_args);

 std::string p_name, q_name, a_name, b_name, c_name, d_name;
 std::string output = "noname.yml";

 if(qp_args.find("quadratic") != qp_args.end())
  p_name = qp_args["quadratic"];
 if(qp_args.find("linear") != qp_args.end())
  q_name = qp_args["linear"];
 if(qp_args.find("equality bounds") != qp_args.end())
  b_name = qp_args["equality bounds"];
 if(qp_args.find("equalities") != qp_args.end())
  a_name = qp_args["equalities"];
 if(qp_args.find("inequalities") != qp_args.end())
  c_name = qp_args["inequalities"];
 if(qp_args.find("inequality bounds") != qp_args.end())
  d_name = qp_args["inequality bounds"];
 if(qp_args.find("output") != qp_args.end())
  output = qp_args["output"];

 auto *ie = load_ie(p_name, q_name, a_name,b_name, c_name, d_name);
 Description test;
 std::string dtes=test.get_desc();
 auto *qfc = new QPFormatConverter(ie);
//qfc->smp_->write("tmp2.yml");
 qfc->ie_to_smp();
 qfc->smp_->write(output);

 delete ie;
 delete qfc;

 return 0;
}

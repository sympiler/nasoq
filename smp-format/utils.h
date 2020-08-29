//
// Created by kazem on 2020-05-16.
//

#ifndef SCO_CONVERTOR_UTILS_H
#define SCO_CONVERTOR_UTILS_H

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <map>
#include <iostream>
#include "cxxopts.hpp"

namespace format{


bool parse_args_bounded(int argc, const char *argv[], std::map<std::string, std::string> &qp_args) {
  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  try
  {
   cxxopts::Options options(argv[0], " - example command line options");
   options
     .positional_help("[optional args]")
     .show_positional_help();

   options
     .allow_unrecognised_options()
     .add_options()
       ("p,quadratic", "Quadratic matrix",  cxxopts::value<std::string>())
       ("q,linear", "linear vector",  cxxopts::value<std::string>())
       ("a,constraints", "Constraint matrix",  cxxopts::value<std::string>() )
       ("l,lbounds", "Lower bounds", cxxopts::value<std::string>())
       ("u,ubounds", "Upper bounds", cxxopts::value<std::string>())
       ("o,output", "Output path of SMP", cxxopts::value<std::string>())
       ("help", "Print help")
     ;

   auto result = options.parse(argc, argv);
   if(result.count("p"))
    qp_args.insert(std::pair<std::string, std::string>("quadratic",
      result["p"].as<std::string>()));

   if(result.count("q"))
    qp_args.insert(std::pair<std::string, std::string>("linear",
      result["q"].as<std::string>()));

   if(result.count("a"))
    qp_args.insert(std::pair<std::string, std::string>("constraints",
      result["a"].as<std::string>()));

   if(result.count("l"))
    qp_args.insert(std::pair<std::string, std::string>("l-bounds",
                                                       result["l"].as<std::string>()));

   if(result.count("u"))
    qp_args.insert(std::pair<std::string, std::string>("u-bounds",
                                                       result["u"].as<std::string>()));

   if(result.count("o"))
    qp_args.insert(std::pair<std::string, std::string>("output",
                                                       result["o"].as<std::string>()));

   if (result.count("h"))
   {
    std::cout << "needs to write it! :|" << std::endl;
    exit(0);
   }
  }
  catch (const cxxopts::OptionException& e)
  {
   std::cout << "error parsing options: " << e.what() << std::endl;
   exit(1);
  }
  return true;
 }


/*
 void parse_args_bounded(int argc, char **argv, std::map<std::string, std::string> &qp_args){
  const char* const short_opts = "p:q:a:l:u:o:h";
  const option long_opts[] = {
    {"quadratic", required_argument, nullptr, 'p'},
    {"linear", required_argument, nullptr, 'q'},
    {"constraints", required_argument, nullptr, 'a'},
    {"lbounds", required_argument, nullptr, 'l'},
    {"ubounds", required_argument, nullptr, 'u'},
    {"output", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };

  while (true)
  {
   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

   if (-1 == opt)
    break;

   switch (opt)
   {
    case 'p':
     qp_args.insert(std::pair<std::string,std::string>("quadratic", optarg));
     break;
    case 'q':
     qp_args.insert(std::pair<std::string,std::string>("linear", optarg));
     break;
    case 'l':
     qp_args.insert(std::pair<std::string,std::string>("l-bounds", optarg));
     break;
    case 'a':
     qp_args.insert(std::pair<std::string,std::string>("constraints", optarg));
     break;
    case 'u':
     qp_args.insert(std::pair<std::string,std::string>("u-bounds", optarg));
     break;
    case 'o':
     qp_args.insert(std::pair<std::string,std::string>("output", optarg));
     break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
     std::cout<<"Help !\n";
     break;
   }
  }
 }
*/


 bool parse_args_ie(int argc, const char *argv[], std::map<std::string, std::string> &qp_args) {
  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  try
  {
   cxxopts::Options options(argv[0], " - example command line options");
   options
     .positional_help("[optional args]")
     .show_positional_help();

   options
     .allow_unrecognised_options()
     .add_options()
       ("p,quadratic", "Quadratic matrix",  cxxopts::value<std::string>())
       ("q,linear", "linear vector",  cxxopts::value<std::string>())
       ("a,eq", "Eq constraint matrix",  cxxopts::value<std::string>() )
       ("b,eqbounds", "Eq bounds", cxxopts::value<std::string>())
       ("c,ineq", "Inequality matrix", cxxopts::value<std::string>())
       ("d,ineqbounds", "Ineq bounds", cxxopts::value<std::string>())
       ("o,output", "Output path of SMP", cxxopts::value<std::string>())
       ("help", "Print help")
     ;

   auto result = options.parse(argc, argv);
   if(result.count("p"))
    qp_args.insert(std::pair<std::string, std::string>("quadratic",
      result["p"].as<std::string>()));

   if(result.count("q"))
    qp_args.insert(std::pair<std::string, std::string>("linear",
      result["q"].as<std::string>()));

   if(result.count("a"))
    qp_args.insert(std::pair<std::string, std::string>("equalities",
      result["a"].as<std::string>()));

   if(result.count("b"))
    qp_args.insert(std::pair<std::string, std::string>("equality bounds",
      result["b"].as<std::string>()));

   if(result.count("c"))
    qp_args.insert(std::pair<std::string, std::string>("inequalities",
      result["c"].as<std::string>()));

   if(result.count("d"))
    qp_args.insert(std::pair<std::string, std::string>("inequality bounds",
      result["d"].as<std::string>()));

   if(result.count("o"))
    qp_args.insert(std::pair<std::string, std::string>("output",
                                                       result["o"].as<std::string>()));

   if (result.count("h"))
   {
    std::cout << "needs to write it! :|" << std::endl;
    exit(0);
   }
  }
  catch (const cxxopts::OptionException& e)
  {
   std::cout << "error parsing options: " << e.what() << std::endl;
   exit(1);
  }
  return true;
 }


/*
 void parse_args_ie(int argc, char **argv, std::map<std::string, std::string> &qp_args){
  const char* const short_opts = "p:q:a:b:c:d:o:h";
  const option long_opts[] = {
    {"quadratic", required_argument, nullptr, 'p'},
    {"linear", required_argument, nullptr, 'q'},
    {"eq", required_argument, nullptr, 'a'},
    {"eqbounds", required_argument, nullptr, 'b'},
    {"ineq", required_argument, nullptr, 'c'},
    {"ineqbounds", required_argument, nullptr, 'd'},
    {"output", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };

  while (true)
  {
   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

   if (-1 == opt)
    break;

   switch (opt)
   {
    case 'p':
     qp_args.insert(std::pair<std::string,std::string>("quadratic", optarg));
     break;
    case 'q':
     qp_args.insert(std::pair<std::string,std::string>("linear", optarg));
     break;
    case 'a':
     qp_args.insert(std::pair<std::string,std::string>("equalities", optarg));
     break;
    case 'b':
     qp_args.insert(std::pair<std::string,std::string>("equality bounds", optarg));
     break;
    case 'c':
     qp_args.insert(std::pair<std::string,std::string>("inequalities", optarg));
     break;
    case 'd':
     qp_args.insert(std::pair<std::string,std::string>("inequality bounds", optarg));
     break;
    case 'o':
     qp_args.insert(std::pair<std::string,std::string>("output", optarg));
     break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
     std::cout<<"Help !\n";
     break;
   }
  }
 }
*/

 bool parse_args(int argc, const char *argv[], std::map<std::string, std::string> &qp_args) {
  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  try
  {
   cxxopts::Options options(argv[0], " - example command line options");
   options
     .positional_help("[optional args]")
     .show_positional_help();

   options
     .allow_unrecognised_options()
     .add_options()
       ("v,variant", "NASOQ Variant", cxxopts::value<std::string>())
       ("i,input", "Input SMP format", cxxopts::value<std::string>()) //nasoq driver
       ("u,output", "Output SMP format", cxxopts::value<std::string>()) // not imp
       ("d,header", "Print CSV header", cxxopts::value<std::string>()) // only for nasoq driver
       ("o,objective", "Quadratic matrix",  cxxopts::value<std::string>())
       ("l,linear", "linear vector",  cxxopts::value<std::string>())
       ("a,eq", "Equality matrix",  cxxopts::value<std::string>() )
       ("b,eqb", "Equality vector", cxxopts::value<std::string>())
       ("c,ineq", "Inequality matrix", cxxopts::value<std::string>())
       ("g,ineqb", "Inequality vector", cxxopts::value<std::string>())
       ("p,perturb", "Pertubation power value",
        cxxopts::value<std::string>())
       ("r,refinement", "Refinement iterations",
        cxxopts::value<std::string>())
       ("e,epsilon", "Accuracy threshold power",
        cxxopts::value<std::string>())
       ("t,toli", "Iterative refinement threshold",
        cxxopts::value<std::string>())
       ("help", "Print help")
#ifdef CXXOPTS_USE_UNICODE
    ("unicode", u8"A help option with non-ascii: à. Here the size of the"
        " string should be correct")
#endif
     ;

   auto result = options.parse(argc, argv);
   if(result.count("v"))
    qp_args.insert(std::pair<std::string, std::string>("variant",
                                                       result["v"].as<std::string>()));

   if(result.count("i"))
    qp_args.insert(std::pair<std::string, std::string>("input",
                                                       result["i"].as<std::string>()));

   if(result.count("o"))
    qp_args.insert(std::pair<std::string, std::string>("H",
                                                       result["o"].as<std::string>()));

   if(result.count("l"))
    qp_args.insert(std::pair<std::string, std::string>("q",
                                                       result["l"].as<std::string>()));

   if(result.count("a"))
    qp_args.insert(std::pair<std::string, std::string>("A",
                                                       result["a"].as<std::string>()));

   if(result.count("b"))
    qp_args.insert(std::pair<std::string, std::string>("b",
                                                       result["b"].as<std::string>()));

   if(result.count("c"))
    qp_args.insert(std::pair<std::string, std::string>("C",
                                                       result["c"].as<std::string>()));

   if(result.count("g"))
    qp_args.insert(std::pair<std::string, std::string>("d",
                                                       result["g"].as<std::string>()));

   if(result.count("p"))
    qp_args.insert(std::pair<std::string, std::string>("perturbation",
                                                       result["p"].as<std::string>()));

   if(result.count("e"))
    qp_args.insert(std::pair<std::string, std::string>("epsilon",
                                                       result["e"].as<std::string>()));

   if(result.count("r"))
    qp_args.insert(std::pair<std::string, std::string>("iterations",
                                                       result["r"].as<std::string>()));

   if(result.count("t"))
    qp_args.insert(std::pair<std::string, std::string>("tolerance",
                                                       result["t"].as<std::string>()));

   if(result.count("d"))
    qp_args.insert(std::pair<std::string, std::string>("header",
                                                       result["d"].as<std::string>()));

   if(result.count("u"))
    qp_args.insert(std::pair<std::string, std::string>("output",
                                                       result["u"].as<std::string>()));
   if (result.count("h"))
   {
    std::cout << "needs to write it! :|" << std::endl;
    exit(0);
   }
  }
  catch (const cxxopts::OptionException& e)
  {
   std::cout << "error parsing options: " << e.what() << std::endl;
   exit(1);
  }
  return true;
 }

/*

 void parse_args(int argc, char **argv, std::map<std::string, std::string> &qp_args){
  const char* const short_opts = "m:o:l:a:b:c:d:n:p:r:e:t:h";
  const option long_opts[] = {
    {"mode", required_argument, nullptr, 'm'},
    {"objective", required_argument, nullptr, 'o'},
    {"linear", required_argument, nullptr, 'l'},
    {"eq", required_argument, nullptr, 'a'},
    {"eqb", required_argument, nullptr, 'b'},
    {"ineq", required_argument, nullptr, 'c'},
    {"ineqb", required_argument, nullptr, 'd'},
    {"nasoq", required_argument, nullptr, 'n'},
    {"perturb", required_argument, nullptr, 'p'},
    {"refinement", required_argument, nullptr, 'r'},
    {"epsilon", required_argument, nullptr, 'e'},
    {"toli", required_argument, nullptr, 't'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };

  while (true)
  {
   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

   if (-1 == opt)
    break;

   switch (opt)
   {
    case 'm':
     qp_args.insert(std::pair<std::string,std::string>("mode", optarg));
     break;
    case 'o':
     qp_args.insert(std::pair<std::string,std::string>("H", optarg));
     break;
    case 'l':
     qp_args.insert(std::pair<std::string,std::string>("q", optarg));
     break;
    case 'a':
     qp_args.insert(std::pair<std::string,std::string>("A", optarg));
     break;
    case 'b':
     qp_args.insert(std::pair<std::string,std::string>("b", optarg));
     break;
    case 'c':
     qp_args.insert(std::pair<std::string,std::string>("C", optarg));
     break;
    case 'd':
     qp_args.insert(std::pair<std::string,std::string>("d", optarg));
     break;
    case 'n':
     qp_args.insert(std::pair<std::string,std::string>("nasoq", optarg));
     break;
    case 'p':
     qp_args.insert(std::pair<std::string,std::string>("perturbation", optarg));
     break;
    case 'r':
     qp_args.insert(std::pair<std::string,std::string>("iterations", optarg));
     break;
    case 'e':
     qp_args.insert(std::pair<std::string,std::string>("epsilon", optarg));
     break;
    case 't':
     qp_args.insert(std::pair<std::string,std::string>("tolerance", optarg));
     break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
     std::cout<<"Help !\n";
     break;
   }
  }
 }
*/

 // trim from start (in place)
 static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
   return !std::isspace(ch);
  }));
 }

// trim from end (in place)
 static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
   return !std::isspace(ch);
  }).base(), s.end());
 }

// trim from both ends (in place)
 static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
 }

 /// find the name of a file from a path
 std::string strip_name(std::string name){
  auto p1 = name.rfind("/");
  auto p2 = name.rfind(".");
  std::string slide = name.substr( p1 != std::string::npos ? p1+1 : 0,
                                   p2 != std::string::npos ? p2-p1-1 : std::string::npos);
  return slide;
 }




}


#endif //SCO_CONVERTOR_UTILS_H

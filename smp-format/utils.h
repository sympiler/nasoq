//
// Created by kazem on 2020-05-16.
//

#ifndef SCO_CONVERTOR_UTILS_H
#define SCO_CONVERTOR_UTILS_H

#include <getopt.h>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <map>
#include <iostream>

namespace format{

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

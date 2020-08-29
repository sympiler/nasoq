#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

#include "../eigen_interface/nasoq_eigen.h"
using namespace nasoq;
void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  // Reading input matrices.
  Eigen::SparseMatrix<double,Eigen::ColMajor,int> H, A, C;
  Eigen::VectorXd q, b, d;
  Eigen::VectorXd x,y,z;
  mexErrMsgTxt(nrhs>=6,"nrhs should be >= 6");
  const auto parse_sparse  = [&H](
    const mxArray * prhs[], 
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> & S)
  {
    if(mxGetM(prhs[0]) == 0 || mxGetN(prhs[0]) == 0)
    {
      S.resize(0,H.cols());
    }else if(mxIsSparse(prhs[0]))
    {
      igl::matlab::parse_rhs(prhs+0,S);
    }else
    {
      Eigen::MatrixXd Sd;
      igl::matlab::parse_rhs_double(prhs+0,Sd);
      S = Sd.sparseView();
    }
  };
  const auto parse_vector = [](
    const mxArray * prhs[],
    Eigen::VectorXd & V)
  {
    if(mxGetM(prhs[0]) == 0 || mxGetN(prhs[0]) == 0)
    {
      V.resize(0);
    }else
    {
      igl::matlab::parse_rhs_double(prhs,V);
    }
  };
  parse_sparse(prhs+0,H);
  parse_vector(prhs+1,q);
  parse_sparse(prhs+2,A);
  parse_vector(prhs+3,b);
  parse_sparse(prhs+4,C);
  parse_vector(prhs+5,d);
  mexErrMsgTxt(H.cols()==H.rows(),"H should be symmetric");
  mexErrMsgTxt(H.rows()==q.rows(),"H and q should have same # rows");
  mexErrMsgTxt(H.cols()==A.cols(),"H and A should have same # cols");
  mexErrMsgTxt(A.rows()==b.rows(),"A and b should have same # rows");
  mexErrMsgTxt(H.cols()==C.cols(),"H and C should have same # cols");
  mexErrMsgTxt(C.rows()==d.rows(),"C and d should have same # rows");
  auto *qs = new QPSettings();
  {
    int i = 6;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Eps",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        qs->eps = (double)*mxGetPr(prhs[++i]);
        mexErrMsgTxt(qs->eps>0,"Eps should be positive");
      }else if(strcmp("Tol",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        qs->stop_tol = (double)*mxGetPr(prhs[++i]);
        mexErrMsgTxt(qs->stop_tol>0,"Tol should be positive");
      }
//      else if(strcmp("Iter",name) == 0)
//      {
//       validate_arg_double(i,nrhs,prhs,name);
//       validate_arg_scalar(i,nrhs,prhs,name);
//       qs->inner_iter_ref = (double)*mxGetPr(prhs[++i]);
//       qs->nasoq_mode = "predet";
//       mexErrMsgTxt(qs->inner_iter_ref>0,"Iter should be zero or positive");
//      }
//      else if(strcmp("Mode",name) == 0)
//     {
//      validate_arg_double(i,nrhs,prhs,name);
//      validate_arg_scalar(i,nrhs,prhs,name);
//      qs->nasoq_mode = (int)*mxGetPr(prhs[++i]);
//     }
      else{
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }


  Eigen::MatrixXd ret(1,1);
  ret(0,0) = nasoq::quadprog(
    H.triangularView<Eigen::Lower>(),q,A,b,C,d,x,y,z,qs);

  switch(nlhs)
  {
    case 4:
      prepare_lhs_double(z,plhs+3);
    case 3:
      prepare_lhs_double(y,plhs+2);
    case 2:
      prepare_lhs_double(x,plhs+1);
    case 1:
      prepare_lhs_double(ret,plhs+0);
    default:break;
  }
  
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  delete qs;
  return;
}


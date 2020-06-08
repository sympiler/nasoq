% NASOQ Solve a convex quadratic program of the form
%
%   minimize        0.5 x' H x + q' x
%  subject to       A x = b; Cx <= d;
%
% [r,x,y,z] = nasoq(H,q,A,b,C,d)
% [r,x,y,z] = nasoq(H,q,A,b,C,d,'ParameterName',ParameterValue, ...)
% 
%  Inputs:
%    H  n by n sparse Hessian matrix **lower triangle only** (see
%      .triangularView<Eigen::Lower>() )
%    q  n by 1 linear coefficients
%    A  m1 by n linear equality constraint coefficients
%    b  m1 by 1 linear equality upper bounds
%    C  m2 by n linear inequality constraint coefficients
%    d  m2 by 1 linear inequality upper bounds
%    Optional:
%      'Eps'  followed by nasoq's eps paramter {1e-6?}
%      'Tol'  followed by nasoq's tol_def paramter {?}
%
%  Outputs:
%    x  n by 1 solution vector of primal vars
%    y  m1 by 1 solution vector of primal vars
%    z  m2 by 1 solution vector of primal vars
%

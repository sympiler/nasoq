% https://docs.mosek.com/9.2/capi/tutorial-qo-shared.html
%
% min x₁² + 0.1 x₂² + x₃² - x₁x₃ - x₂
% subject to 1 ≤ x₁+x₂+x₃
%            0 ≤ x
%
% or
%
% min ½ [x₁ x₂ x₃] / 2 0   -1 \ /x₁\ + [x₁ x₂ x₃] / 0 \
%                  | 0 0.2  0 | |x₂|              | -1|
%                  \-1 0    2 / \x₃/              \ 0 /
% subject to /-1 -1 -1 \ /x₁\ ≤ /-1\
%            |-1  0  0 | |x₂|   | 0|
%            | 0 -1  0 | \x₃/   | 0|
%            \ 0  0 -1 /        \ 0/
%
H = [2 0 -1;0 0.2 0;-1 0 2];
C = [-1 -1 -1;-eye(3)];
A = [];
b = [];
d = [-1;0;0;0];
q = [0;-1;0];
[r,x,y,z] = nasoq(H,q,A,b,C,d);
assert(r==1,'NASOQ failed');
assert(max(abs(x-[0;5;0]))<eps,'NASOQ gave wrong solution');
fprintf('nasoq wrapper appears to function correctly\n');

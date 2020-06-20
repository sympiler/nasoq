[V,F] = create_regular_grid(100);
L = cotmatrix(V,F);
M = massmatrix(V,F);
H = L'*(M\L);
b = snap_points([0.25 0.25;0.75 0.75],V);
%lx = repmat(0,size(H,1),1);
%ux = repmat(1,size(H,1),1);
d = zeros(2*size(H,1),1);
C = sparse(2*size(H,1),size(H,1));
for i=1:size(H,1)
    d(2*i-1) = 1;
    C(2*i-1,i)=1;
    d(2*i) = 0;
    C(2*i,i)=-1;
end

tic;
epsilon = 1e-3;
[r,x,y,z] = nasoq( ...
  H,zeros(size(H,1),1),sparse(1:numel(b),b,1,numel(b),size(H,1)),[0;1], ...
  C,d,'Eps',epsilon);
toc

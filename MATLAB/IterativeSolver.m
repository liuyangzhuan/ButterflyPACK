clear all
clc 
close all
maxit = 100;
tol = 1e-10;

Z = load('fort.333');
V = load('fort.334');
nn = sqrt(size(Z,1));
NumVect = size(V,1)/nn;

ZZ = zeros(nn,nn);
VV = zeros(nn,NumVect);

iijj=0;
for ii=1:nn
    for jj=1:nn
        iijj = iijj+1;
        ZZ(ii,jj) = Z(iijj,1)+1i*Z(iijj,2);
    end
end

iijj=0;
for ii=1:nn
    for jj=1:NumVect
        iijj = iijj+1;
        VV(ii,jj) = V(iijj,1)+1i*V(iijj,2);
    end
end
[x,flag,relres,iter,resvec] = tfqmr(ZZ,VV(:,1),tol,maxit);
[iter,relres]
clear all;
clc;
close all;

load A.txt
load b.txt


%% 
N = 7165;
ntrain=2000;
diags = 0;
ntest=N-ntrain;
A = reshape(A,N,N);

p = linspace(1,N,N);
p = randperm(N);


A = A(p,p);
b = b(p);

Areg = A(1:ntrain,1:ntrain);
Areg = Areg + diag(ones(ntrain,1)*diags);
Atest = A(1+ntrain:N,1:ntrain);
b1 = b(1:ntrain);
b2 = b(1+ntrain:N);

% norm(b2-Atest*pinv(Areg)*b1,2)/sqrt(ntest)


% x =randn(ntrain,1);
% b1 = Areg*x;
% norm(x-pinv(Areg)*b1,2)/sqrt(ntest)


n1 = ntrain/2;
tol = 1e-4;

A11 = Areg(1:n1,1:n1);
A22 = Areg(1+n1:ntrain,1+n1:ntrain);
A12 = Areg(1:n1,1+n1:ntrain);
A21 = Areg(1+n1:ntrain,1:n1);


A12 = LRcompress(A12,tol);
A21 = LRcompress(A21,tol);


% A11 = [A11(1:n1/2,1:n1/2),LRcompress(A11(1:n1/2,1+n1/2:n1),tol);LRcompress(A11(1+n1/2:n1,1:n1/2),tol),A11(1+n1/2:n1,1+n1/2:n1)];
% A22 = [A22(1:n1/2,1:n1/2),LRcompress(A22(1:n1/2,1+n1/2:n1),tol);LRcompress(A22(1+n1/2:n1,1:n1/2),tol),A22(1+n1/2:n1,1+n1/2:n1)];


% A21 = A12';


Acompress = [A11,A12;A21,A22];
disp(['prediction error: ',num2str(norm(b2-Atest*pinv(Acompress)*b1,2)/sqrt(ntest))])


x =randn(ntrain,1);
b1 = Acompress*x;
disp(['||X_t-H\(H*X_t)||_F/||X_t||_F: ',num2str(norm(x-pinv(Acompress)*b1,2)/sqrt(ntest))])


% A12 = pinv(A11)*A12;
% A21 = pinv(A22)*A21;
% Ahat2 = eye(ntrain);
% Ahat2(1:n1,1+n1:ntrain) = A12;
% Ahat2(1+n1:ntrain,1:n1) = A21;
% Ahat1 = blkdiag(A11,A22);
% Ainv = pinv(Ahat2)*pinv(Ahat1);
% norm(x-Ainv*b1,2)/sqrt(ntest)
% 
% % Adiff = Ainv*Areg-eye(ntrain);
% % 









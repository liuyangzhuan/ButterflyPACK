close all;
clc;
clear all;
addpath('C:\Users\liuya\Desktop\research\codes\eigtool');
mat=load('fort.665');
rmat=mat(:,1);
imat=mat(:,2);
Z = reshape(rmat,500,500)+1i*reshape(imat,500,500);

opts.levels=-4:0;
% opts.ax=[-100 100 -100 100];
eigtool(Z,opts)
eigs(Z,4)
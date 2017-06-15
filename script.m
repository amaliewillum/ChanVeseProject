%%
%Script for running the chan-vese segmentation model
% Amalie Willum 

close all 
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% img: some picture with extension .bmp
% init: the initial contour. Can be: "grid", "square", "circle"
% num_iter: number of total simulations
% j: number of steps between simulation image
% dt: time step
% lamba, mu, nu: scaling parameters from algorithm
% bw: black/white, 0 for off, 1 for on. 
% Noize: apply noise to image. 0 for off, else 'Gauss', 'Speckle', ... 
%   'Salt & Pepper', 'Blur 1', 'Blur 2', 'Blur 3', 'Blur 5', 'Blur 10', ...
%   'Blur 15', 'Blur + Gauss'
% DJ: Dice/Jensen index. 0 for off, 1 for on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    img = 'LPgray.bmp';
    phi0 = init_phi(img,init);  
    num_iter = 5000;
    j = 5000;
    dt = 0.01;
    lambda = 1;
    mu = 0.2;
    nu = 0;
    bw = 1;
    Noize = 0;
    DJ = 1;

cv(img,lambda,lambda,num_iter,mu,nu,dt,bw,j,phi0,Noize,DJ)

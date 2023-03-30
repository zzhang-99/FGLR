%%  Name: FGLR
%   [1]Su X, Zhang Z, Yang F."Fast hyperspectral image denoising and destriping method based on graph
%   Laplacian regularization",IEEE Transactions on Geoscience and Remote Sensing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT NOTE:
% This demo is applicable to the noise case given in the FGLR paper. 
% We have provided an experimental case of the ROSIS Pavia city center data case 2 in the FGLR_demo.
% The default parameters are fit for the noise cases in the paper.
% It is worth noting that when the noise intensity becomes very high, there should be some subtle modification of the parameters.
% For example, when Gaussian noise \sigma = 0.05, impulse noise with the percentage \o = 0.1, 
% the optimal parameters should be: IterMax = 10, \lambda = 0.2, \beta = 4.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;
load Pavia_80.mat;
load Pavia_case2.mat;
[noipsnr, noissim, noimsam] = MSIQA3(OriData3*255, oriData3_noise*255);
tic;
mu = 1;lambda = 0.21; beta = 3.7;
FGLR_outpimg = FGLR(oriData3_noise,mu,lambda,beta);
time(1) = toc;
[FGLRpsnr, FGLRssim, FGLRmsam] = MSIQA3(OriData3*255, FGLR_outpimg*255);
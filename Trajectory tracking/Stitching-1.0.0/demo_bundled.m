% This code implements the method described in paper "Bundled Camera Paths 	for Video Stabilization" (Liu et al. 2013).
% Written by Su Tan
% email: tsu@cse.cuhk.edu.hk

%% Set up vlfeat
cd vlfeat-0.9.19/toolbox;
feval('vl_setup');
cd ../..;
%%

clear;
addpath('mesh');
addpath('RANSAC');
addpath('Path');
meshSize = 3;
lamda = 1;
asaplamda = 3;
maxIte = 20;
stableW = 20;
% More parameters can be modified in "BundledPath.m"

 
bundled = BundledPath('../case14/right.mp4', '../case14/right_stab.mp4', meshSize, lamda, asaplamda, stableW, maxIte);
%bundled.optFoo();
bundled.calcOmega();
bundled.optPath();
bundled.render();
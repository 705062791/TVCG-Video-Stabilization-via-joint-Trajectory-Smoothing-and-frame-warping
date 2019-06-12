%Partial Matlab code for "Bundled Camera Paths for Video Stabilization" (SIGGRAPH 2013)
%Implementation of motion model estimation.
%1. As-similar-as-possible warping.
%2. Local homography estimation on mesh quads.
%require vision tool box for detectSURFFeatures, or you may want to use
%your own features. (N x 2)


clear all;
clc;

addpath('mesh');
addpath('RANSAC');

I1 = imread('D:/img0.png');
I2 = imread('D:/img1.png');
fprintf('detect surf features...');
[I1_features,I2_features]=SURF(I1,I2);
fprintf('[DONE]');

[H,~] = EstimateHomographyByRANSAC(I1_features',I2_features', 0.001);

% H = inv(H);
% H = H./H(3,3);
% 
% H = eye(3);
% H(1, 3) = -100;
H = H./H(3,3);
tform = projective2d(H');
J = imwarp(I1, tform, 'OutputView', imref2d(size(I1)));
figure; imshow(J);
[I1_features,I2_features]=SURF(J,I2);







%%

if length(I1_features) < 20
    error('not enough matched features');
    return;
end

[height,width,~] = size(I1);
%3x3 mesh
quadWidth = width/8;
quadHeight = height/8;
% %4x4 mesh
% quadWidth = width/(2^4);
% quadHeight = height/(2^4);

lamda = 0.2; %mesh more rigid if larger value. [0.2~5]
asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,lamda);
asap.SetControlPts(I1_features,I2_features);%set matched features
asap.Solve();            %solve Ax=b for as similar as possible
homos = asap.CalcHomos();% calc local hommograph transform

gap = 50;
I1warp = asap.Warp(J,gap);                     %warp source image to target image
I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
imshow(I1warpmesh);

gap = 0;
I1warp = asap.Warp(J,gap);                     
imwrite(I1warp,'D:/warp.png');

%access local homography
[h,w,~,~] = size(homos);
for i=1:h-1
    for j=1:w-1
       H(:,:) = homos(i,j,:,:);
       fprintf('Quad=[%d %d]\n',i,j);
       H
    end
end

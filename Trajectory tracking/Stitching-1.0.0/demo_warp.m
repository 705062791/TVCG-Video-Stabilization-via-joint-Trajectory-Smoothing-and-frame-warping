%Partial Matlab code for "Bundled Camera Paths for Video Stabilization" (SIGGRAPH 2013)
%Implementation of motion model estimation.
%1. As-similar-as-possible warping.
%2. Local homography estimation on mesh quads.
%require vision tool box for detectSURFFeatures, or you may want to use
%your own features. (N x 2)

%%
clc;

addpath('mesh');
addpath('RANSAC');
addpath('tan');

I1 = imread('A.png');
I2 = imread('B.png');
%%
[I1_features, I2_features] = SURF(I1, I2);
[height, width, ~] = size(I1);
I1 = insertMarker(I1, I1_features, 'x', 'color', 'red');
I2 = insertMarker(I2, I2_features, 'x', 'color', 'red');
quadWidth = width/8;
quadHeight = height/8;
homos = NewWarping(I1_features, I2_features, height, width, quadHeight, quadWidth, 1);
pa = [I1_features(50, :) 1];
pb = [I2_features(50, :) 1];
col = ceil(pa(1) / quadWidth);
row = ceil(pa(2) / quadHeight);
H = squeeze(homos(row, col, :, :));
pa_ = H * pa';
pa_ = pa_ ./ pa_(3); 
pa_ = pa_';
pb - pa_
% tform = projective2d(H');
% J = imwarp(I1, tform, 'OutputView', imref2d(size(I1)));
figure(1);
imshow(I1);
figure(2);
imshow(I2);
% figure(3);
% imshow(J);


asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,1);
asap.SetControlPts(I1_features,I2_features);%set matched features
asap.Solve();            %solve Ax=b for as similar as possible
homos = asap.CalcHomos();% calc local hommograph transform
H2 = squeeze(homos(row, col, :, :));
pa_ = H2 * pa';
pa_ = pa_ ./ pa_(3); 
pa_ = pa_';
pb - pa_
% tform = projective2d(H2');
% J2 = imwarp(I1, tform, 'OutputView', imref2d(size(I1)));
% J3 = asap.Warp(I1, 1000);
% figure(6);
% imshow(J2);
% figure(7);
% imshow(J3);
return
%%
fprintf('detect surf features...');
[I2_features,I1_features, preH]=SURF(I2,I1);
fprintf('[DONE]');
preH = preH./preH(3,3);
tform = projective2d(preH');
J = imwarp(I2, tform, 'OutputView', imref2d(size(I1)));
figure(3);
imshow(J);
[I1_features,I2_features]=SURF(I1,J);
if length(I1_features) < 20
    error('not enough matched features');
    return;
end
%shFtr(I1, I2, I1_features, I2_features);
[height,width,~] = size(I1);
%1
quadWidth = width/8;
quadHeight = height/8;


%3x3 mesh
%quadWidth = width/(2^3);
%quadHeight = height/(2^3);

%4x4 mesh
% quadWidth = width/(2^4);
% quadHeight = height/(2^4);
 
% temp = [I1_features ones(length(I1_features), 1)] * prehomo;
% I1p_features = temp(:,1:2);


lamda = 0.2; %mesh more rigid if larger value. [0.2~5]
asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,lamda);
asap.ADAPTIVE_WEIGHT = 0;
asap.SetControlPts(I1_features,I2_features);%set matched features
asap.Solve();            %solve Ax=b for as similar as possible
homos = asap.CalcHomos();% calc local hommograph transform
% beta = 0.01;
% error = asap.CalcError(homos, beta)
%%
gap = 50;
I1warp = asap.Warp(I1,gap);                     %warp source image to target image
I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
figure(2);
imshow(I1warpmesh);

gap = 0;
I1warp = asap.Warp(I1,gap);                     
imwrite(I1warp,'WarpResult0.png');

tform = projective2d(inv(preH)');
I1warp = imwarp(I1warp, tform, 'OutputView', imref2d(size(I1)));

imwrite(I1warp,'WarpResult1.png');
%access local homography
[h,w,~,~] = size(homos);
for i=1:h-1
    for j=1:w-1
       H(:,:) = homos(i,j,:,:);
       fprintf('Quad=[%d %d]\n',i,j);
       H
    end
end

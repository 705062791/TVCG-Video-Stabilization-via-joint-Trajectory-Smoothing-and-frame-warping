clc;
clear;
%% load feature trajectories
nFeatures = 3000; % to adjust

nFrames = 500;

folderPath = 'D:/dataset/';
featuretablePath = sprintf('datasets/ret_%d_features.ft', nFeatures);
maskPath = sprintf('datasets/ret_%d_mask.png', nFeatures);

[X, Y, Val] = klt_read_featuretable(featuretablePath);
mask = imread(maskPath);
imshow(mask);

%% show feature trajectories
for i = 1:nFrames
    tmpIdx = sprintf('%03d', i);
    filePath = [folderPath tmpIdx '.png'];
    img = imread(filePath);
    
    pos = [X(:,i), Y(:,i)];
    img = insertMarker(img, pos);
    
    imshow(img);
end
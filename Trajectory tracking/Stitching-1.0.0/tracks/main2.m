clc;
clear;

addpath 'helpers';
epsilon = logspace(-5,3,101); %distortion level

% clean sequence
% [rawData, trueLabels] = load_sequence('arm');
% processedData = process_sequence(rawData, true);
% result = try_sequence('arm', processedData, epsilon);
% computedLabels = find_best_segmentation(result, processedData, 2, epsilon);
% err = compare_labels(trueLabels, computedLabels);

% incomplete sequence
% [rawData, trueLabels, mask] = load_sequence('oc1R2RC');
% processedData = process_sequence(rawData, 'sparse', 'incomplete', mask);
% result = try_sequence('oc1R2RC', processedData, epsilon);
% computedLabels = find_best_segmentation(result, processedData, 3, epsilon);
% err = compare_labels(trueLabels, computedLabels);

% corrupted sequence
% [rawData, trueLabels] = load_sequence('oc1R2RC');
% processedData = process_sequence(rawData, 'sparse', 'corrupted');
% result = try_sequence('oc1R2RC', processedData, epsilon);
% computedLabels = find_best_segmentation(result, processedData, 3, epsilon);
% err = compare_labels(trueLabels, computedLabels);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting
data = 'D:/MATLABs/stitching/case17/';
input_A = 'left/';
input_B = 'right/';
meshSize = 20;
%% get track
trackA = GetTracks([data input_A], meshSize, 400);
save([data 'trackA.mat'], 'trackA');

load([data 'trackA.mat']);
tracks = trackA;
window_size = 40;
window_index = 15;

tracks.setWindowSize(window_size);
[rawData, ID] = tracks.getM(window_index);
mask = (rawData > 0);
mask = double(squeeze(mask(1,:,:) & mask(2,:,:)));



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% complete the motion matrix and segmentation
tic
processedData = process_sequence(rawData, 'sparse', 'incomplete', mask);
result = try_sequence('oc1R2RC', processedData, epsilon);
computedLabels = find_best_segmentation(result, processedData, 5, epsilon);
save('test');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% show the results

load('test');
label_num = max(computedLabels);
feature_num = numel(ID);

colors = {'red', 'green', 'magenta', 'blue', 'cyan', 'yellow'};
shapes = {'o', 'x', '+', '*', 's', 's'};

imgFolder = 'D:/MATLABs/Stitching/case17/left_sd/';
imgSaveFolder = 'D:/dataset/left/';

startIdx = window_size / 2 * (window_index - 1);
endIdx = startIdx + window_size - 1;

% figure;
for idx=startIdx:endIdx
%     disp(idx);
    tmpIdx = sprintf('%03d', idx+1);
    curImgPath = [imgFolder tmpIdx '.png'];
    curImgSavePath = [imgSaveFolder tmpIdx '.png'];
    curImg = imread(curImgPath);
    
    frameIdx = idx - startIdx + 1;
    for labelIndex = 1:label_num        
        pos = rawData(:, computedLabels == labelIndex & mask(:, frameIdx) > 0, frameIdx);
        curImg = insertMarker(curImg, pos', shapes{labelIndex}, 'color', colors{labelIndex});
    end
    imshow(curImg);
    imwrite(curImg, curImgSavePath);
end


function [ path ] = GetPathFromMask( input, MeshSize, tracks, maskPath)
%GETPATH Summary of this function goes here
%   Detailed explanation goes here
    disp(['Computing camera path of ' input]);
    if input(length(input)) ~= '/'
        input = [input '/'] ;
    end
    
    window_size = 40;
    tracks.setWindowSize(window_size);
%     [M, M2, ID] = tracks.getM(10);
    
    
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    fileListMask = dir([maskPath '/*.png']);    
    nFrames = length(fileList);
    if nFrames < 2
        error('Wrong inputs') ;
    end
    path = zeros(nFrames, MeshSize, MeshSize, 3, 3);
    for row = 1:MeshSize
        for col = 1:MeshSize
            path(1, row, col, :, :) = eye(3);
        end
    end
    fprintf('%5d', 1);
    for frameIndex = 2:length(fileList)
        fprintf('%5d', frameIndex);
        if mod(frameIndex, 20) == 0
            fprintf('\n') ;
        end
        fileName = fileList(frameIndex - 1).name;
        I1 = imread([input fileName]);
        fileName = fileList(frameIndex).name;
        I2 = imread([input fileName]);
        [H, W, ~] = size(I1);
        quadH = H / MeshSize;
        quadW = W / MeshSize;
        % TODO: padding
%         asap = AsSimilarAsPossibleWarping(H,W,quadW,quadH, 0.3);
% %         [I1_features,I2_features]=SIFT2(I1,I2);
%         [I1_features,I2_features] = tracks.getGoodF(frameIndex - 1, backList);
%         [I1_features,I2_features] = tracks.getF(frameIndex - 1);
        [I1_features,I2_features] = tracks.getAllF(frameIndex - 1);
        
        % filter by mask
        mask = imread([maskPath '/' fileListMask(frameIndex - 1).name]);
        if size(mask, 3) == 3
            mask = rgb2gray(mask);            
        end
        binmask = mask == 0;
        valid = I1_features(:, 1) > 0.5 & I1_features(:, 2) > 0.5 & I1_features(:, 1) < W - 0.5 & I1_features(:, 2) < H - 0.5...
            & I2_features(:, 1) > 0.5 & I2_features(:, 2) > 0.5 & I2_features(:, 1) < W - 0.5 & I2_features(:, 2) < H - 0.5;
        I1_features = I1_features(valid, :);
        I2_features = I2_features(valid, :);
        valid = binmask(round(I1_features(:, 1) - 1) * H + round(I1_features(:, 2)));
        I1_features = I1_features(valid, :);
        I2_features = I2_features(valid, :);
%         asap.SetControlPts(I1_features,I2_features);
%         asap.Solve();            
%         homos = asap.CalcHomos();
        homos = NewWarping(I1_features, I2_features, H, W, quadH, quadW, 1);
        markedI1 = insertMarker(I1, I1_features, 's', 'color', 'blue');
        markedI2 = insertMarker(I2, I2_features, 's', 'color', 'blue');
%         WarpIm = asap.Warp(markedI1, 100);
%         imwrite(WarpIm, ['../foo_reference/' int2str(frameIndex) '_warp.png'])
        imwrite(markedI1, ['../foo_reference/' int2str(frameIndex) '_1.png'])
        imwrite(markedI2, ['../foo_reference/' int2str(frameIndex) '_2.png'])
        if length(I1_features) < 4
            error('not enough matched features');
        end
        for i = 1:MeshSize
            for j = 1:MeshSize
                path(frameIndex, i, j, :, :) = squeeze(homos(i, j, :, :)) * squeeze(path(frameIndex - 1, i, j, :, :));
                path(frameIndex, i, j, :, :) = squeeze(path(frameIndex, i, j, :, :)) / path(frameIndex, i, j, 3, 3);
            end
        end            
    end
end


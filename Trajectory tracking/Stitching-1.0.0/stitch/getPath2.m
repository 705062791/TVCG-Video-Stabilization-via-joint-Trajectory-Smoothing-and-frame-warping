function [ path, backMask ] = getPath2( input, MeshSize, UseImageSeq )
%GETPATH Summary of this function goes here
%   Detailed explanation goes here
    disp(['Computing camera path of ' input]);
    if UseImageSeq
        if input(length(input)) ~= '/'
            input = [input '/'] ;
        end
        fileList = dir(input);
        fileList = fileList(3:length(fileList));
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
        
        backMask = zeros(200, nFrames);
        backMask(:, 1) = 1;
        
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
            
            if frameIndex > 2
                validMask = backMask(101:200, frameIndex - 1);
            else
                validMask = ones(100, 1);
            end
            
            [homos, maskOut] = getHomo(H, W, quadW, quadH, I1, I2, 1, 1, frameIndex - 1, validMask); 
            backMask(:, frameIndex) = maskOut;
            
%             asap = AsSimilarAsPossibleWarping(H,W,quadW,quadH, 1);
%             [I1_features,I2_features]=SIFT2(I1,I2);
%             asap.SetControlPts(I1_features,I2_features);
%             asap.Solve();            
%             homos = asap.CalcHomos();
            
            % compute C from F
            for i = 1:MeshSize
                for j = 1:MeshSize
                    path(frameIndex, i, j, :, :) = squeeze(homos(i, j, :, :)) * squeeze(path(frameIndex - 1, i, j, :, :));
                    path(frameIndex, i, j, :, :) = squeeze(path(frameIndex, i, j, :, :)) / path(frameIndex, i, j, 3, 3);
                end
            end            
        end
    else
        nFrames = 101;
        path = zeros(nFrmaes, MeshSize, MeshSize, 3, 3); 
    end

end

function [homo, backMask] = getHomo(H, W, quadW, quadH, I1, I2, asaplamda, ADAPTIVE_WEIGHT, frameIndex, maskIn)
    [I1_features, I2_features, score, I1_all, I2_all]=SIFT2(I1,I2,maskIn);
    fprintf('Frame: %d - %d Score: %f\n', frameIndex, frameIndex + 1, score);
    mask1 = getValidMask(H, W, I1_features, I1_all);
    mask2 = getValidMask(H, W, I2_features, I2_all);

    if sum(mask1) < 4
        error('not enough matched features');
    end

    asap = AsSimilarAsPossibleWarping(H, W, quadW, quadH, asaplamda);
    asap.ADAPTIVE_WEIGHT = ADAPTIVE_WEIGHT;
    asap.SetControlPts(I1_features,I2_features);%set matched features
    asap.Solve();            %solve Ax=b for as similar as possible
    qW = floor(W / 10);
    qH = floor(H / 10);
    if sum(maskIn) ~= 101
        markedI1 = drawFeatures(I1, I1_all, 'red');
        markedI2 = drawFeatures(I2, I2_all, 'red');
        markedI1 = drawFeatures(markedI1, I1_features, 'green');
        markedI2 = drawFeatures(markedI2, I2_features, 'green');
        WarpIm = asap.Warp(markedI1, 100);
        for i = 1:10
            for j = 1:10
                if (mask1((i-1)*10 + j) == 1)
                     markedI1 = insertShape(markedI1, 'Rectangle', [(i-1) * qW + 2 (j-1) * qH + 2 qW - 4 qH - 4]);
                end
                if (mask2((i-1)*10 + j) == 1)
                     markedI2 = insertShape(markedI2, 'Rectangle', [(i-1) * qW + 2 (j-1) * qH + 2 qW - 4 qH - 4]);
                end
            end
        end


        imwrite(WarpIm, ['../foo_reference/' int2str(frameIndex) '_warp.png'])
        imwrite(markedI1, ['../foo_reference/' int2str(frameIndex) '_1.png'])
        imwrite(markedI2, ['../foo_reference/' int2str(frameIndex) '_2.png'])
    end
    homo = asap.CalcHomos();            
    backMask(1:100) = mask1;
    backMask(101:200) = mask2;
end

function mask = getValidMask(H, W, features, all)
    [m, n] = size(features);
    if m == 2 && n > 10
        features = features';
    end
    qH = H / 10;
    qW = W / 10;

    index = floor(features(:, 1) / qW) * 10 + (floor(features(:, 2) / qH)) + 1;
    index_all = floor(all(:, 1) / qW) * 10 + (floor(all(:, 2) / qH)) + 1;
    voting = histcounts([index; 1; 100], 100);  
    voting(1) = voting(1) - 1;
    voting(100) = voting(100) - 1;       

    voting_all = histcounts([index_all; 1; 100], 100);  
    voting_all(1) = voting_all(1) - 1;
    voting_all(100) = voting_all(100) - 1;

    voting(voting < voting_all * 0.5) = 0;
    voting(voting < 2) = 0;
    voting(voting > 0) = 1;
    mask = voting;
end

function imout = drawFeatures(im, features, color)
    imout = insertMarker(im, features, 'color', color);
end


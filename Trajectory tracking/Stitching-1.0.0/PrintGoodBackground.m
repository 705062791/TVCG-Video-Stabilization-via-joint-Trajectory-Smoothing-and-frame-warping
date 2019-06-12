function PrintGoodBackground( input, tracks, good, savePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    nFrames = length(fileList);    
    if nFrames ~= tracks.nFrame
        error('frame number not match');
    end
    % pad good with 1 column for label 0
    good = [ones(tracks.nWindow, 1) * (-1) good];
    for frameIndex = 1:nFrames
        fileName = fileList(frameIndex).name;         
        frame = imread([input fileName]);
        windowIndex = floor((frameIndex - 1) / (tracks.wSize / 2));
        hasPoints = tracks.points(:, frameIndex, 1) > 0;
        if windowIndex > 0 && windowIndex < tracks.nWindow - 1            
            backPoint = hasPoints & (good(windowIndex + 1, tracks.labels(:, frameIndex, 2) + 1)' == 1 | good(windowIndex, tracks.labels(:, frameIndex, 1) + 1)' == 1);
            forePoint = hasPoints & (good(windowIndex + 1, tracks.labels(:, frameIndex, 2) + 1)' == 0 & good(windowIndex, tracks.labels(:, frameIndex, 1) + 1)' == 0);
        end
           
        if windowIndex == 0            
            backPoint = hasPoints & good(windowIndex + 1, tracks.labels(:, frameIndex, 2) + 1)' == 1;
            forePoint = hasPoints & good(windowIndex + 1, tracks.labels(:, frameIndex, 2) + 1)' == 0;
        end
        
        if windowIndex == tracks.nWindow            
            backPoint = hasPoints & good(windowIndex, tracks.labels(:, frameIndex, 1) + 1)' == 1;
            forePoint = hasPoints & good(windowIndex, tracks.labels(:, frameIndex, 1) + 1)' == 0;
        end
        backPos = squeeze(tracks.points(backPoint, frameIndex, :));           
        forePos = squeeze(tracks.points(forePoint, frameIndex, :));  
        frame = insertMarker(frame, backPos, 's', 'color', 'green');        
        frame = insertMarker(frame, forePos, 'o', 'color', 'red');        
        imwrite(frame, [savePath int2str(frameIndex) '.png']);
        imshow(frame);
    end

end


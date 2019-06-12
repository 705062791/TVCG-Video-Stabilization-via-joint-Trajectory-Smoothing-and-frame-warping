function PrintBackground( input, tracks, backList, savePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    nFrames = length(fileList);
    if nFrames ~= tracks.nFrame
        error('frame number not match');
    end
    for frameIndex = 1:nFrames
        fileName = fileList(frameIndex).name;         
        frame = imread([input fileName]);        
        hasPoints = tracks.points(:, frameIndex, 1) > 0;
        backPoints = hasPoints & backList == 1;
        forePoints = hasPoints & backList == -1;
        
        backPos = squeeze(tracks.points(backPoints, frameIndex, :));           
        forePos = squeeze(tracks.points(forePoints, frameIndex, :));
        if size(backPos, 2) == 1
            backPos = backPos'; 
        end
        if size(forePos, 2) == 1
            forePos = forePos'; 
        end
        frame = insertMarker(frame, backPos, 's', 'color', 'green');        
        frame = insertMarker(frame, forePos, 's', 'color', 'red');        
        imwrite(frame, [savePath int2str(frameIndex) '.png']);
%         imshow(frame);
    end

end


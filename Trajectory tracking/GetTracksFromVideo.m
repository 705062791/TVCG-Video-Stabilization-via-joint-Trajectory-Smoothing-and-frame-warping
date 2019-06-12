function [ tracks ] = GetTracksFromVideo( trackFile, input, meshSize, demand )
%GETTRACKSFROMVIDEO Summary of this function goes here
%   Detailed explanation goes here
%GetTracks Compute tracks by KLT
%   Use KLT to track evenly fistributed track points
%   input: the path to images
%   meshSize: the meshSize of video stitching

    videoPath = input;
    
    

    xyloObj = VideoReader(videoPath);
    nFrames = xyloObj.NumberOfFrames;
    
    tracks = TrackLib(nFrames);
    
    tracks.videoHeight = xyloObj.Height;
    tracks.videoWidth = xyloObj.Width;
    
    tracker = vision.PointTracker('MaxBidirectionalError', 1);
    
    frame = read(xyloObj, 1); 
    livePoints = getMorePoints(frame, meshSize, 0, [], demand);
    livePoints = filtermask(frame, livePoints);  %
    
    initialize(tracker, livePoints, frame);
    tracks.addPoints(livePoints, 1);
    for frameIndex = 2:nFrames
        fprintf('%5d', frameIndex);
        if mod(frameIndex, 20) == 0
            fprintf('\n') ;
        end        
        frame = read(xyloObj, frameIndex); 
        
        [livePoints, validity] = step(tracker, frame);
        tracks.endPoints(validity, frameIndex);
        tracks.updatePoints(livePoints(validity == true, :), frameIndex);
        
        morePoints = getMorePoints(frame, meshSize, length(tracks.live), livePoints(validity == true, :), demand);
        morePoints = filtermask(frame, morePoints); %
        
        tracks.addPoints(morePoints, frameIndex);
        livePoints = [livePoints(validity == true, :); morePoints];
        setPoints(tracker, livePoints);
        marked = insertMarker(frame, livePoints);
        imshow(marked);
        
    end
    tracks.endPoints(false(length(tracks.live), 1), nFrames + 1);
    
    save(trackFile, 'tracks');
end

function pointsMore = getMorePoints(frame, meshSize, nP, oldpoints, demand)
    demand = demand / (meshSize * meshSize);
    votes = zeros(meshSize);
    [H, W, ~] = size(frame);
    threshold = 0.5;
    if nP > 0
        votes = getVotes(frame, meshSize, oldpoints);
    end
    points = [];
    
    for row = 1:meshSize
        for col = 1:meshSize
            if votes(row, col) < demand * 0.7
                nMore = demand - votes(row, col);
                roi = [1 + (col - 1) * W / meshSize, 1 + (row - 1) * H / meshSize, W / meshSize - 1, H / meshSize - 1];                
                pNew = detectMinEigenFeatures(rgb2gray(frame), 'ROI', roi, 'MinQuality', threshold); 
                while (size(pNew, 1) < nMore) && threshold > 0.1
                    threshold = threshold - 0.1; 
                    threshold = max(threshold, 0);
                    pNew = detectMinEigenFeatures(rgb2gray(frame), 'ROI', roi, 'MinQuality', threshold); 
                end
                if nMore < size(pNew, 1)
                    ordering = randperm(length(pNew));
                    pNew = pNew(ordering);
                    pNew = pNew(1:nMore);
                end
                points = [points; pNew.Location];
            end
        end
    end
    pointsMore = points;
    
end

function votes = getVotes(frame, meshSize, points)
    [H, W, ~] = size(frame);    
    qH = H / meshSize;
    qW = W / meshSize;
    index = floor(points(:, 1) / qW) * meshSize + (floor(points(:, 2) / qH)) + 1; 
    voting = histcounts([index; 1; meshSize*meshSize], meshSize*meshSize);  
    voting(1) = voting(1) - 1;
    voting(meshSize * meshSize) = voting(meshSize*meshSize) - 1;    
    votes = reshape(voting, [meshSize meshSize]);
%     votes = votes';    
end

function points = filtermask(frame, points_)

%     points = points_;
%     return;

    mask = frame(:, :, 1) < 20 & frame(:, :, 2) < 20 & frame(:, :, 3) < 20;
    mask_ = ~mask;
    mask = imgaussfilt(double(mask), 50);
%     mask_ = imgaussfilt(double(mask_), 10);
%     mask = double(mask_ > mask);
    videoH = size(frame, 1);
%     imshow(mask * 100);
    mask(mask > 0.2) = 1;
    mask(mask ~= 1) = 0;
%     imshow(mask);

    if numel(points_) > 0
        valid = mask(round(points_(:, 1) - 1) * videoH + round(points_(:, 2))) == 0;
        points = points_(valid, :);
    else
        points = points_;
    end
end

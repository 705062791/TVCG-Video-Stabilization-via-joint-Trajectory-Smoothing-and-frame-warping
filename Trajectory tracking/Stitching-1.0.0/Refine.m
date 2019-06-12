function idBack = Refine(track_old, backList)
%GETREST Refinement of remaining trajectories    
    track = track_old;
    tempBack = zeros(track.nWindow, track.nTrack);
    trackLength = zeros(track.nWindow, track.nTrack);
    idBack = ones(track.nTrack, 1);
    backList = backList + 1;
    
    maxIte = 10;
    for ite = 1:maxIte
        for windowIndex = 1:track.nWindow
            fprintf('%d\n', windowIndex);
            left = (windowIndex - 1) * track.wSize / 2 + 1;
            right= (windowIndex + 1) * track.wSize / 2;
            homos = zeros(3, 3, track.wSize - 1);
            errors = zeros(track.wSize - 1, 1);

            for frameIndex = left:right - 1               
                if frameIndex < left + track.wSize / 2 - 1 % 1-19
                    hasPoints = track.head <= frameIndex & track.tail > frameIndex & track.labels(:, frameIndex, 2) > 0;
                    backPoints = track.labels(:, frameIndex, 2) == backList(windowIndex) & track.labels(:, frameIndex + 1, 2) == backList(windowIndex) ...
                        & track.points(:, frameIndex, 1) > 0 & track.points(:, frameIndex, 2) > 0 ...
                        & track.points(:, frameIndex + 1, 1) > 0 & track.points(:, frameIndex + 1, 2) > 0;
                    pa = squeeze(track.points(backPoints, frameIndex, :));
                    pb = squeeze(track.points(backPoints, frameIndex + 1, :));                
                    [homos(:, :, frameIndex - left + 1), errors(frameIndex - left + 1)] = getHomo(pa, pb);
                end
                if frameIndex == left + track.wSize / 2 - 1 % 20
                    backPoints = track.labels(:, frameIndex, 2) == backList(windowIndex) & track.labels(:, frameIndex + 1, 1) == backList(windowIndex) ...
                        & track.points(:, frameIndex, 1) > 0 & track.points(:, frameIndex, 2) > 0 ...
                        & track.points(:, frameIndex + 1, 1) > 0 & track.points(:, frameIndex + 1, 2) > 0;
                    hasPoints = track.head <= frameIndex & track.tail > frameIndex & track.labels(:, frameIndex, 2) > 0;
                    pa = squeeze(track.points(backPoints, frameIndex, :));
                    pb = squeeze(track.points(backPoints, frameIndex + 1, :));                
                    [homos(:, :, track.wSize / 2), errors(track.wSize / 2)] = getHomo(pa, pb);                
                end
                if frameIndex > left + track.wSize / 2 - 1 % 21-39
                    backPoints = track.labels(:, frameIndex, 1) == backList(windowIndex) & track.labels(:, frameIndex + 1, 1) == backList(windowIndex) ...
                        & track.points(:, frameIndex, 1) > 0 & track.points(:, frameIndex, 2) > 0 ...
                        & track.points(:, frameIndex + 1, 1) > 0 & track.points(:, frameIndex + 1, 2) > 0; 
                    hasPoints = track.head <= frameIndex & track.tail > frameIndex & track.labels(:, frameIndex, 1) > 0;
                    pa = squeeze(track.points(backPoints, frameIndex, :));
                    pb = squeeze(track.points(backPoints, frameIndex + 1, :));                
                    [homos(:, :, frameIndex - left + 1), errors(frameIndex - left + 1)] = getHomo(pa, pb);
                end            
                trackLength(windowIndex, hasPoints) = 1;
    %             tempBack(windowIndex, backPoints) = 1;
            end

            refError = mean(errors);
            errors = zeros(track.nTrack, 1);
            count = zeros(track.nTrack, 1);
            for frameIndex = left:right - 1
    %             if frameIndex < left + track.wSize / 2 - 1 % 1-19
    %                 backPoints = track.labels(:, frameIndex, 2) == backList(windowIndex) & track.labels(:, frameIndex + 1, 2) == backList(windowIndex) ...
    %                     & track.head <= frameIndex & track.tail > frameIndex                
    %             end
    %             if frameIndex == left + track.wSize / 2 - 1 % 20
    %                 backPoints = track.labels(:, frameIndex, 2) == backList(windowIndex) & track.labels(:, frameIndex + 1, 1) == backList(windowIndex) ...
    %                     & track.points(:, frameIndex, 1) > 0 & track.points(:, frameIndex, 2) > 0 ...
    %                     & track.points(:, frameIndex + 1, 1) > 0 & track.points(:, frameIndex + 1, 2) > 0;                
    %             end
    %             if frameIndex > left + track.wSize / 2 - 1 % 21-39
    %                 backPoints = track.labels(:, frameIndex, 1) == backList(windowIndex) & track.labels(:, frameIndex + 1, 1) == backList(windowIndex) ...
    %                     & track.points(:, frameIndex, 1) > 0 & track.points(:, frameIndex, 2) > 0 ...
    %                     & track.points(:, frameIndex + 1, 1) > 0 & track.points(:, frameIndex + 1, 2) > 0;                 
    %             end
                selected = track.head <= frameIndex & track.tail > frameIndex & (track.labels(:, frameIndex, 1) > 0 | track.labels(:, frameIndex, 2) > 0);
                pa = squeeze(track.points(selected, frameIndex, :));
                pb = squeeze(track.points(selected, frameIndex + 1, :));
                error1 = getError(pa, pb, squeeze(homos(:, :, frameIndex - left + 1)));
                errors(selected) = errors(selected) + error1;
                count(selected) = count(selected) + 1;
            end
            errors(count == 0) = 100;
            errors(count ~=0) = errors(count ~=0) ./ count(count ~=0);
            tempBack(windowIndex, errors < max(2 * refError, 2)) = 1;
        end
    end
    idBack = sum(tempBack) >= ceil(sum(trackLength) * 0.8);
    idBack = idBack';
    
end

function error = getError(pa, pb, H)
    if size(pa, 2) == 1
        pa = pa';
        pb = pb';
    end
    nP = size(pa, 1);
    pa = [pa'; ones(1,nP)];
    pb = [pb'; ones(1,nP)];
    paWarp = H * pa;
    paWarp(1, :) = paWarp(1, :) ./ paWarp(3, :);
    paWarp(2, :) = paWarp(2, :) ./ paWarp(3, :);
    paWarp(3, :) = 1;
    error = sqrt(sum((paWarp - pb) .^2))';    
end

function [H, error] = getHomo(pa, pb)
    pa = [pa'; ones(1,size(pa, 1))];
    pb = [pb'; ones(1,size(pb, 1))];
    [x1, T1] = NormalisePoints(pa);
    [x2, T2] = NormalisePoints(pb);
    H = DLT(x1, x2);
    H = T2\H*T1;
    
    paWarp = H * pa;
    paWarp(1, :) = paWarp(1, :) ./ paWarp(3, :);
    paWarp(2, :) = paWarp(2, :) ./ paWarp(3, :);
    paWarp(3, :) = 1;
    error = mean(sqrt(sum((paWarp - pb) .^2)));    
end

function error = getAvgError(pa, pb, H)
    if size(pa, 2) == 1
        pa = pa';
        pb = pb';
    end
    pa = [pa'; ones(1,size(pa, 1))];
    pb = [pb'; ones(1,size(pb, 1))];
    paWarp = H * pa;
    paWarp(1, :) = paWarp(1, :) ./ paWarp(3, :);
    paWarp(2, :) = paWarp(2, :) ./ paWarp(3, :);
    paWarp(3, :) = 1;
    error = mean(sqrt(sum((paWarp - pb) .^2)));    
end


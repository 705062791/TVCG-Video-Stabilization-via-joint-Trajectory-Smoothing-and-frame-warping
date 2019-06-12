function checkTrack(track)
    for frameIndex = 2:track.nFrame - 1
        bad = track.points(:, frameIndex - 1, 1) > 0 ...
            & track.points(:, frameIndex, 1) == 0 ...
            & track.points(:, frameIndex + 1, 1) > 0;
        
        if sum(bad) > 0
            fprintf('?') ;
        end
    end
end
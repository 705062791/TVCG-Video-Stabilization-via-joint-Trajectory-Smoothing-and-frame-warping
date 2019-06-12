function GetRest( input, name, track, backList)
%GETREST Refinement of remaining trajectories    
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    h = fspecial('average', 11);
    parfor windowIndex = 1:track.nWindow
        fprintf('%d\n', windowIndex);
%         file1 = fopen(['tmp/' name int2str(windowIndex) '.txt'], 'w');
%         file2 = fopen(['tmp/' name int2str(windowIndex) '_id.txt'], 'w');
        [M, id] = track.getM(windowIndex);
        left = (windowIndex - 1) * track.wSize / 2 + 1;
        right= (windowIndex + 1) * track.wSize / 2;
        mask = (M > 0);   mask = double(squeeze(mask(1,:,:) & mask(2,:,:)));
        [~, complete] = process_sequence(M, 'sparse', 'incomplete', mask);
        x = complete(1:2:end, :)'; % N * 40
        y = complete(2:2:end, :)'; % N * 40        
        % print out the track ids
        restMask = track.labels(id, left + track.wSize / 2, 1) ~= backList(windowIndex) + 1;
        id = id(restMask);
        x = x(restMask, :);
        y = y(restMask, :);
        nTrack = sum(restMask);
        rest = zeros(nTrack, 7);
        rest(:, 4) = x(:, 1) ./ track.videoWidth;
        rest(:, 5) = y(:, 1) ./ track.videoHeight;
        rest(:, 6) = x(:, end) ./ track.videoWidth;
        rest(:, 7) = y(:, end) ./ track.videoWidth;
        L = zeros(nTrack, 1);
        U = zeros(nTrack, 1);
        V = zeros(nTrack, 1);
        
        Lcount = L;
        Ucount = U;
        Vcount = V;
        
        for frameIndex = left:right
            fileName = fileList(frameIndex).name;
            frame = imread([input fileName]);
%             frame = colorspace('rgb->luv', frame);
            frame = im2double(frame);
            frame = imfilter(frame, h, 'replicate');
            inWindowIndex = frameIndex - left + 1;
            
            frameL = frame(:, :, 1);
            frameU = frame(:, :, 2);
            frameV = frame(:, :, 3);
            
                        
            valid = x(:, inWindowIndex) > 0.5 & x(:, inWindowIndex) < track.videoWidth - 0.5 ...
                & y(:, inWindowIndex) > 0.5 & y(:, inWindowIndex) < track.videoHeight - 0.5;
            
            L(valid) = L(valid) + frameL(round(x(valid, inWindowIndex) - 1) * track.videoHeight + round(y(valid, inWindowIndex) - 1));
            Lcount(valid) = Lcount(valid) + 1;
            U(valid) = U(valid) + frameU(round(x(valid, inWindowIndex) - 1) * track.videoHeight + round(y(valid, inWindowIndex) - 1));
            Ucount(valid) = Ucount(valid) + 1;
            V(valid) = V(valid) + frameV(round(x(valid, inWindowIndex) - 1) * track.videoHeight + round(y(valid, inWindowIndex) - 1));
            Vcount(valid) = Vcount(valid) + 1;
            
        end
        rest(:, 1) = L ./ Lcount;
        rest(:, 2) = U ./ Ucount;
        rest(:, 3) = V ./ Vcount;
%         fprintf(file1, '%d 7\n', nTrack);
%         for i = 1:nTrack
%             fprintf(file2, '%d\n', id(i));
%             fprintf(file1, '%f %f %f %f %f %f %f\n', ...
%                 rest(i, 1), rest(i, 2), rest(i, 3), rest(i, 4), rest(i, 5), rest(i, 6), rest(i, 7));
%         end
        
        % kmeans
        
        prm.nTrial=4; prm.display=1; prm.outFrac=0.05;
        IDX = kmeans2( rest, 5, prm ); 
        
        
        
%         fclose(file1);
%         fclose(file2);
    end
end


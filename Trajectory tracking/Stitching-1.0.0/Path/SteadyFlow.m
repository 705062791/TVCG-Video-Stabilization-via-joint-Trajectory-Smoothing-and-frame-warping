classdef SteadyFlow <handle

	properties
		
		% resource variables
		source;		% input image sequence
		srcName;	% input file name
		outputPath;	% output file path
		imwarp;		% output image sequence
		
		num_frame;
		videoWidth;
		videoHeight;
		
		%As-similar-as-possible warp parameters
		quadWidth;
		quadHeight;
		gap;
		
		%steadyflow parameters
		span;
		w;
		w_c;
		maxIte;
		data;
		gamma;
		
        %data terms
        of;
        dataC;
		mask;
        
	end
	
	methods (Access = public)
		function obj = SteadyFlow(videoFileName, outputPath)
			obj.srcName = videoFileName;
            obj.outputPath = outputPath;
            if ~exist(outputPath, 'dir')
                mkdir(outputPath);
            end
        end
        
        function init(obj)
            fprintf('Loading frames...\n');
            video = VideoReader(obj.srcName);
            k = 1;
            frame = struct('cdata',zeros(video.Height,video.Width,3,'uint8'), 'colormap',[]);            
            while hasFrame(video);
                fprintf('%4d', k);
                if mod(k, 20) == 1
                    fprintf('\n');
                end
                frame(k).cdata = readFrame(video);
                k = k+1;
            end
            obj.num_frame = k-1;
            disp(['- Total number of frames:' int2str(obj.num_frame)]);
            obj.source = frame; 
            obj.videoHeight = video.Height;
            obj.videoWidth = video.Width;
            
        end
        
        function calcOpticalFlow(obj)
            disp('Calculating Optical Flow...');
            if exist('flow_single.mat', 'file')
                load('flow_single.mat');
            else
                opticalFlow = struct('vx', zeros(obj.videoHeight, obj.videoWidth), 'vy', zeros(obj.videoHeight, obj.videoWidth));
                alpha = 0.012;
                ratio = 0.75;
                minWidth = 20;
                nOuterFPIterations = 7;
                nInnerFPIterations = 1;
                nSORIterations = 30;

                para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
                for k = 1:obj.num_frame - 1
                    fprintf('%d/%d\n', k, obj.num_frame-1);
                    % Things to do here:
                    % estimate global homography before calculating the
                    % optical flow.
                    im1 = obj.source(k).cdata;
                    im2 = obj.source(k+1).cdata;
                    [homo, im1warp] = obj.getSingleHomo(im1, im2);
                    vx1 = zeros(obj.videoHeight, obj.videoWidth);
                    vy1 = zeros(obj.videoHeight, obj.videoWidth);
                    for i = 1:obj.videoHeight
                        for j = 1:obj.videoWidth
                            [x,y] = obj.transform(j,i,homo);
                            vy1(i,j) = y - i;
                            vx1(i,j) = x - j;
                        end
                    end
                    %[vx,vy,~] = Coarse2FineTwoFrames(im1,im2,para);
                    [vx2,vy2,~] = Coarse2FineTwoFrames(im1warp,im2,para);
                    vx3 = vx1 + vx2;
                    vy3 = vy1 + vy2;
                    opticalFlow(k).vx = single(vx3);
                    opticalFlow(k).vy = single(vy3);
                end
                save('flow_single.mat','opticalFlow');
            end
            obj.of = opticalFlow; 
        end
        function calcDataC(obj)
%             obj.dataC = struct();
%             disp('Calculating data C...');
%             obj.C_timeline = cell(obj.videoHeight, obj.videoWidth);
%             [obj.C_timeline{1:obj.videoHeight, 1:obj.videoWidth}] = deal(zeros(obj.num_frame, 2));
%             for frameIndex = 1:obj.num_frame;
%                 fprintf('%4d', frameIndex);
%                 if frameIndex == 1
%                     obj.dataC(frameIndex).x = zeros(obj.videoHeight, obj.videoWidth); 
%                     obj.dataC(frameIndex).y = zeros(obj.videoHeight, obj.videoWidth);
%                 else
%                     obj.dataC(frameIndex).x = obj.dataC(frameIndex - 1).x + obj.of(frameIndex - 1).vx;  
%                     obj.dataC(frameIndex).y = obj.dataC(frameIndex - 1).y + obj.of(frameIndex - 1).vy;
%                     for i = 1:obj.videoHeight
%                         for j = 1:obj.videoWidth
%                             obj.C_timeline{i,j}(frameIndex,1) = obj.dataC(frameIndex).x(i,j);
%                             obj.C_timeline{i,j}(frameIndex,2) = obj.dataC(frameIndex).y(i,j);
%                         end
%                     end
%                 end
%                 
%             end
            obj.dataC = zeros(obj.num_frame, obj.videoHeight, obj.videoWidth, 2);
            for frameIndex = 2:obj.num_frame
                fprintf('%4d', frameIndex);
                if mod(frameIndex, 20) == 1
                    fprintf('\n');
                end
                obj.dataC(frameIndex, :, :, 1) = squeeze(obj.dataC(frameIndex-1, :, :, 1)) + obj.of(frameIndex - 1).vx; 
                obj.dataC(frameIndex, :, :, 2) = squeeze(obj.dataC(frameIndex-1, :, :, 2)) + obj.of(frameIndex - 1).vy; 
            end
        end
        function getMask(obj, Gfilter, e)
            fprintf('\nEstimating Masks...\n')
            obj.mask = ones(obj.num_frame, obj.videoHeight, obj.videoWidth);
            % gaussian filter over time
            lenG = length(Gfilter);
            for i = 1:obj.videoHeight
                for j = 1: obj.videoWidth
                    Cx = squeeze(obj.dataC(:,i,j,1));
                    Cy = squeeze(obj.dataC(:,i,j,2));
                    GCx = conv(Cx, Gfilter);
                    GCy = conv(Cy, Gfilter);
                    Dx = zeros(size(Cx));
                    Dy = zeros(size(Cy));
                    Dx(fix(lenG/2)+1 : length(Cx) - fix(lenG/2)) = ...
                        Cx(fix(lenG/2)+1 : length(Cx) - fix(lenG/2)) ...
                        - GCx(lenG:length(GCx) - lenG + 1);
                    Dy(fix(lenG/2)+1 : length(Cy) - fix(lenG/2)) = ...
                        Cy(fix(lenG/2)+1 : length(Cy) - fix(lenG/2)) ...
                        - GCy(lenG:length(GCy) - lenG + 1);
                    D = Dx.^2 + Dy.^2 - e;
                    D(D<=0) = -1;
                    D(D>0) = 0;
                    D = -D;
                    obj.mask(:,i,j) = D;
                end
            end
            for k = 1:obj.num_frame
                figure(1);imshow(squeeze(obj.mask(k,:,:)));
                pause
            end
            %gaussian filter over space
%             for k = 1:obj.num_frame
%                 fprintf('%4d', k);
%                 if mod(k, 20) == 1
%                     fprintf('\n');
%                 end
%                 Cx = squeeze(obj.dataC(k, :, :, 1));
%                 Cy = squeeze(obj.dataC(k, :, :, 2));
%                 Dx = Cx - imfilter(Cx,Gfilter);
%                 Dy = Cy - imfilter(Cy,Gfilter);
%                 D = Dx.^2 + Dy.^2 - e;
%                 D(D<=0) = -1;
%                 D(D>0) = 0;
%                 D = -D;
%                 obj.mask(k,:,:) = D;
%             end
%             fprintf('\n');
        end
        
        function asap(obj, gridSize, asapLambda)
            for frameIndex = 1:obj.num_frame - 1;
                %show the original optical flow
                flow(:,:,1) = obj.of(frameIndex).vx;
                flow(:,:,2) = obj.of(frameIndex).vy;
                imflow = flowToColor(flow);
                figure(1);imshow(imflow);
                
                %asap
                mask1 = squeeze(obj.mask(frameIndex,:,:));
                [fx, fy] = gradient(mask1);
                boundary = abs(fx)+abs(fy);
                boundary(mask1==0) = 0;
                boundary(1:gridSize,:) = 0;
                boundary(obj.videoHeight - gridSize +1:obj.videoHeight, :) = 0;
                boundary(:,1:gridSize) = 0;
                boundary(:, obj.videoWidth -gridSize + 1:obj.videoWidth) = 0;
                
                [row, col] = find(boundary);
                I1_features = [col row];
                I2_features = I1_features;
                if length(row) < 20
                    continue; 
                end
                for k = 1:length(I1_features)
                     I2_features(k,1) = I1_features(k,1) + obj.of(frameIndex).vy(I1_features(k,2), I1_features(k,1));
                     I2_features(k,2) = I1_features(k,2) + obj.of(frameIndex).vx(I1_features(k,2), I1_features(k,1));
                end
                %shFtr(boundary, boundary, I1_features, I2_features);
                asap = AsSimilarAsPossibleWarping(obj.videoHeight, obj.videoWidth, gridSize, gridSize, asapLambda);
                asap.ADAPTIVE_WEIGHT = 1;
                asap.SetControlPts(I1_features, I2_features);   
                asap.Solve();
                homos = asap.CalcHomos();
                
                gap = 50;
                I1warp = asap.Warp(obj.source(frameIndex).cdata,gap);                     %warp source image to target image
                I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
                figure(2);
                imshow(I1warpmesh);


                [row, col] = find(1-mask1);
                missing = [row col];
                for k = 1:length(missing)
                    x = missing(k,1);
                    y = missing(k,2);
                    homo = squeeze(homos(floor((x-1)/gridSize) + 1,floor((y-1)/gridSize) + 1, :, :));
                    [yy, xx] = obj.transform(y, x, homo);
                    obj.of(frameIndex).vx(x, y) = xx - x;
                    obj.of(frameIndex).vy(x, y) = yy - y;
                end
                
                %show the original optical flow
                flow(:,:,1) = obj.of(frameIndex).vx;
                flow(:,:,2) = obj.of(frameIndex).vy;
                imflow = flowToColor(flow);
                figure(3);imshow(imflow);
                
                %revise dataC
                obj.dataC(frameIndex+1, :, :, 1) = squeeze(obj.dataC(frameIndex, :, :, 1)) + obj.of(frameIndex).vx; 
                obj.dataC(frameIndex+1, :, :, 2) = squeeze(obj.dataC(frameIndex, :, :, 2)) + obj.of(frameIndex).vy;
            end
        end
    end
    
    methods (Access = private)
        function [homo, imwarp] = getSingleHomo(obj, I1, I2)
            [I1_features,I2_features]=SURF(I1,I2);

            if length(I1_features) < 20
                error('not enough matched features');
            end
            
            asap = AsSimilarAsPossibleWarping(obj.videoHeight,obj.videoWidth,obj.videoWidth,obj.videoHeight,1);
            asap.ADAPTIVE_WEIGHT = 1;
            asap.SetControlPts(I1_features,I2_features);
            asap.Solve();          
            homo = asap.CalcHomos(); 
            imwarp = asap.Warp(I1, 0);
            homo = squeeze(homo(1,1,:,:)); 
        end
        function [x,y] = transform(~, xx, yy, B)
            res = B * [xx;yy;1];
            x = res(1)/res(3);
            y = res(2)/res(3);
        end     
    end
end
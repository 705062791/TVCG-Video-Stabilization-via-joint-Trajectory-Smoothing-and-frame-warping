classdef BundledPath <handle
    % This code implements the method described in paper "Bundled Camera Paths for Video Stabilization" (Liu et al. 2013).
    % Written by Su Tan
	% email: tsu@cse.cuhk.edu.hk
    
    
    properties
        % source;     % The raw video
        srcName;    % file name of the input video
        outputPath;
        % imwarp;
        
        num_frame;  % number of frames
        width;      % mesh width
        height;     % mesh height
        videoWidth;
        videoHeight;
        
        quadWidth;
        quadHeight;
        gap;
        maxIte;
        span;
        
        data;
        w; 			% weight
        w_c;
        lamda;      % parameter lamda
        gamma;
        stableW;
        
        asaplamda;
        
        validMask;

        
    end
    
    methods (Access = public)
        function obj = BundledPath(videoFileName, outputPath, meshSize, lamda, asaplamda, stableW, maxIte)
            % Load video file
            video = VideoReader(videoFileName);
            k = 1;
            frame = struct('cdata',zeros(video.Height,video.Width,3,'uint8'), 'colormap',[]);            
            while hasFrame(video)
                fprintf('%5d', k);
                if mod(k, 20) == 0
                    fprintf('\n') ;
                end
                readFrame(video);
                k = k+1;
            end
            obj.num_frame = k-1;
            disp(['total number of frames:' int2str(obj.num_frame)]);
            % obj.source = frame;
            
            % Parameters initialization
            obj.outputPath = outputPath;
            if ~exist(outputPath, 'dir')
                mkdir(outputPath);
            end
            obj.gap = 100;
            obj.stableW = stableW;
            % obj.imwarp = struct('cdata', zeros(video.Height+obj.gap*2,video.Width+obj.gap*2,3,'double'), 'colormap', []);
            obj.srcName = videoFileName;
            obj.width = meshSize*meshSize;
            obj.height = obj.width;
            obj.videoWidth = video.Width;
            obj.videoHeight = video.Height;
            obj.quadHeight = video.Height/obj.height;
            obj.quadWidth = video.Width/obj.width;
            obj.lamda = lamda;
            obj.maxIte = maxIte;
            obj.span = 30;
            obj.asaplamda = asaplamda;
            obj.w = zeros(obj.num_frame, obj.num_frame, obj.height, obj.width);
            obj.w_c = cell(obj.num_frame, obj.height, obj.width);
            obj.gamma = zeros(obj.num_frame, obj.height, obj.width);
            
            obj.validMask = ones(200, obj.num_frame);
            % frame registration
            repath = 1;
            if exist('homos.mat', 'file')
                load('homos.mat'); 
                repath = 0;
                if length(homos) ~= obj.num_frame - 1
                    disp('homos.mat doesn"t match with the video');
                    repath = 1;
                end
            end
            if repath == 1
                images = cell(obj.num_frame, 1);
                homos = cell(obj.num_frame - 1, 1);
                error = zeros(obj.num_frame - 1, 1);
                disp('calculating camera path');
                video = VideoReader(videoFileName);
                I2 = readFrame(video);
                images{1} = I2;
                for k = 1:obj.num_frame - 1
                    fprintf('%5d', k);
                    if mod(k, 20) == 0
                        fprintf('\n') ;
                    end
                    I1 = I2;
                    I2 = readFrame(video);
                    images{k+1} = I2;
%                     [homos{k}, error(k)] = obj.getHomo(I1, I2, asaplamda, 1, k, ones(100, 1));
%                     disp(error(k))
                end
                % find the reference trustable frame
%                 [~, bestFrame] = min(error);
                bestFrame = 1;
                fprintf('Best frame warped: %d\n', bestFrame);
                [homos{bestFrame}, error(bestFrame)] = obj.getHomo(images{bestFrame}, images{bestFrame+1}, asaplamda, 1, bestFrame, ones(100, 1));
                % warp again                
                if bestFrame < obj.num_frame - 1
                    for k = bestFrame + 1:obj.num_frame - 1
                        [homos{k}, error(k)] = obj.getHomo(images{k}, images{k+1}, asaplamda, 1, k, squeeze(obj.validMask(101:200, k - 1))); 
                    end
                end
                if bestFrame > 1
                    for k = bestFrame - 1 : 1
                        [homos{k}, error(k)] = obj.getHomo(images{k}, images{k+1}, asaplamda, 1, k, squeeze(obj.validMask(1:100, k + 1))); 
                    end    
                end                                
                save('homos.mat', 'homos');
            end
            obj.data = struct('F', zeros(obj.height, obj.width, 3, 3), 'C', zeros(obj.height, obj.width, 3, 3),'C_inv', zeros(obj.height, obj.width, 3, 3), 'P', cell(maxIte, 1), 'index', 0);
            
            for i = 1:obj.height
                for j = 1:obj.width
                    obj.data(1).C(i,j,:,:) = eye(3);
                end
            end
            for k = 1:obj.num_frame
                obj.data(k).index = k;
                if k<obj.num_frame
                    obj.data(k).F = homos{k};
                end
                if k > 1
                    for i = 1:obj.height
                        for j = 1:obj.width
                            obj.data(k).C(i,j,:,:) = obj.normH(squeeze(obj.data(k-1).F(i,j,:,:)) * squeeze(obj.data(k-1).C(i,j,:,:)));
                            obj.data(k).C_inv(i,j,:,:) = inv(squeeze(obj.data(k).C(i,j,:,:)));
                            %obj.data(k).C(i,j,:,:) = squeeze(obj.data(k-1).C(i,j,:,:)) * squeeze(obj.data(k-1).F(i,j,:,:));
                        end                        
                    end                    
                end
                obj.data(k).P{1} = obj.data(k).C;
                % obj.imwarp(k).cdata = zeros(video.Height+obj.gap*2,video.Width+obj.gap*2,3,'double');
            end
            
        end
        function calcOmega(obj)
            stableW = obj.stableW;
            for i = 1:obj.height
                disp(['calculating parameters: ' int2str(i) '/' int2str(obj.height) ]);
                for j = 1:obj.width
                    for t = obj.span+1:obj.num_frame-obj.span
                        obj.w_c{t,i,j} = zeros(obj.num_frame, 1);
                        for r = t-obj.span:t+obj.span
                            dC = abs(obj.data(t).C(i,j,1,3) - obj.data(r).C(i,j,1,3)) + abs(obj.data(t).C(i,j,2,3) - obj.data(r).C(i,j,2,3));
                            obj.w(t,r,i,j) = gaussmf(abs(t-r), [10 0]) * gaussmf(dC, [200 0]);
                            if t == r
                                obj.w(t,r,i,j) = 0;
                            end
                            obj.w_c{t,i,j}(r) = obj.w(t,r,i,j);
                        end
                        obj.gamma(t,i,j) = sum(obj.w_c{t,i,j}) * 2 * obj.lamda;
                        if (1 == obj.width) || stableW == 0
                            obj.gamma(t,i,j) = obj.gamma(t,i,j) + 1;
                            continue;
                        end
                        if ((i==1)||(i==obj.height))&&((j==1)||(j==obj.width))
                            obj.gamma(t,i,j) = obj.gamma(t,i,j) + 2 * 3 * stableW + 1;
                        end
                        if ((i==1)||(i==obj.height))&&((j>1)&&(j<obj.width))
                            obj.gamma(t,i,j) = obj.gamma(t,i,j) + 2 * 5 * stableW + 1;
                        end
                        if ((i>1)&&(i<obj.height))&&((j==1)||(j==obj.width))
                            obj.gamma(t,i,j) = obj.gamma(t,i,j) + 2 * 5 * stableW + 1;
                        end
                        if ((i>1)&&(i<obj.height))&&((j>1)&&(j<obj.width))
                            obj.gamma(t,i,j) = obj.gamma(t,i,j) + 2 * 8 * stableW + 1;
                        end
                    end
                end
            end
            
        end
        function optFoo(obj)
            for t = obj.span+1:obj.num_frame-obj.span
                for i = 1:obj.height
                    for j = 1:obj.width
                        obj.data(t).P{obj.maxIte}(i,j,:,:) = obj.data(t+1).C(i,j,:,:);
                    end
                end
            end
        end
        function optPath(obj)
            % Fill up A and b
            % A is like
            % [1 0 0 0 ... 0 0]
            % [0 1 0 0 ... 0 1]
            % .          .
            % .          .
            % .          .
            % [0 0 0 0 ... 0 1]
            % [w -w 0 .. 0 0 0]
            % [w 0 -w 0 .. 0 0]
            % .          .
            % .          .
            % [0 0 .. 0 0 -w w]
            % while b is like
            % [c1 c2 ... cn 0 0 0 0 0 0 0 0 0 0 0 0 0]T
            % calculate the gamma and omega
            % use iteration to solve the problem.
            stableW = obj.stableW;
            
            
            for iteIndex = 2:obj.maxIte
                disp(iteIndex);
                for t = 1:obj.span
                    obj.data(t).P{iteIndex} =  obj.data(t).P{iteIndex-1};
                end
                for t = obj.num_frame - obj.span + 1:obj.num_frame
                    obj.data(t).P{iteIndex} =  obj.data(t).P{iteIndex-1};
                end
                for i = 1:obj.height
                    disp([int2str(i) '/' int2str(obj.height)]);
                    for j = 1:obj.width
                        for t = obj.span+1:obj.num_frame - obj.span 
                            value = squeeze(obj.data(t).C(i,j,:,:));
                            % value = squeeze(obj.data(t).P{iteIndex - 1}(i,j,:,:));
                            for r = t - obj.span : t + obj.span
                                if t ~= r
                                    value = value + 2 * obj.lamda* obj.w(t,r,i,j) * squeeze(obj.data(r).P{iteIndex-1}(i,j,:,:));
                                end
                            end
                            if obj.width == 1 || stableW == 0
                                obj.data(t).P{iteIndex}(i,j,:,:) = value / obj.gamma(t,i,j);
                                continue;
                            end
                            
                            
                            % i=1, j=1
                            if i == 1 && j == 1
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j, :, :)) * squeeze(obj.data(t).C_inv(i+1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j+1, :, :)) * squeeze(obj.data(t).C_inv(i, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=1, j=2:15
                            if i == 1 && j < obj.width && j > 1
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j, :, :)) * squeeze(obj.data(t).C_inv(i+1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j+1, :, :)) * squeeze(obj.data(t).C_inv(i, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j-1, :, :)) * squeeze(obj.data(t).C_inv(i, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=1, j=16
                            if i == 1 && j == obj.width
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j, :, :)) * squeeze(obj.data(t).C_inv(i+1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j-1, :, :)) * squeeze(obj.data(t).C_inv(i, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=2:15, j=1
                            if j == 1 && i < obj.height && i > 1
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j, :, :)) * squeeze(obj.data(t).C_inv(i+1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j, :, :)) * squeeze(obj.data(t).C_inv(i-1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j+1, :, :)) * squeeze(obj.data(t).C_inv(i, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=2:15, j=2:15
                            if i > 1 && i < obj.height && j > 1 && j < obj.width
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j, :, :)) * squeeze(obj.data(t).C_inv(i+1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j, :, :)) * squeeze(obj.data(t).C_inv(i-1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j+1, :, :)) * squeeze(obj.data(t).C_inv(i, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j-1, :, :)) * squeeze(obj.data(t).C_inv(i, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :)); 
                            end
                            % i=2:15, j=16
                            if j == obj.width && i < obj.height && i > 1
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j, :, :)) * squeeze(obj.data(t).C_inv(i+1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j-1, :, :)) * squeeze(obj.data(t).C_inv(i, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j, :, :)) * squeeze(obj.data(t).C_inv(i-1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i+1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i+1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=16, j=1
                            if i == obj.height && j == 1
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j, :, :)) * squeeze(obj.data(t).C_inv(i-1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j+1, :, :)) * squeeze(obj.data(t).C_inv(i, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=16, j=2:15
                            if i == obj.height && j < obj.width && j > 1
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j, :, :)) * squeeze(obj.data(t).C_inv(i-1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j+1, :, :)) * squeeze(obj.data(t).C_inv(i, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j+1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j+1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j-1, :, :)) * squeeze(obj.data(t).C_inv(i, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end
                            % i=16, j=16
                            if i == obj.height && j == obj.width
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j, :, :)) * squeeze(obj.data(t).C_inv(i-1, j, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i, j-1, :, :)) * squeeze(obj.data(t).C_inv(i, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                                value = value + 2 * stableW * squeeze(obj.data(t).P{iteIndex-1}(i-1, j-1, :, :)) * squeeze(obj.data(t).C_inv(i-1, j-1, :, :)) * squeeze(obj.data(t).C(i, j, :, :));
                            end

                            obj.data(t).P{iteIndex}(i,j,:,:) = value / obj.gamma(t,i,j);
                            
%                             if obj.data(t).P{iteIndex}(i,j,3,3) - 1 > 1e-6
%                                 disp(obj.data(t).P{iteIndex}(i,j,3,3) - 1);
%                             end
                        end
                    end
                end
            end
            
            
        end
        
        
        function render(obj)
            video = VideoReader(obj.srcName);
            for i = 1:obj.span
                readFrame(video); 
            end
            for frameIndex = obj.span + 1:obj.num_frame - obj.span
                disp(['rendering: # ' int2str(frameIndex)]);
                source = readFrame(video);
                imwarp = zeros(obj.videoHeight + 2 * obj.gap, obj.videoWidth + 2 * obj.gap, 3);
                for i = 1:obj.height
                    for j = 1:obj.width
                        imwarp = obj.warp1(frameIndex, i, j, source, imwarp); %imwarp
                    end
                end
                % obj.drawFeatures(obj.imwarp(frameIndex).cdata, obj.features(frameIndex));
                imwrite(uint8(imwarp), [obj.outputPath '/' int2str(frameIndex) '.png']);
            end
        end
        
        % Rewrite render
        function render2(obj)
            video = VideoReader(obj.srcName);
            for frameIndex = 1: obj.span
                readFrame(video);
            end
            for frameIndex = obj.span+1: obj.num_frame - obj.span
                % determine mesh
                disp(['rendering: ' int2str(frameIndex - obj.span) '/' int2str(obj.num_frame - 2* obj.span)]);
                src = Mesh(obj.videoHeight, obj.videoWidth, obj.quadWidth, obj.quadHeight);
                des = Mesh(obj.videoHeight, obj.videoWidth, obj.quadWidth, obj.quadHeight);
                %
                for i = 0 : obj.height
                    for j = 0 : obj.width
                        x = i * obj.quadHeight + 1;
                        y = j * obj.quadWidth + 1;
                        if i == 0 && j == 0
                            B11 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j+1, :, :)) ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j+1, :, :));
                            %B11 = B11;
                            
                            [xx11, yy11] = obj.transform(y, x, B11);
                            des.setVertex(i,j, myPoint(xx11, yy11));
                            continue;
                        end
                        if i == 0 && j == obj.width
                            B10 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j, :, :))...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j, :, :)) ;
                            %B10 = B10;
                            
                            [xx10, yy10] = obj.transform(y, x, B10);
                            des.setVertex(i,j, myPoint(xx10, yy10));
                            continue;
                        end
                        if i == 0 && j > 0 && j < obj.width
                            B11 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j+1, :, :));
                            %B11 = B11;
                            B10 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j, :, :));
                            %B10 = B10;
                            [xx11, yy11] = obj.transform(y, x, B11);
                            [xx10, yy10] = obj.transform(y, x, B10);
                            des.setVertex(i,j, myPoint(mean([xx10 xx11]), mean([yy10 yy11])));
                            continue;
                        end
                        if i>0 && i < obj.height && j == 0
                            B11 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j+1, :, :));
                            %B11 = B11;
                            B01 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j+1, :, :));
                            %B01 = B01;
                            [xx11, yy11] = obj.transform(y, x, B11);
                            [xx01, yy01] = obj.transform(y, x, B01);
                            des.setVertex(i,j, myPoint(mean([xx01 xx11]), mean([yy01 yy11])));
                            continue;
                        end
                        if i > 0 && i < obj.height && j > 0 && j < obj.width
                            B11 = squeeze(obj.data(frameIndex).C(i+1, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j+1, :, :));
                            %B11 = B11;
                            B01 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j+1, :, :));
                            %B01 = B01;
                            B00 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j, :, :));
                            %B00 = B00;
                            B10 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j, :, :));
                            %B10 = B10;
                            [xx11, yy11] = obj.transform(y, x, B11);
                            [xx01, yy01] = obj.transform(y, x, B01);
                            [xx10, yy10] = obj.transform(y, x, B10);
                            [xx00, yy00] = obj.transform(y, x, B00);
                            [xx, yy] = obj.mergepoints([xx01, xx11, xx10, xx00], [yy01, yy11, yy10, yy00]);
                            des.setVertex(i,j, myPoint(xx, yy));
                            continue;
                        end
                        if i>0 && i < obj.height && j == obj.width
                            B00 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j, :, :));
                            %B00 = B00;
                            B00 = B00 ./B00(3,3); 
                            B10 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i+1, j, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i+1, j, :, :));
                            %B10 = B10;
                            B10 = B10 ./B10(3,3); 
                            [xx10, yy10] = obj.transform(y, x, B10);
                            [xx00, yy00] = obj.transform(y, x, B00);
                            des.setVertex(i,j, myPoint(mean([xx10 xx00]), mean([yy10 yy00])));
                            continue;
                        end
                        if i == obj.height && j == 0
                            B01 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j+1, :, :));
                            %B01 = B01;
                            B01 = B01 ./B01(3,3); 
                            [xx01, yy01] = obj.transform(y, x, B01);
                            des.setVertex(i,j, myPoint(xx01, yy01));
                            continue;
                        end
                        if i == obj.height && j>0 && j < obj.width
                            B00 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j, :, :));
                            %B00 = B00;
                            B00 = B00 ./B00(3,3); 
                            B01 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j+1, :, :))  ...
                * squeeze(obj.data(frameIndex).C_inv(i, j+1, :, :));
                            %B01 = B01;
                            B01 = B01 ./B01(3,3);
                            [xx00, yy00] = obj.transform(y, x, B00);
                            [xx01, yy01] = obj.transform(y, x, B01);
                            des.setVertex(i,j, myPoint(mean([xx01 xx00]), mean([yy01 yy00])));
                            continue;
                        end
                        if i == obj.height && j == obj.width
                            B00 = squeeze(obj.data(frameIndex).P{obj.maxIte}(i, j, :, :)) ...
                * squeeze(obj.data(frameIndex).C_inv(i, j, :, :));
                            %B00 = B00;
                            B00 = B00 ./B00(3,3); 
                            [xx00, yy00] = obj.transform(y, x, B00); 
                            des.setVertex(i,j, myPoint(xx00, yy00));
                            continue;
                        end
                    end
                end
                Img = readFrame(video);
                imwarp = zeros(obj.videoHeight+obj.gap*2,obj.videoWidth+obj.gap*2,3);

                for i=1:obj.height
                    for j=1:obj.width
%                         disp([i,j]);
                        p0 = src.getVertex(i-1,j-1);
                        p1 = src.getVertex(i-1,j);
                        p2 = src.getVertex(i,j-1);
                        p3 = src.getVertex(i,j);

                        q0 = des.getVertex(i-1,j-1);
                        q1 = des.getVertex(i-1,j);
                        q2 = des.getVertex(i,j-1);
                        q3 = des.getVertex(i,j);

                        qd1 = Quad(p0,p1,p2,p3);
                        qd2 = Quad(q0,q1,q2,q3);
                        imwarp = quadWarp(obj,Img,qd1,qd2, imwarp);
                    end
                end
                
                imwrite(uint8(imwarp), [obj.outputPath '/' int2str(frameIndex) '.png']);
            end
        end
        
        function [xx, yy] = mergepoints(~, x, y)
            d = ones(4, 4) * 1e8;
            for i = 1:3
                for j = i+1:4
                    d(i, j) = (x(i) - x(j))^2 + (y(i) - y(j))^2;
                end
            end
            [minrow, min_i] = min(d);
            [~, min_j] = min(minrow);
            xx = 0.5 * (x(min_i(min_j))+ x(min_j));
            yy = 0.5 * (y(min_i(min_j))+ y(min_j));            
        end
        
        % Rewrite render use asap itself
        % Use feature points as controls. 
        function render3(obj)
            % add gap round
            % obj.imwarp(1).cdata(obj.gap+1:obj.gap+obj.videoHeight, obj.gap+1:obj.gap+obj.videoWidth, :) = double(obj.source(1).cdata);
            video = VideoReader(obj.srcName);
            for i = 1:obj.span
                I2 = readFrame(video); 
            end
            for frameIndex = obj.span+1:obj.num_frame - obj.span
                disp(['rendering: # ' int2str(frameIndex)]);
                I1 = I2;
                I2 = readFrame(video);
                [features, ~] = SURF(I2, I1);
                % [features, ~] = saliencyFilter(features, 0, frameIndex);
                % compute new location of feature point
                features_2 = zeros(size(features));
                for featureIndex = 1:size(features, 1)
                    x = features(featureIndex, 1);
                    y = features(featureIndex, 2);
                    location_i = floor(y / obj.quadHeight) + 1;
                    location_j = floor(y / obj.quadWidth) + 1;
                    %B = squeeze(obj.data(frameIndex).C(location_i, location_j, :, :))  * squeeze(obj.data(frameIndex).P{obj.maxIte}(location_i, location_j, :, :));
                    B = squeeze(obj.data(frameIndex).P{obj.maxIte}(location_i, location_j, :, :)) * squeeze(obj.data(frameIndex).C_inv(location_i, location_j, :, :));
                    [x_2, y_2 ] = obj.transform(x, y, B);
                    features_2(featureIndex, 1) = x_2;
                    features_2(featureIndex, 2) = y_2;
                end
                % asap solve and warp
                asap = AsSimilarAsPossibleWarping(obj.videoHeight,obj.videoWidth,obj.quadWidth * 2,obj.quadHeight * 2, 0.001);
                asap.ADAPTIVE_WEIGHT = 0;
                asap.SetControlPts(features, features_2);%set matched features
                asap.Solve();       
%                 marked = obj.drawFeatures(I2, features);
                imwarp = asap.Warp(marked, obj.gap);
                % save image file
                
                imwrite(uint8(imwarp), [obj.outputPath '/' int2str(frameIndex) '.png']);
            end 
        end
        

        function saveSeq(obj)
            for index = 1 : obj.num_frame
               % imwrite(obj.imwarp(index).cdata, ['result/' int2str(index) '.png']); 
            end 
        end
    end
    
    methods (Access = private)
        function [homo, err] = getHomo(obj, I1, I2, asaplamda, ADPTIVE_WEIGHT, frameIndex, maskIn)
            
            % I1 = frames(frameIndex).cdata;
            % I2 = frames(frameIndex+1).cdata;
            [I1_features,I2_features, score, I1_all, I2_all]=SIFT2(I1,I2,maskIn);
            fprintf('Frame: %d, Score: %f\n', frameIndex, score);
%             if score < 0.5
%                 asaplamda = 5;
%                 [I1_features,I2_features, score] = SIFT2(I1,I2,ones(100, 1));
%                 fprintf('Frame: %d, Score: %f\n', frameIndex, score);
%             end
            if score < 0.1
                error('?');
            end
            mask1 = obj.getValidMask(I1_features, I1_all);
            mask2 = obj.getValidMask(I2_features, I2_all);
            
            % [I1_features, I2_features] = saliencyFilter(I1_features, I2_features, frameIndex);
            
            if sum(mask1) < 4
                error('not enough matched features');
            end
            
            %obj.features(frameIndex) = I1_features;
            
            asap = AsSimilarAsPossibleWarping(obj.videoHeight,obj.videoWidth,obj.quadWidth,obj.quadHeight,asaplamda);
            asap.ADAPTIVE_WEIGHT = 0;
            asap.SetControlPts(I1_features,I2_features);%set matched features
            asap.Solve();            %solve Ax=b for as similar as possible
%             qW = floor(obj.videoWidth / 10);
%             qH = floor(obj.videoHeight / 10);
%             if sum(maskIn) ~= 101
%                 markedI1 = obj.drawFeatures(I1, I1_all, 'red');
%                 markedI2 = obj.drawFeatures(I2, I2_all, 'red');
%                 markedI1 = obj.drawFeatures(markedI1, I1_features, 'green');
%                 markedI2 = obj.drawFeatures(markedI2, I2_features, 'green');
%                 WarpIm = asap.Warp(markedI1, 100);
%                 for i = 1:10
%                     for j = 1:10
%                         if (mask1((i-1)*10 + j) == 1)
%                              markedI1 = insertShape(markedI1, 'Rectangle', [(i-1) * qW + 2 (j-1) * qH + 2 qW - 4 qH - 4]);
%                         end
%                         if (mask2((i-1)*10 + j) == 1)
%                              markedI2 = insertShape(markedI2, 'Rectangle', [(i-1) * qW + 2 (j-1) * qH + 2 qW - 4 qH - 4]);
%                         end
%                     end
%                 end
%                 
%                 
%                 imwrite(WarpIm, ['../foo_reference/' int2str(frameIndex) '_warp.png'])
%                 imwrite(markedI1, ['../foo_reference/' int2str(frameIndex) '_1.png'])
%                 imwrite(markedI2, ['../foo_reference/' int2str(frameIndex) '_2.png'])
%             end
            homo = asap.CalcHomos();            
%             err = asap.CalcError();
            err = 0;
            obj.validMask(1:100, frameIndex) = mask1;
            obj.validMask(101:200, frameIndex) = mask2;
        end
		
        function imwarp = warp1(obj, frameIndex, location_i, location_j, source, imwarp)
%             if location_i == 16 && location_j == 7
%                 fooooo = 1; 
%             end
            % disp([location_i location_j])
            %B = squeeze(obj.data(frameIndex).C(location_i, location_j, :, :)) * squeeze(obj.data(frameIndex).P{obj.maxIte}(location_i, location_j, :, :));
            B = squeeze(obj.data(frameIndex).P{obj.maxIte}(location_i, location_j, :, :)) * squeeze(obj.data(frameIndex).C_inv(location_i, location_j, :, :));
            %B = B ./B(3,3);
            minx = (location_i - 1)*obj.quadHeight+1;
            maxx = minx + obj.quadHeight;
            miny = (location_j - 1)*obj.quadWidth+1;
            maxy = miny + obj.quadWidth;
            [x00, y00 ] = obj.transform(miny, minx, B);
            [x01, y01 ] = obj.transform(miny, maxx, B);
            [x10, y10 ] = obj.transform(maxy, minx, B);
            [x11, y11 ] = obj.transform(maxy, maxx, B);
            minx = min(x00, x01);minx = min(minx, x10);minx = min(minx, x11);
            miny = min(y00, y10);miny = min(miny, y01);miny = min(miny, y11);
            maxx = max(x10, x11);maxx = max(maxx, x01);maxx = max(maxx, x00);
            maxy = max(y01, y11);maxy = max(maxy, y10);maxy = max(maxy, y00);
            % imwarp = myWarp(minx, maxx, miny, maxy, double(source), obj.imwarp(frameIndex).cdata, inv(B), obj.gap);
            imwarp = myWarp(minx, maxx, miny, maxy, double(source), imwarp, inv(B), obj.gap);
        end
        function [x,y] = transform(~, xx, yy, B)
            res = B * [xx;yy;1];
            x = res(1)/res(3);
            y = res(2)/res(3);
        end
        
        function H = normH(~, H)
            H = H./H(3, 3); 
        end
    end
    methods (Access = private) %functions of warping
        function imwarp = quadWarp(obj,im,q1,q2, imwarp)
            
            minx = q2.getMinX();
            maxx = q2.getMaxX();
            miny = q2.getMinY();
            maxy = q2.getMaxY();
            
            source = zeros(4,2);
            target = zeros(4,2);
            
            source(1,1) = q2.V00.x;source(1,2) = q2.V00.y;
            source(2,1) = q2.V01.x;source(2,2) = q2.V01.y;
            source(3,1) = q2.V10.x;source(3,2) = q2.V10.y;
            source(4,1) = q2.V11.x;source(4,2) = q2.V11.y;
            
            target(1,1) = q1.V00.x;target(1,2) = q1.V00.y;
            target(2,1) = q1.V01.x;target(2,2) = q1.V01.y;
            target(3,1) = q1.V10.x;target(3,2) = q1.V10.y;
            target(4,1) = q1.V11.x;target(4,2) = q1.V11.y;
            
            H = homography_4pts(source',target');
            H = H./H(3,3);
            
            %qd = Quad(q2.V00,q2.V01,q2.V10,q2.V11);
            imwarp = myWarp(minx,maxx,miny,maxy,double(im),imwarp,H,obj.gap);
            imwarp = uint8(imwarp);
            
        end
        
        function imout = drawFeatures(~, im, features, color)
            imout = insertMarker(im, features, 'color', color);
        end
        
        function mask = getValidMask(obj, features, all)
            [m, n] = size(features);
            if m == 2 && n > 10
                features = features';
            end
            qH = obj.videoHeight / 10;
            qW = obj.videoWidth / 10;
    
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
    end
end


classdef VideoStitch < handle
    %VIDEOSTITCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % inputs
        seqA;
        seqB;
        
        % Dimensions
        nFrames;
        
        videoHeight;
        videoWidth;
        
        meshSize;
        
        quadHeight;
        quadWidth;
        
        % correspondence
        maxppf;
        ppf; % points per frame
        CP; % Control Points
        validCP; 
        nCP;
        
        % optimization parameters
        span;
        smoothness;
        cropping;
        stitchness;
        stableW;
        gap;
        
        % Optimization data
        Pa;
        Pb;
        Pa_inv;
        Pb_inv;
        Ca;
        Cb;
        Ca_inv;
        Cb_inv;
        H;
        
        Offset;
        
        % Tuples for building sparse matrix
        % 3 types of quadratic term
        % Pt = Ct
        % 
        % 
        % Pt = Pr
        % Pa*CPa = H*Pb*CPb
        % <x, y, z> = H*P*<x, y, z>
        
        % options
        useImage;
        
        % Others 
        
        % a2b2 for Pt - Ct
        % C^(-1) * corners
        CaCorner; %zeros(obj.nFrames, obj.meshSize, obj.meshSize, 4, 2);
        CbCorner;
        PaCorner;
        PbCorner;
        % C^(-1) * CP
        CaCP;
        CbCP;
        
        w_a;
        w_b;
       
    end
    
    methods
        function obj = VideoStitch(seqA, seqB, pathA, pathB, ControlPoints, ppf, UseImageSeq, smoothness, cropping, stitchness)
            obj.seqA = seqA;
            obj.seqB = seqB;
            obj.Ca = pathA;
            obj.Cb = pathB;
            obj.Pa = obj.Ca;
            obj.Pb = obj.Cb;
            obj.CP = ControlPoints;
            obj.smoothness = smoothness;
            obj.cropping = cropping;
            obj.stitchness = stitchness;
            obj.useImage = UseImageSeq;
            obj.stableW = 40;
            obj.span = 30;
            obj.H = eye(3);
            if UseImageSeq
                fileList = dir(seqA);
                fileList = fileList(3:length(fileList));
                obj.nFrames = length(fileList);
                if obj.nFrames < 2
                    error('Wrong inputs directory') ;
                end
                fileList = dir(seqB);
                fileList = fileList(3:length(fileList));
                if obj.nFrames ~= length(fileList)
                    error('Input length  doesn"t match') ;
                end
                frame = imread([seqB fileList(1).name]);
            else
                video = VideoReader(seqA);
                if hasFrame(video)
                    frame = readFrame(video);
                end
            end
            [obj.videoHeight, obj.videoWidth, ~] = size(frame);
            [~, obj.meshSize, ~, ~, ~] = size(pathA);
            [~, obj.meshSize, ~, ~, ~] = size(pathB);
            [obj.nFrames, obj.maxppf, ~] = size(obj.CP);
            obj.ppf = ppf;
            obj.validCP = zeros(obj.nFrames, obj.maxppf);
            obj.nCP = sum(ppf);
            obj.quadHeight = obj.videoHeight / obj.meshSize;
            obj.quadWidth = obj.videoWidth / obj.meshSize;
            
            
            
            %compute C_inv
            obj.Ca_inv = zeros(size(obj.Ca));
            obj.Cb_inv = zeros(size(obj.Cb));
            
            
            obj.CaCP = zeros(obj.nFrames, obj.maxppf, 3);
            obj.CbCP = obj.CaCP;
            
            obj.w_a = zeros(obj.nFrames, obj.nFrames, obj.meshSize, obj.meshSize);
            obj.w_b = zeros(obj.nFrames, obj.nFrames, obj.meshSize, obj.meshSize);
        end
        
        function [Pa_inv, Pb_inv] = inverse(obj, Pa, Pb)
            Pa_inv = zeros(size(Pa));
            Pb_inv = zeros(size(Pb));
            for frameIndex = 1:obj.nFrames
                for i = 1:obj.meshSize
                    for j = 1:obj.meshSize
                        Pa_inv(frameIndex, i, j, :, :) = squeeze(Pa(frameIndex, i, j, :, :))^(-1);                         
                        Pb_inv(frameIndex, i, j, :, :) = squeeze(Pb(frameIndex, i, j, :, :))^(-1);   
                    end
                end
            end
        end
        
        function [PaCorner, PbCorner] = PxCorner(obj, Pa_inv, Pb_inv)
            PaCorner = ones(obj.nFrames, obj.meshSize, obj.meshSize, 4, 3);
            PbCorner = ones(obj.nFrames, obj.meshSize, obj.meshSize, 4, 3);
            for frameIndex = 1:obj.nFrames
                for i = 1:obj.meshSize
                    for j = 1:obj.meshSize
                        [PaCorner(frameIndex, i, j, 1, 1), PaCorner(frameIndex, i, j, 1, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i - 1) * obj.quadHeight + 1], squeeze(Pa_inv(frameIndex, i, j, :, :))); 
                        [PaCorner(frameIndex, i, j, 2, 1), PaCorner(frameIndex, i, j, 2, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i - 1) * obj.quadHeight + 1], squeeze(Pa_inv(frameIndex, i, j, :, :)));
                        [PaCorner(frameIndex, i, j, 3, 1), PaCorner(frameIndex, i, j, 3, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i) * obj.quadHeight], squeeze(Pa_inv(frameIndex, i, j, :, :)));
                        [PaCorner(frameIndex, i, j, 4, 1), PaCorner(frameIndex, i, j, 4, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i) * obj.quadHeight], squeeze(Pa_inv(frameIndex, i, j, :, :)));
                        [PbCorner(frameIndex, i, j, 1, 1), PbCorner(frameIndex, i, j, 1, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i - 1) * obj.quadHeight + 1], squeeze(Pb_inv(frameIndex, i, j, :, :))); 
                        [PbCorner(frameIndex, i, j, 2, 1), PbCorner(frameIndex, i, j, 2, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i - 1) * obj.quadHeight + 1], squeeze(Pb_inv(frameIndex, i, j, :, :)));
                        [PbCorner(frameIndex, i, j, 3, 1), PbCorner(frameIndex, i, j, 3, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i) * obj.quadHeight], squeeze(Pb_inv(frameIndex, i, j, :, :)));
                        [PbCorner(frameIndex, i, j, 4, 1), PbCorner(frameIndex, i, j, 4, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i) * obj.quadHeight], squeeze(Pb_inv(frameIndex, i, j, :, :)));
                    end
                end
            end
        end
        
        function init(obj)
            obj.Pa = obj.Ca;
            obj.Pb = obj.Cb;
            
            [obj.Ca_inv, obj.Cb_inv] = obj.inverse(obj.Ca, obj.Cb);
            [obj.Pa_inv, obj.Pb_inv] = obj.inverse(obj.Pa, obj.Pb);    
            
            [obj.CaCorner, obj.CbCorner] = obj.PxCorner(obj.Ca_inv, obj.Cb_inv);
            [obj.PaCorner, obj.PbCorner] = obj.PxCorner(obj.Pa_inv, obj.Pb_inv);
            
            for frameIndex = 1:obj.nFrames
                for k = 1:obj.ppf(frameIndex)
                    obj.CaCP(frameIndex, k, :) = obj.timesCa([obj.CP(frameIndex, k, 1) obj.CP(frameIndex, k, 2) frameIndex]);
                    obj.CbCP(frameIndex, k, :) = obj.timesCb([obj.CP(frameIndex, k, 3) obj.CP(frameIndex, k, 4) frameIndex]);
                    obj.validCP(frameIndex, k) = 1;
                end
            end
            obj.calcOmega();
            obj.H = eye(3);
            obj.nCP = sum(obj.ppf);
            fprintf('Number of Valid Control Points :%5d\n', obj.nCP);
        end
        
        function updateAll(obj)
            
            obj.updateH_simple();
            obj.updateCP();
            [obj.Pa_inv, obj.Pb_inv] = obj.inverse(obj.Pa, obj.Pb);
            [obj.PaCorner, obj.PbCorner] = obj.PxCorner(obj.Pa_inv, obj.Pb_inv);            
        end
        
        function updateOffset(obj)
            B = zeros(obj.nFrames, 3, 3);
            tPa = zeros(obj.nFrames, obj.meshSize, obj.meshSize, 3, 3);
            tPb = tPa;
            oPa = obj.Pa;
            oPb = obj.Pb;
            oCa = obj.Ca;
            ms = obj.meshSize;
            for frameIndex = 1:obj.nFrames
                for row = 1:ms
                    for col = 1:ms 
                        B1 = squeeze(oPa(frameIndex, row, col, :, :)) / squeeze(oCa(frameIndex, row, col, :, :));
                        B(frameIndex, :, :) = squeeze(B(frameIndex, :, :)) + B1 ./ B1(3,3); 
                    end
                end
            end
            nf = obj.nFrames;
%             for frameIndex = 1:nf
%                 head = max(frameIndex - window, 1);
%                 tail = min(frameIndex + window, nf);
%                 offset = sum(B(head:tail, :, :)) / (tail - head + 1) / (ms * ms);
%                 for row = 1:ms
%                     for col = 1:ms 
%                         tPa(frameIndex, row, col, :, :) = squeeze(offset) \ squeeze(oPa(frameIndex, row, col, :, :));
%                         tPb(frameIndex, row, col, :, :) = squeeze(offset) \ squeeze(oPb(frameIndex, row, col, :, :));
%                     end
%                 end
%                 
%             end
            offset = squeeze(sum(B(:, :, :)) / nf / (ms * ms));
            t = 2 * obj.span + 1;
            [ax00, ay00] = obj.transform([1, 1], offset \ squeeze(obj.Pa(t, obj.meshSize, 1, :, :)) / squeeze(obj.Ca(t, obj.meshSize, 1, :, :)));
            [ax01, ay01] = obj.transform([1, obj.videoHeight], offset \ squeeze(obj.Pa(t, 1, 1, :, :)) / squeeze(obj.Ca(t, 1, 1, :, :)));
            [ax10, ay10] = obj.transform([obj.videoWidth, 1], offset \ squeeze(obj.Pa(t, obj.meshSize, obj.meshSize, :, :)) / squeeze(obj.Ca(t, obj.meshSize, obj.meshSize, :, :)));
            [ax11, ay11] = obj.transform([obj.videoWidth, obj.videoHeight], offset \ squeeze(obj.Pa(t, obj.meshSize, 1, :, :)) / squeeze(obj.Ca(t, obj.meshSize, 1, :, :)));
            
            [bx00, by00] = obj.transform([1, 1], offset \ obj.H * squeeze(obj.Pb(t, obj.meshSize, 1, :, :)) / squeeze(obj.Cb(t, obj.meshSize, 1, :, :)));
            [bx01, by01] = obj.transform([1, obj.videoHeight], offset \ obj.H * squeeze(obj.Pb(t, 1, 1, :, :)) / squeeze(obj.Cb(t, 1, 1, :, :)));
            [bx10, by10] = obj.transform([obj.videoWidth, 1], offset \ obj.H * squeeze(obj.Pb(t, obj.meshSize, obj.meshSize, :, :)) / squeeze(obj.Cb(t, obj.meshSize, obj.meshSize, :, :)));
            [bx11, by11] = obj.transform([obj.videoWidth, obj.videoHeight], offset \ obj.H * squeeze(obj.Pb(t, obj.meshSize, 1, :, :)) / squeeze(obj.Cb(t, obj.meshSize, 1, :, :)));
            
            a = [ax00 ay00; ax01 ay01; ax10 ay10; ax11 ay11];
            b = [bx00 by00; bx01 by01; bx10 by10; bx11 by11];
            c = (a + b) * 0.5;
            
            ha = homography_4pts(c', a');
            hb = homography_4pts(b', c');            
            
%             for frameIndex = 1:nf
%                 for row = 1:ms
%                     for col = 1:ms 
%                         tPa(frameIndex, row, col, :, :) = ha * squeeze(tPa(frameIndex, row, col, :, :));
%                         tPb(frameIndex, row, col, :, :) = ha * squeeze(tPb(frameIndex, row, col, :, :));
%                     end
%                 end
%                 
%             end
%             
%             obj.Pa = tPa;
%             obj.Pb = tPb;
            obj.Offset = ha / offset;
        end
        
        function value = getCoherenceTerm(obj, ab, P, row, col, frameIndex)
            value = 0;
            if ab == 'a'
                C = obj.Ca;
            else
                C = obj.Cb;                
            end
            if obj.meshSize == 1
                return 
            end
            
            % i=1, j=1
            if row == 1 && col == 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=1, col=2:15
            if row == 1 && col < obj.meshSize && col > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=1, col=16
            if row == 1 && col == obj.meshSize
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=2:15, col=1
            if col == 1 && row < obj.meshSize && row > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=2:15, col=2:15
            if row > 1 && row < obj.meshSize && col > 1 && col < obj.meshSize
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :)); 
            end
            % row=2:15, col=16
            if col == obj.meshSize && row < obj.meshSize && row > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=16, col=1
            if row == obj.meshSize && col == 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=16, col=2:15
            if row == obj.meshSize && col < obj.meshSize && col > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=16, col=16
            if row == obj.meshSize && col == obj.meshSize
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            
        end
        
        function optPath(obj, maxIte)
            % first round
            
            obj.optPa(1);
            obj.optPb(1);
            obj.updateAll();
            ite = 0;
            while ite < maxIte 
                ite = ite + 1;
                disp(['round#' int2str(ite)]);
                done1 = obj.optPa(ite - maxIte);
                done2 = obj.optPb(ite - maxIte);
                obj.updateAll();
                if done1 * done2 == 1
                    break; 
                end
            end
        end
        
        
     
        function done = optPa(obj, is1st)
            % get J and r
            % iteratively get updated Pa and Pb
            if is1st == 0
                maxIte = 30;
            else
                maxIte = 20;
            end
            
            boost = 1;
            done = 0;
            for ite = 1: maxIte
                fprintf('.');
                J = obj.computeJa(is1st);
                r = obj.computeRa(is1st);
                b = - J' * r;
                A = J' * J;
                dx = A \ b;
                obj.update('a', dx * boost);
                rr = r' * r
                if ite == 1
                    minimum = rr;
                else
                    if abs(rr - minimum) < 1e4
                        disp('Pa opted');
                        if ite == 2
                            done = 1;
                        end
                        break; 
                    end
                    if minimum <= rr
                        boost = boost * 0.5;
                        obj.update('a', - dx * boost);
                    else
                        minimum = rr;
                    end                    
                end      
            end             
            if ite == maxIte
                disp('reached max iteration'); 
            end
        end
        
        function done = optPb(obj, is1st)
            % get J and r
            % iteratively get updated Pa and Pb
            if is1st == 0
                maxIte = 30;
            else
                maxIte = 20;
            end
            boost = 1;
            done = 0;
            for ite = 1: maxIte
                fprintf('.');
                J = obj.computeJb(is1st);
                r = obj.computeRb(is1st);
                b = - J' * r;
                A = J' * J;
                dx = A \ b;
                obj.update('b', dx * boost);
                rr = r' * r
                if ite == 1
                    minimum = rr;
                else
                    if abs(rr - minimum) < 1e4
                        disp('Pb opted');
                        if ite == 2
                            done = 1;
                        end
                        break; 
                    end
                    if minimum <= rr
                        boost = boost * 0.5;
                        obj.update('b', - dx * boost);
                    else
                        minimum = rr;
                    end                    
                end      
            end             
            if ite == maxIte
                disp('reached max iteration'); 
            end
        end
        
        function update(obj, path, dx)
            if path == 'a'
                P = obj.Pa; 
            else
                P = obj.Pb;
            end
            for frameIndex = 1 : obj.nFrames
                for row = 1:obj.meshSize
                    for col = 1:obj.meshSize 
                        unknownIndex = ((frameIndex - 1) * obj.meshSize * obj.meshSize + (row - 1) * obj.meshSize + col) * 8;
                        dx33 = reshape([dx(unknownIndex - 7 : unknownIndex) ; 0], [3,3]);
                        dx33 = dx33';
                        P(frameIndex, row, col, :, :) = squeeze(P(frameIndex, row, col, :, :)) + dx33; 
                    end
                end
            end
            if path == 'a'
                obj.Pa = P; 
            else
                obj.Pb = P;
            end
            [obj.Pa_inv, obj.Pb_inv] = obj.inverse(obj.Pa, obj.Pb);
            [obj.PaCorner, obj.PbCorner] = obj.PxCorner(obj.Pa_inv, obj.Pb_inv);
        end
        
        function updateH(obj)
            % compute H from Pa Ca Pb Cb and CP 
            PA = zeros(3, obj.nCP);
            PB = zeros(3, obj.nCP);
            CPcount = 0;
            for frameIndex = obj.span + 1 : obj.nFrames - obj.span
                cpa = squeeze(obj.CaCP(frameIndex, obj.validCP(frameIndex, :) == 1, :));
                cpb = squeeze(obj.CbCP(frameIndex, obj.validCP(frameIndex, :) == 1, :));
                % TODO
                PA(:, CPcount + 1 : CPcount + obj.ppf(frameIndex)) = squeeze(obj.Pa(frameIndex, 1, 1, :, :)) * cpa';
                PB(:, CPcount + 1 : CPcount + obj.ppf(frameIndex)) = squeeze(obj.Pb(frameIndex, 1, 1, :, :)) * cpb';
                CPcount = CPcount + obj.ppf(frameIndex);
            end
            [obj.H, ~] = EstimateHomographyByRANSAC(PA, PB, 0.001);
            obj.H = obj.H / obj.H(3, 3);
            obj.H = obj.H ^ (-1);
        end
        
        function updateH_simple(obj)
            % compute H from Pa Ca Pb Cb and CP 
            PA = zeros(2, obj.ppf(obj.span+1) + obj.ppf(obj.span+2));
            PB = zeros(2, obj.ppf(obj.span+1) + obj.ppf(obj.span+2));
            CPcount = 0;
            for frameIndex = floor(obj.nFrames * 0.5) : floor(obj.nFrames * 0.5) + 2
                for k = 1:obj.maxppf
                    if obj.validCP(frameIndex, k) == 0
                        break                        
                    end
                    if obj.validCP(frameIndex, k) == 1
                        %xa, ya, xb, yb;
                        pb = obj.CP(frameIndex, k, 3:4);
                        colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                        rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                        Bb = squeeze(obj.Pb(frameIndex, rowb, colb, :, :)) / squeeze(obj.Cb(frameIndex, rowb, colb, :, :));
                        pa = obj.CP(frameIndex, k, 1:2);
                        cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                        rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1; 
                        Ba = squeeze(obj.Pa(frameIndex, rowa, cola, :, :)) / squeeze(obj.Ca(frameIndex, rowa, cola, :, :));
                        CPcount = CPcount + 1;
                        [PA(1, CPcount), PA(2, CPcount)] = obj.transform(pa, Ba);
                        [PB(1, CPcount), PB(2, CPcount)] = obj.transform(pb, Bb);
                    end                        
                end
            end
            [obj.H, ~] = EstimateHomographyByRANSAC(PA, PB, 0.001);
            obj.H = obj.H / obj.H(3, 3);
            obj.H = obj.H ^ (-1);
        end
        
        function updateCP(obj)
            % set validCp and ppf according to H
%             threshold = 0.06 * obj.videoWidth / ite;
            obj.ppf = zeros(obj.nFrames, 1);
            for frameIndex = 1:obj.nFrames
                d = zeros(obj.maxppf, 1);
                for k = 1:obj.maxppf
                    if obj.validCP(frameIndex, k) == 0
                        break;
                    end
                    pb = obj.CP(frameIndex, k, 3:4);
                    colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                    rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                    pa = obj.CP(frameIndex, k, 1:2);
                    cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                    rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                    [a2, b2] = obj.transform(pb, obj.H * squeeze(obj.Pb(frameIndex, rowb, colb, :, :)) / squeeze(obj.Cb(frameIndex, rowb, colb, :, :)));
                    [a1, b1] = obj.transform(pa, squeeze(obj.Pa(frameIndex, rowa, cola, :, :)) / squeeze(obj.Ca(frameIndex, rowa, cola, :, :)));
                    d(k) = sqrt((a1 - a2)^2 + (b1 - b2)^2);                            
                end
                dsort = sort(d(d ~= 0));
                threshold = mean(dsort(floor(length(dsort) / 4):floor(length(dsort) * 3/ 4)));
                threshold = max(threshold, 0.02 * obj.videoWidth);
                threshold = min(threshold, 0.08 * obj.videoWidth);
                obj.validCP(frameIndex, d <= threshold) = 1;
                obj.validCP(frameIndex, d == 0) = 0;
                obj.validCP(frameIndex, d > threshold) = 0;
                obj.ppf(frameIndex) = sum(obj.validCP(frameIndex, :));
                obj.validCP(frameIndex, d > threshold) = -1;
                if obj.ppf(frameIndex) < 20
                    obj.validCP(frameIndex, obj.validCP(frameIndex, :) ~= 0) = -1; 
                    obj.ppf(frameIndex) = 0;
                end
            end
            obj.nCP = sum(obj.ppf);
            fprintf('\nNumber of Valid Control Points :%5d\n', obj.nCP);
        end
        
        function J = computeJa(obj, is1st)
            % get the last updated J matrix 
            rowCount = 0;
            itemCount = 0;
            if is1st == 1
                nItems = obj.nFrames * (1 + 2 * obj.span) * obj.meshSize * obj.meshSize * 4 * 2 * 8;
            else
                nItems = (obj.nCP + obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4) * 2 * 8;
            end
            Jbuilder = zeros(nItems, 3);
            for frameIndex = 1 : obj.nFrames
                for row = 1: obj.meshSize
                    for col = 1 : obj.meshSize
                        unknownIndex = ((frameIndex - 1) * obj.meshSize * obj.meshSize + (row - 1) * obj.meshSize + col) * 8;
                        for r = frameIndex - obj.span : frameIndex + obj.span
                            if r == frameIndex
                                weight = sqrt(obj.cropping);
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 1, 1), obj.CaCorner(frameIndex, row, col, 1, 2));                              
                                
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 2, 1), obj.CaCorner(frameIndex, row, col, 2, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 3, 1), obj.CaCorner(frameIndex, row, col, 3, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 4, 1), obj.CaCorner(frameIndex, row, col, 4, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                            elseif r > 0 && r <= obj.nFrames
                                weight = sqrt(obj.smoothness * obj.w_a(frameIndex, r, row, col));
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 1, 1), obj.PaCorner(r, row, col, 1, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 2, 1), obj.PaCorner(r, row, col, 2, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 3, 1), obj.PaCorner(r, row, col, 3, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 4, 1), obj.PaCorner(r, row, col, 4, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                            end
                        end
                    end
                end
                if is1st ~= 1
                    for k = 1:obj.maxppf
                        if obj.validCP(frameIndex, k) == 1
                            weight = 1;
                            pa = obj.CP(frameIndex, k, 1:2);
                            cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                            rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                            unknownIndex = ((frameIndex - 1) * obj.meshSize * obj.meshSize + (rowa - 1) * obj.meshSize + cola) * 8;

                            values = obj.getSubJType1(squeeze(obj.Pa(frameIndex, rowa, cola, :, :)), ...
                                obj.CaCP(frameIndex, k, 1), obj.CaCP(frameIndex, k, 2)) * weight;    
                            rowCount = rowCount + 1;
                            Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                            Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                            Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :);
                            itemCount = itemCount + 8;
                            rowCount = rowCount + 1;
                            Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                            Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                            Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :);
                            itemCount = itemCount + 8;
                        end                    
                    end
                end
            end
            Jbuilder(itemCount+1:length(Jbuilder), :) = [];
            if itemCount ~= size(Jbuilder, 1)
                error('?'); 
            end
            J = sparse(Jbuilder(:, 1), Jbuilder(:, 2), Jbuilder(:, 3), rowCount, obj.nFrames * 8 * obj.meshSize * obj.meshSize, itemCount);
        end
        
        function J = computeJb(obj, is1st)
            % get the last updated J matrix 
            rowCount = 0;
            itemCount = 0;
            if is1st == 1
                nItems = obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4 * 2 * 8;
            else
                nItems = (obj.nCP + obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4) * 2 * 8;
            end
            Jbuilder = zeros(nItems, 3);
            for frameIndex = 1 : obj.nFrames 
                for row = 1: obj.meshSize
                    for col = 1 : obj.meshSize
                        unknownIndex = ((frameIndex - 1) * obj.meshSize * obj.meshSize + (row - 1) * obj.meshSize + col) * 8;
                        for r = frameIndex - obj.span : frameIndex + obj.span
                            if r == frameIndex
                                weight = sqrt(obj.cropping);
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 1, 1), obj.CbCorner(frameIndex, row, col, 1, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 2, 1), obj.CbCorner(frameIndex, row, col, 2, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 3, 1), obj.CbCorner(frameIndex, row, col, 3, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 4, 1), obj.CbCorner(frameIndex, row, col, 4, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                            elseif r > 0 && r <= obj.nFrames
                                weight = sqrt(obj.smoothness * obj.w_b(frameIndex, r, row, col));
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 1, 1), obj.PbCorner(r, row, col, 1, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 2, 1), obj.PbCorner(r, row, col, 2, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 3, 1), obj.PbCorner(r, row, col, 3, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                                
                                values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 4, 1), obj.PbCorner(r, row, col, 4, 2));
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :) * weight;
                                itemCount = itemCount + 8;
                                rowCount  = rowCount + 1;
                                Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                                Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                                Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :) * weight;
                                itemCount = itemCount + 8;
                            end
                        end
                    end
                end
                if is1st ~= 1
                    for k = 1:obj.maxppf
                        if obj.validCP(frameIndex, k) == 1
                            weight = 1;
                            pb = obj.CP(frameIndex, k, 3:4);
                            colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                            rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                            unknownIndex = ((frameIndex - 1) * obj.meshSize * obj.meshSize + (rowb - 1) * obj.meshSize + colb) * 8;

                            values = obj.getSubJType1(squeeze(obj.Pb(frameIndex, rowb, colb, :, :)), ...
                                obj.CbCP(frameIndex, k, 1), obj.CbCP(frameIndex, k, 2)) * weight;    

                            rowCount = rowCount + 1;
                            Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                            Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                            Jbuilder(itemCount + 1:itemCount + 8, 3) = values(1, :);
                            itemCount = itemCount + 8;
                            rowCount = rowCount + 1;
                            Jbuilder(itemCount + 1:itemCount + 8, 1) = rowCount;
                            Jbuilder(itemCount + 1:itemCount + 8, 2) = unknownIndex - 7 : unknownIndex;
                            Jbuilder(itemCount + 1:itemCount + 8, 3) = values(2, :);
                            itemCount = itemCount + 8;
                        end
                    end
                end
            end
            Jbuilder(itemCount+1:length(Jbuilder), :) = [];
            if itemCount ~= size(Jbuilder, 1)
                error('?'); 
            end
            J = sparse(Jbuilder(:, 1), Jbuilder(:, 2), Jbuilder(:, 3), rowCount, obj.nFrames * 8 * obj.meshSize * obj.meshSize, nItems);
        end
                
        function subJ = getSubJType1(~, H, a, b)
            % for Pt - Ct and Pt - H
            H = H./ H(3,3);
            x = reshape(H', [9 1]);
            % return 2*8
            deno = a * x(7) + b * x(8) + 1;
            
            nume = a * x(1) + b * x(2) + x(3);
            subJ = zeros(2, 8);
            subJ(1, 1) = a / deno; subJ(1, 2) = b / deno; subJ(1, 3) = 1 / deno;
            subJ(1, 7) = - a * nume / deno^2; subJ(1, 8) = - b * nume / deno^2;
            
            nume = a * x(4) + b * x(5) + x(6);
            subJ(2, 4) = a / deno; subJ(2, 5) = b / deno; subJ(2, 6) = 1 / deno;
            subJ(2, 7) = - a * nume / deno^2; subJ(2, 8) = - b * nume / deno^2;
        end
        
        function subJ = getSubJType3(~, Ht, Hr, a, b)
            Ht = Ht./Ht(3,3);
            xt = reshape(Ht', [9 1]);
            Hr = Hr./Hr(3,3);
            xr = reshape(Hr', [9 1]);
            
            
            
            % return 2 * 2 * 8 
            subJ = zeros(2, 2, 8);
            denot = a * xt(7) + b * xt(8) + 1;
            denor = a * xr(7) + b * xr(8) + 1;
            
            numet = a * xt(1) + b * xt(2) + xt(3);
            numer = a * xr(1) + b * xr(2) + xr(3);
            subJ(1, 1, 1) = a / denot; subJ(1, 1, 2) = b / denot; subJ(1, 1, 3) = 1 / denot;
            subJ(1, 2, 1) = - a / denor; subJ(1, 2, 2) = - b / denor; subJ(1, 2, 3) = - 1 / denor;
            subJ(1, 1, 7) = - a * numet / denot^2; subJ(1, 1, 8) = - b * numet / denot^2;
            subJ(1, 2, 7) = a * numer / denor^2; subJ(1, 2, 8) = b * numer / denor^2;
            
            numet = a * xt(4) + b * xt(5) + xt(6);
            subJ(2, 1, 4) = a / denot; subJ(2, 1, 5) = b / denot; subJ(2, 1, 6) = 1 / denot;
            subJ(2, 1, 7) = - a * numet / denot^2; subJ(2, 1, 8) = - b * numet / denot^2;
            
            numer = a * xr(4) + b * xr(5) + xr(6);
            subJ(2, 2, 4) = - a / denor; subJ(2, 2, 5) = - b / denor; subJ(2, 2, 6) = -1 / denor;
            subJ(2, 2, 7) = a * numer / denor^2; subJ(2, 2, 8) = b * numer / denor^2;
            
        end
        
        function subJ = getSubJType2(~, Ht, Hr, at, bt, ar, br)
            Ht = Ht./Ht(3,3);
            xt = reshape(Ht', [9 1]);
            Hr = Hr./Hr(3,3);
            xr = reshape(Hr', [9 1]);
            
            
            
            % return 2 * 2 * 8 
            subJ = zeros(2, 2, 8);
            denot =at * xt(7) + bt * xt(8) + 1;
            denor =ar * xr(7) + br * xr(8) + 1;
            
            numet =at * xt(1) + bt * xt(2) + xt(3);
            numer =ar * xr(1) + br * xr(2) + xr(3);
            subJ(1, 1, 1) = at / denot; subJ(1, 1, 2) = bt / denot; subJ(1, 1, 3) = 1 / denot;
            subJ(1, 2, 1) = - ar / denor; subJ(1, 2, 2) = - br / denor; subJ(1, 2, 3) = - 1 / denor;
            subJ(1, 1, 7) = - at * numet / denot^2; subJ(1, 1, 8) = - bt * numet / denot^2;
            subJ(1, 2, 7) = ar * numer / denor^2; subJ(1, 2, 8) = br * numer / denor^2;
            
            numet = at * xt(4) + bt * xt(5) + xt(6);
            subJ(2, 1, 4) = at / denot; subJ(2, 1, 5) = bt / denot; subJ(2, 1, 6) = 1 / denot;
            subJ(2, 1, 7) = - at * numet / denot^2; subJ(2, 1, 8) = - bt * numet / denot^2;
            
            numer = ar * xr(4) + br * xr(5) + xr(6);
            subJ(2, 2, 4) = - ar / denor; subJ(2, 2, 5) = - br / denor; subJ(2, 2, 6) = -1 / denor;
            subJ(2, 2, 7) = ar * numer / denor^2; subJ(2, 2, 8) = br * numer / denor^2;
            
        end
        
        function Ra = computeRa(obj, is1st)
            % return a vector with the length of number of quadratic terms 
            if is1st
                Ra = zeros(obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4 * 2, 1);
            else
                Ra = zeros((obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4 + obj.nCP)* 2, 1);
            end
            rowCount = 0;
            for frameIndex = 1 : obj.nFrames
                for row = 1: obj.meshSize
                    for col = 1 : obj.meshSize
                        for r = frameIndex - obj.span : frameIndex + obj.span
                            if r == frameIndex
                                % P -C
                                weight = sqrt(obj.cropping);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 1, 1), obj.CaCorner(frameIndex, row, col, 1, 2), ...
                                    (col - 1) * obj.quadWidth + 1, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 2, 1), obj.CaCorner(frameIndex, row, col, 2, 2), ...
                                    col * obj.quadWidth, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 3, 1), obj.CaCorner(frameIndex, row, col, 3, 2), ...
                                    (col - 1) * obj.quadWidth + 1, row * obj.quadHeight);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.CaCorner(frameIndex, row, col, 4, 1), obj.CaCorner(frameIndex, row, col, 4, 2), ...
                                    col * obj.quadWidth, row * obj.quadHeight);
                            elseif r > 0 && r <= obj.nFrames
                                weight = sqrt(obj.smoothness * obj.w_a(frameIndex, r, row, col));
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 1, 1), obj.PaCorner(r, row, col, 1, 2), ...
                                    (col - 1) * obj.quadWidth + 1, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 2, 1), obj.PaCorner(r, row, col, 2, 2), ...
                                    col * obj.quadWidth, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 3, 1), obj.PaCorner(r, row, col, 3, 2), ...
                                    (col - 1) * obj.quadWidth + 1, row * obj.quadHeight);
                                rowCount = rowCount + 2;
                                Ra(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pa(frameIndex, row, col, :, :)), ...
                                    obj.PaCorner(r, row, col, 4, 1), obj.PaCorner(r, row, col, 4, 2), ...
                                    col * obj.quadWidth, row * obj.quadHeight);
                            end
                        end
                
                    end
                end
                % PaCaCP - HPbCaC
                if is1st ~= 1
                    for k = 1:obj.maxppf
                        if obj.validCP(frameIndex, k) == 1
                            weight = 1;
                            pb = obj.CP(frameIndex, k, 3:4);
                            colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                            rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                            pa = obj.CP(frameIndex, k, 1:2);
                            cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                            rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                            rowCount = rowCount + 2;
                            [a2, b2] = obj.transform(obj.CbCP(frameIndex, k, 1:2), obj.H * squeeze(obj.Pb(frameIndex, rowb, colb, :, :)));
                            Ra(rowCount - 1:rowCount) = obj.get1rType1(squeeze(obj.Pa(frameIndex, rowa, cola, :, :)), ...
                                obj.CaCP(frameIndex, k, 1), obj.CaCP(frameIndex, k, 2), a2, b2) * weight;                            
                        end                        
                    end
                end
            end
            Ra(rowCount + 1:length(Ra)) = [];
        end
        
        function Rb = computeRb(obj, is1st)
            % return a vector with the length of number of quadratic terms 
            if is1st
                Rb = zeros(obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4 * 2, 1);
            else
                Rb = zeros((obj.nFrames * obj.meshSize * obj.meshSize * (1 + 2 * obj.span) * 4 + obj.nCP)* 2, 1);
            end
            rowCount = 0;
            for frameIndex = 1 : obj.nFrames
                for row = 1: obj.meshSize
                    for col = 1 : obj.meshSize
                        for r = frameIndex - obj.span : frameIndex + obj.span
                            if r == frameIndex
                                % P -C
                                weight = sqrt(obj.cropping);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 1, 1), obj.CbCorner(frameIndex, row, col, 1, 2), ...
                                    (col - 1) * obj.quadWidth + 1, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 2, 1), obj.CbCorner(frameIndex, row, col, 2, 2), ...
                                    col * obj.quadWidth, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 3, 1), obj.CbCorner(frameIndex, row, col, 3, 2), ...
                                    (col - 1) * obj.quadWidth + 1, row * obj.quadHeight);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.CbCorner(frameIndex, row, col, 4, 1), obj.CbCorner(frameIndex, row, col, 4, 2), ...
                                    col * obj.quadWidth, row * obj.quadHeight);
                            elseif r > 0 && r <= obj.nFrames
                                weight = sqrt(obj.smoothness * obj.w_b(frameIndex, r, row, col));
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 1, 1), obj.PbCorner(r, row, col, 1, 2), ...
                                    (col - 1) * obj.quadWidth + 1, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 2, 1), obj.PbCorner(r, row, col, 2, 2), ...
                                    col * obj.quadWidth, (row - 1) * obj.quadHeight + 1);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 3, 1), obj.PbCorner(r, row, col, 3, 2), ...
                                    (col - 1) * obj.quadWidth + 1, row * obj.quadHeight);
                                rowCount = rowCount + 2;
                                Rb(rowCount - 1:rowCount)= weight * obj.get1rType1(squeeze(obj.Pb(frameIndex, row, col, :, :)), ...
                                    obj.PbCorner(r, row, col, 4, 1), obj.PbCorner(r, row, col, 4, 2), ...
                                    col * obj.quadWidth, row * obj.quadHeight);
                            end
                        end
                    end
                end
                if is1st ~= 1
                    for k = 1:obj.maxppf
                        if obj.validCP(frameIndex, k) == 1
                            weight = 1;
                            pa = obj.CP(frameIndex, k, 1:2);
                            cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                            rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                            pb = obj.CP(frameIndex, k, 3:4);
                            colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                            rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                            rowCount = rowCount + 2;
                            [a2, b2] = obj.transform(obj.CaCP(frameIndex, k, 1:2), obj.H \ squeeze(obj.Pa(frameIndex, rowa, cola, :, :)));
                            Rb(rowCount - 1:rowCount) = obj.get1rType1(squeeze(obj.Pb(frameIndex, rowb, colb, :, :)), ...
                                obj.CbCP(frameIndex, k, 1), obj.CbCP(frameIndex, k, 2), a2, b2) * weight;
                        end                        
                    end
                end
            end
            Rb(rowCount + 1:length(Rb)) = [];
        end
        
        function dadb = get1rType1(obj, P, a1, b1, a2, b2)
            % P [a1, b1, 1] - [a2 b2 1] 
            [a3, b3] = obj.transform([a1, b1], P);
            da = a3 - a2;
            db = b3 - b2;
            dadb = [da; db];
        end
        
        function dadb = get1rType2(obj, Pt, Pr, at, bt, ar, br)
            % Pt [a1 b1 1] - Pr [a2 b2 1] 
            [a3, b3] = obj.transform([at, bt], Pt);
            [a4, b4] = obj.transform([ar, br], Pr);
            da = a3 - a4;
            db = b3 - b4;
            dadb = [da; db];
        end
        
                
        function render(obj, outPath, gap, marker, blend)
            if ~exist(outPath, 'dir')
                mkdir(outPath);
            end
            if obj.useImage
                 
            else
%                 videoA = VideoReader(obj.seqA);
%                 videoB = VideoReader(obj.seqB);                
%                 for i = 1:obj.span
%                     readFrame(videoA); 
%                     readFrame(videoB); 
%                 end
            end
            obj.gap = gap;
            obj.updateOffset();
%             obj.gap = 0;
            parfor frameIndex = 1 + obj.span:obj.nFrames - obj.span %parfor
                disp(['rendering: # ' int2str(frameIndex)]);
                if obj.useImage
                    fileListA = dir(obj.seqA);
                    fileListA = fileListA(3:length(fileListA));
                    fileListB = dir(obj.seqB);
                    fileListB = fileListB(3:length(fileListB));
                    fileNameA = fileListA(frameIndex).name;
                    fileNameB = fileListB(frameIndex).name;
                    IA = imread([obj.seqA fileNameA]);
                    IB = imread([obj.seqB fileNameB]);
                else
                    IA = readFrame(videoA);
                    IB = readFrame(videoB);
                end
                valid = squeeze(obj.CP(frameIndex, obj.validCP(frameIndex, :) == 1, :));
                notvalid = squeeze(obj.CP(frameIndex, obj.validCP(frameIndex, :) == -1, :));
                if marker
                    if size(valid, 2) ~= 1
                        IA = insertMarker(IA, valid(:, 1:2), 'x', 'color', 'green');
                        IB = insertMarker(IB, valid(:, 3:4), 'o', 'color', 'green');
                    end
                    if size(notvalid, 2) ~= 1
                        IA = insertMarker(IA, notvalid(: ,1:2), 'x', 'color', 'red');
                        IB = insertMarker(IB, notvalid(: ,3:4), 'o', 'color', 'red');
                    end
                end
%                 warpA = zeros(obj.videoHeight + 2 * obj.gap, obj.videoWidth + 2 * obj.gap, 3);
%                 warpB = zeros(obj.videoHeight + 2 * obj.gap, obj.videoWidth + 2 * obj.gap, 3);
%                 for i = 1:obj.meshSize
%                     for j = 1:obj.meshSize
%                         warpA = obj.warp1('a', frameIndex, i, j, IA, warpA); 
%                         warpB = obj.warp1('b', frameIndex, i, j, IB, warpB); 
%                     end
%                 end
                warpA = obj.render1('a', IA, frameIndex, obj.Pa, obj.Ca);
                warpB = obj.render1('b', IB, frameIndex, obj.Pb, obj.Cb);
                imwrite(uint8(warpA), [outPath '/A' int2str(frameIndex) '.png']);
                imwrite(uint8(warpB), [outPath '/B' int2str(frameIndex) '.png']);

                warpA = uint8(warpA); warpB = uint8(warpB);

                if ~blend
                    warp = obj.imageBlending(warpA, warpB);
                else
                    warp = mblend3(warpA, warpB, 0.5);
                end
                imwrite(warp, [outPath '/' int2str(frameIndex) '.png']);
            end
        end
        
        function imwarp = render1(obj, ab, I, frameIndex, P, C)
            src = Mesh(obj.videoHeight, obj.videoWidth, obj.quadWidth, obj.quadHeight);
            des = Mesh(obj.videoHeight, obj.videoWidth, obj.quadWidth, obj.quadHeight);
            %
            if ab == 'a'
                HH = obj.Offset;
            else
                HH = obj.Offset * obj.H;
            end
            for i = 0 : obj.meshSize
                for j = 0 : obj.meshSize
                    x = i * obj.quadHeight + 1;
                    y = j * obj.quadWidth + 1;
                    if i == 0 && j == 0
                        B11 = HH * squeeze(P(frameIndex, i+1, j+1, :, :)) ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        des.setVertex(i,j, myPoint(xx11, yy11));
                        continue;
                    end
                    if i == 0 && j == obj.meshSize
                        B10 = HH * squeeze(P(frameIndex, i+1, j, :, :))...
            / squeeze(C(frameIndex, i+1, j, :, :)) ;
                        %B10 = B10;

                        [xx10, yy10] = obj.transform([y x], B10);
                        des.setVertex(i,j, myPoint(xx10, yy10));
                        continue;
                    end
                    if i == 0 && j > 0 && j < obj.meshSize
                        B11 = HH * squeeze(P(frameIndex, i+1, j+1, :, :))  ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        B10 = HH * squeeze(P(frameIndex, i+1, j, :, :))  ...
            / squeeze(C(frameIndex, i+1, j, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        [xx10, yy10] = obj.transform([y x], B10);
                        des.setVertex(i,j, myPoint(mean([xx10 xx11]), mean([yy10 yy11])));
                        continue;
                    end
                    if i>0 && i < obj.meshSize && j == 0
                        B11 = HH * squeeze(P(frameIndex, i+1, j+1, :, :))  ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        B01 = HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        [xx01, yy01] = obj.transform([y x], B01);
                        des.setVertex(i,j, myPoint(mean([xx01 xx11]), mean([yy01 yy11])));
                        continue;
                    end
                    if i > 0 && i < obj.meshSize && j > 0 && j < obj.meshSize
                        B11 = HH * squeeze(P(frameIndex, i+1, j+1, :, :))  ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        B01 = HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        B00 = HH * squeeze(P(frameIndex, i, j, :, :))  ...
            / squeeze(C(frameIndex, i, j, :, :));
                        B10 = HH * squeeze(P(frameIndex, i+1, j, :, :))  ...
            / squeeze(C(frameIndex, i+1, j, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        [xx01, yy01] = obj.transform([y x], B01);
                        [xx10, yy10] = obj.transform([y x], B10);
                        [xx00, yy00] = obj.transform([y x], B00);
                        [xx, yy] = obj.mergepoints([xx01, xx11, xx10, xx00], [yy01, yy11, yy10, yy00]);
                        des.setVertex(i,j, myPoint(xx, yy));
                        continue;
                    end
                    if i>0 && i < obj.meshSize && j == obj.meshSize
                        B00 = HH * squeeze(P(frameIndex, i, j, :, :))  ...
            / squeeze(C(frameIndex, i, j, :, :));
                        B10 = HH * squeeze(P(frameIndex, i+1, j, :, :))  ...
            / squeeze(C(frameIndex, i+1, j, :, :));
                        [xx10, yy10] = obj.transform([y x], B10);
                        [xx00, yy00] = obj.transform([y x], B00);
                        des.setVertex(i,j, myPoint(mean([xx10 xx00]), mean([yy10 yy00])));
                        continue;
                    end
                    if i == obj.meshSize && j == 0
                        B01 = HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        [xx01, yy01] = obj.transform([y x], B01);
                        des.setVertex(i,j, myPoint(xx01, yy01));
                        continue;
                    end
                    if i == obj.meshSize && j>0 && j < obj.meshSize
                        B00 = HH * squeeze(P(frameIndex, i, j, :, :))  ...
            / squeeze(C(frameIndex, i, j, :, :));
                        B01 = HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        [xx00, yy00] = obj.transform([y x], B00);
                        [xx01, yy01] = obj.transform([y x], B01);
                        des.setVertex(i,j, myPoint(mean([xx01 xx00]), mean([yy01 yy00])));
                        continue;
                    end
                    if i == obj.meshSize && j == obj.meshSize
                        B00 = HH * squeeze(P(frameIndex, i, j, :, :)) ...
            / squeeze(C(frameIndex, i, j, :, :));
                        [xx00, yy00] = obj.transform([y x], B00); 
                        des.setVertex(i,j, myPoint(xx00, yy00));
                        continue;
                    end
                end
            end
            imwarp = zeros(obj.videoHeight+obj.gap*2,obj.videoWidth+obj.gap*2,3);

            for i=1:obj.meshSize
                for j=1:obj.meshSize
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
                    imwarp = quadWarp(obj,I,qd1,qd2, imwarp);
                end
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
            
            HH = homography_4pts(source',target');
            HH = HH./HH(3,3);
            
            %qd = Quad(q2.V00,q2.V01,q2.V10,q2.V11);
            imwarp = myWarp(minx,maxx,miny,maxy,double(im),imwarp,HH,obj.gap);
            imwarp = uint8(imwarp);
            
        end
        
        function output_canvas = imageBlending(~, warped_img1,warped_img2)
    
            w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
            w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');

            w1 = mat2gray(w1);
            w2 = mat2gray(w2);

            warped_img1 = double(warped_img1);
            warped_img2 = double(warped_img2);
            output_canvas(:,:,1) = ((warped_img1(:,:,1).*w1)+(warped_img2(:,:,1).*w2))./(w1+w2);
            output_canvas(:,:,2) = ((warped_img1(:,:,2).*w1)+(warped_img2(:,:,2).*w2))./(w1+w2);
            output_canvas(:,:,3) = ((warped_img1(:,:,3).*w1)+(warped_img2(:,:,3).*w2))./(w1+w2);
            output_canvas = uint8(output_canvas);
            
        end
        
        function imwarp = warp1(obj, aorb, frameIndex, location_i, location_j, source, imwarp)
            if aorb == 'a'
                P = obj.Pa; 
                C = obj.Ca;
            else
                P = obj.Pb;
                C = obj.Cb;
            end
            B = squeeze(P(frameIndex, location_i, location_j, :, :)) / squeeze(C(frameIndex, location_i, location_j, :, :));
            if aorb == 'b'
                B = obj.H * B; 
            end
            B = obj.Offset * B;
            %B = B ./B(3,3);
            minx = (location_i - 1)*obj.quadHeight+1;
            maxx = minx + obj.quadHeight - 1;
            miny = (location_j - 1)*obj.quadWidth+1;
            maxy = miny + obj.quadWidth - 1;
            [x00, y00 ] = obj.transform([miny, minx], B);
            [x01, y01 ] = obj.transform([miny, maxx], B);
            [x10, y10 ] = obj.transform([maxy, minx], B);
            [x11, y11 ] = obj.transform([maxy, maxx], B);
            minx = min(x00, x01);minx = min(minx, x10);minx = min(minx, x11);
            miny = min(y00, y10);miny = min(miny, y01);miny = min(miny, y11);
            maxx = max(x10, x11);maxx = max(maxx, x01);maxx = max(maxx, x00);
            maxy = max(y01, y11);maxy = max(maxy, y10);maxy = max(maxy, y00);
            % imwarp = myWarp(minx, maxx, miny, maxy, double(source), obj.imwarp(frameIndex).cdata, inv(B), obj.gap);
            imwarp = myWarp(minx, maxx, miny, maxy, double(source), imwarp, inv(B), obj.gap);
        end
   
        
        
        function calcOmega(obj)
            disp('computing omega...')
            for i = 1:obj.meshSize
                for j = 1:obj.meshSize
                    for t = 1:obj.nFrames                        
                        for r = t-obj.span:t+obj.span
                            if r > 0 && r < obj.nFrames
                                dPa = abs(obj.Pa(t,i,j,1,3) - obj.Pa(r,i,j,1,3)) + abs(obj.Pa(t,i,j,2,3) - obj.Pa(r,i,j,2,3));
                                dPb = abs(obj.Pb(t,i,j,1,3) - obj.Pb(r,i,j,1,3)) + abs(obj.Pb(t,i,j,2,3) - obj.Pb(r,i,j,2,3));
                                obj.w_a(t,r,i,j) = gaussmf(abs(t-r), [10 0]) * gaussmf(dPa, [200 0]);
                                obj.w_a(t,t,i,j) = 0;
                                obj.w_b(t,r,i,j) = gaussmf(abs(t-r), [10 0]) * gaussmf(dPb, [200 0]);
                                obj.w_b(t,t,i,j) = 0;
                            end
                        end
                    end
                end
            end            
        end
        
        function [x,y] = transform(~, xxyy, B)
            xx = xxyy(1); yy = xxyy(2);
            res = B * [xx;yy;1];
            x = res(1)/res(3);
            y = res(2)/res(3);
        end
        
        function CaPa = timesCa(obj, Pa)
            %
            x = Pa(1); y = Pa(2); frameIndex = Pa(3);
            if (x * y) ~= 0
                i = floor((x - 0.0001) / obj.quadWidth) + 1;
                j = floor((y - 0.0001) / obj.quadHeight) + 1;
                [x, y] = obj.transform(Pa, squeeze(obj.Ca_inv(frameIndex, i, j, :, :)));
            end
            CaPa = ones(3, 1);
            CaPa(1) = x; CaPa(2) = y;
        end
        
        function CaPa = timesCb(obj, Pa)
            %
            x = Pa(1); y = Pa(2); frameIndex = Pa(3);
            if (x * y) ~= 0
                i = floor((x - 0.0001) / obj.quadWidth) + 1;
                j = floor((y - 0.0001) / obj.quadHeight + 1);
                [x, y] = obj.transform(Pa, squeeze(obj.Cb_inv(frameIndex, i, j, :, :)));
            end
            CaPa = ones(3, 1);
            CaPa(1) = x; CaPa(2) = y;
        end
        
    end
end


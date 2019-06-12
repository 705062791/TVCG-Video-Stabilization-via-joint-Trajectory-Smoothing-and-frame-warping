classdef TrackList < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nTrack;
        list;
        wSize;
        head;
        tail;
        id;
        length;
        points;
    end
    
    methods
        function obj = TrackList(nTrack, list)
            MAXLENGTH = 100;
            obj.nTrack = nTrack;
            obj.list = list;
            obj.head = zeros(nTrack, 1);
            obj.id = zeros(nTrack, 1);
            obj.tail = obj.head;
            obj.length = obj.head;
            obj.points = zeros(MAXLENGTH, 2, nTrack);
            for i = 1:nTrack
                obj.id(i) = obj.list{i}.id;
                obj.head(i) = obj.list{i}.head;
                obj.tail(i) = obj.list{i}.tail;
                obj.length(i) = obj.list{i}.length;                
                obj.points(1:obj.length(i), :, i) = obj.list{i}.p;
            end
        end
        
        
        function setWindowSize(obj, size)
            obj.wSize = size; 
        end
        
        function [M, M2, id] = getM(obj, index)
            left = index * obj.wSize / 2;
            right = left + obj.wSize - 1;
            mid = (left + right) / 2.0;
            select = obj.head < mid & obj.tail > mid & obj.length > obj.wSize;
            nSelect = sum(select);
            M = zeros(obj.wSize * 2, nSelect);
            id = zeros(nSelect, 1);
            count = 0;
            for i = 1:obj.nTrack
                if select(i)
                    count = count + 1;
                    id(count) = obj.id(i);
                    if obj.head(i) < left && obj.tail(i) > right
                        padleft = left - obj.head(i);                        
                        M(:, count) = reshape(squeeze(obj.points(1 + padleft : padleft + obj.wSize, :, i))', [obj.wSize * 2 1]);                        
                    end
                    if obj.head(i) >= left && obj.tail(i) > right
                        padleft = obj.head(i) - left;
                        padright = obj.tail(i) - right;
                        M(padleft * 2 + 1:obj.wSize * 2, count) = reshape(squeeze(obj.points(1:obj.wSize - padleft, :, i))', [(obj.wSize - padleft) * 2 1]);
                    end
                    if obj.head(i) < left && obj.tail(i) <= right
                        padleft = left - obj.head(i); 
                        padright = right - obj.tail(i);
                        M(1:(obj.wSize - padright) * 2 ,count) = reshape(squeeze(obj.points(padleft + 1 : obj.length(i), :, i))', [(obj.wSize - padright) * 2 1]);
                    end
                    if obj.head(i) >= left && obj.tail(i) <= right
                        padleft = obj.head(i) - left;
                        padright = right - obj.tail(i);
                        M(padleft * 2 + 1 : (obj.wSize - padright) * 2) = reshape(squeeze(obj.points(1:obj.length(i), :, i))', [obj.length(i) * 2 1]);                       
                    end
                end
            end
            
            M3 = reshape(M, [2 obj.wSize nSelect]);
            M2 = zeros(2, nSelect, obj.wSize);
            M2(1, :, :) = squeeze(M3(1, :, :))';
            M2(2, :, :) = squeeze(M3(2, :, :))';
        end
        
        function [f1, f2] = getF(obj, frameIndex)
            select = obj.head < frameIndex & obj.tail >= frameIndex;
            N = sum(select);
            f1 = zeros(N, 2);
            f2 = f1;
            count = 0;
            for i = 1:obj.nTrack                
                if select(i) == 1
                    count = count + 1;
                    f1(count, :) = obj.points(frameIndex - obj.head(i), :, i);
                    f2(count, :) = obj.points(frameIndex - obj.head(i) + 1, :, i);
                end
            end
%             [~, inlier] = EstimateHomographyByRANSAC(f1', f2', 0.0002);
%             f1 = f1(inlier, :);
%             f2 = f2(inlier, :);
        end
        
        
        
    end
    
end


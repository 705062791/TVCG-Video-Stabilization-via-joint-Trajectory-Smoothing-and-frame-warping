classdef Track < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id;
        length;
        head;
        tail;
        p;
    end
    
    methods
        function obj = Track(id, length, head, point_str)
            obj.id = id;
            obj.length = length;
            obj.head = head;
            obj.p = zeros(length, 2);
            points = textscan(point_str, '%f %f');
            obj.p = [points{1} points{2}];
            obj.tail = obj.head + obj.length - 1;
        end
    end
    
    
end


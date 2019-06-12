function [ ] = save_path( tracks,save_file,min_length,casename )
%SAVE_PATH 此处显示有关此函数的摘要
%   此处显示详细说明
addpath(genpath('D:/c++/Stitching-1.0.0'));


point_x=tracks.points(:,:,1);
point_y=tracks.points(:,:,2);
num_of_path=tracks.nTrack;
frame_num=tracks.nFrame;
path=[];
%%
for i=1:num_of_path
    start=0;
    
    for t=1:frame_num
        
        if point_x(i,t)~=0
          start=t;
          break;
        end
        
    end
    
    single_path_x=point_x(i,:);
    single_path_x=single_path_x(single_path_x~=0);
    
    single_path_y=point_y(i,:);
    single_path_y=single_path_y(single_path_y~=0);
    
    end_frame=start+max(size(single_path_x))-1;
    
    single_path=[[start,end_frame];[single_path_x',single_path_y']];
    
    length_of_single_path = size(single_path); 
    
    if length_of_single_path(1,1)>=min_length
       path=[path;single_path];
    end
    
     clear single_path_x;
     clear single_path_y;
    
end
%%
save_file=strcat(save_file,casename,'_path','.txt');
save(save_file,'path','-ascii');
% save D:/data/video33_path.txt path -ascii;
display('over');

end


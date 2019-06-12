casename={'video2','video11','video14','video22','video26','video27','video28','video29','video31','video33','video38',...
    'video39','video40','video41','video42','video43','video45','video46','video47','video48','video53',...
    'video54','video55','video56','video70','video71','video72','video73','video78','video80','video81'};

dbstop error;

addpath(genpath('Stitching-1.0.0'));
%%
%参数和路径
%视频所在路径
computeVideoStabilityHandler = @computeVideoStabilityWindow;
caseFile='D:/data/video_test/';
saveFile = 'd:/data/';
min_length = 50;

for i=3:3
  %%
  videoFile=strcat(caseFile,'/',casename{1,i},'.mp4');
  videotraFile=strcat(caseFile,'/',casename{1,i},'.mat');

%参数
caseName=casename{1,i};

nWinSize = 40;
trackNum = 700;
divNum = 10;
%%
%运行代码
track=GetTracksFromVideo(videotraFile, videoFile, divNum, trackNum);
save_path( track,saveFile,min_length,caseName )
end


% Video Stitching
% Written by Tan Su
% tsu@cse.cuhk.edu.hk
addpath('mesh');
addpath('Path');
addpath('RANSAC');
addpath('stitch');
addpath('blend');
addpath('tracks');
addpath('tracks/helpers');
addpath('graph');
addpath('peter');
%% Set up vlfeat
poolobj = gcp('nocreate');
delete(poolobj);
cd vlfeat-0.9.19/toolbox;
feval('vl_setup');
cd ../..;

%% params
UseImageSeq = 1;
data = '../case21/';
input_A = 'left/';
input_B = 'right/';

PointsPerFrame = 1000;
MeshSize = 8;
MaxIte = 30;
Smoothness = 1;
Cropping = 1;
Stitchness = 20;
OutputPath = 'res_20_2stage';



%% Track by KLT
if ~exist([data 'tracks1000.mat'], 'file')
    trackA = GetTracks([data input_A], 10, 1000);
    trackB = GetTracks([data input_B], 10, 1000);    
    save([data 'tracks1000.mat'], 'trackA', 'trackB');
    
else
    load([data 'tracks1000.mat']);
end
%% Clustering trajectories
% if ~exist([data 'tracks_34_500.mat'], 'file')
%     trackA.addLabel(34);
%     trackB.addLabel(34);
%     if trackA.nLabel > trackB.nLabel
%         trackB.nLabel = trackA.nLabel;
%     else
%         trackA.nLabel = trackB.nLabel;       
%     end
% 
%     save([data 'tracks_34_500.mat'], 'trackA', 'trackB');
% else
%     load([data 'tracks_34_500.mat']);
% end
%% print the label
% PrintLabel([data input_A], trackA, [data '/left_label/']);
% PrintLabel([data input_B], trackB, [data '/right_label/']);

%% Matching SIFT in every frame pair
if ~exist([data 'ControlPoints1000.mat'], 'file')
    [CP, ppf] = getControlPoints([data input_A], [data input_B], PointsPerFrame);
    save([data 'ControlPoints1000.mat'], 'CP', 'ppf');
else
    load([data 'ControlPoints1000.mat']);
end

%% Get common background
% if ~exist([data 'graph500.mat'], 'file')
%     alpha = 0.01;
%     beta = 0.01;
%     maxNlabel = max(trackA.nLabel, trackB.nLabel);
%     trackA.nLabel = maxNlabel; trackB.nLabel = maxNlabel;
%     [path, graph, goodA, goodB] = GetGraph(trackA, trackB, CP, ppf, alpha, beta);
%     backListA = refineTrack(trackA, goodA);
%     backListB = refineTrack(trackB, goodB);
%     save([data 'graph500.mat'], 'path', 'graph', 'backListA', 'backListB');
% % PrintBackground([data '/left/'], trackA, backListA, [data '/left_back/']);
% % PrintBackground([data '/right/'], trackB, backListB, [data '/right_back/']);
% % PrintGoodBackground([data '/left/'], trackA, goodA, [data '/left_back/']);
% % PrintGoodBackground([data '/right/'], trackB, goodB, [data '/right_back/']);
% else
%     load([data 'graph500.mat']);
% end
%% Compute original camera path (by As-similar-as-possible Warping)
if ~exist([data 'Path8_track_new.mat'], 'file')    
    [pathA] = getPath([data input_A], MeshSize, trackA);
    [pathB] = getPath([data input_B], MeshSize, trackB);
    save([data 'Path8_track_new.mat'], 'pathA', 'pathB');  
else 
    load([data 'Path8_track_new.mat']);
end

% if ~exist([data 'ControlPoints1000_refine.mat'], 'file')
%     [CP, ppf] = refineCP(CP, ppf, backListA, backListB, trackA, trackB);
%     save([data 'ControlPoints1000_refine.mat'], 'CP', 'ppf');
% else
%     load([data 'ControlPoints1000_refine.mat']);
% end



%% Optimize the paths
stitcher = VideoStitch2([data input_A], [data input_B], pathA, pathB, CP, ppf, UseImageSeq, Smoothness, Cropping, Stitchness);
stitcher.init();
stitcher.optPath(MaxIte);
%stitcher.optPa(1);
%stitcher.optPb(1);
stitcher.render([data OutputPath], 500, 0, 0);
save([data 'stitcher.mat'], 'stitcher');

%% Evaluation Score
padding = 30;
CP_ = stitcher.getStitchedCP(padding + 1, stitcher.nFrames - padding);
save([data 'stitchedCP.mat'], 'CP_');
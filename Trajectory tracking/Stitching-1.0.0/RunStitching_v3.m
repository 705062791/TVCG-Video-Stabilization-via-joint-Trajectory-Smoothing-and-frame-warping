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
%% Set up vlfeat
poolobj = gcp('nocreate');
delete(poolobj);
cd vlfeat-0.9.19/toolbox;
feval('vl_setup');
cd ../..;

%% params
UseImageSeq = 1;
data = '../case3/';
input_A = 'left/';
input_B = 'right/';
maskA = '/maskleft/';
maskB = '/maskright/';

PointsPerFrame = 2000;
MeshSize = 8;
MaxIte = 30;
Smoothness = 1;
Cropping = 1;
Stitchness = 20;
OutputPath = 'res_20_2stage';



%% Track by KLT
if ~exist([data 'tracks500.mat'], 'file')
    trackA = GetTracks([data input_A], 10, 500);
    trackB = GetTracks([data input_B], 10, 500);    
    save([data 'tracks500.mat'], 'trackA', 'trackB');
    
else
    load([data 'tracks500.mat']);
end
%% Matching SIFT in every frame pair
if ~exist([data 'ControlPoints1000.mat'], 'file')
    [CP, ppf] = getControlPoints([data input_A], [data input_B], PointsPerFrame);
    save([data 'ControlPoints1000.mat'], 'CP', 'ppf');
else
    load([data 'ControlPoints1000.mat']);
end
%% Compute original camera path (by As-similar-as-possible Warping)
if ~exist([data 'Path8_track_new.mat'], 'file')    
    [pathA] = GetPathFromMask([data input_A], MeshSize, trackA, [data '/maskleft/']);
    [pathB] = GetPathFromMask([data input_B], MeshSize, trackB, [data '/maskright/']);
    save([data 'Path8_track_new.mat'], 'pathA', 'pathB');  
else 
    load([data 'Path8_track_new.mat']);
end

if ~exist([data 'ControlPoints2000_refine.mat'], 'file')
    [CP, ppf] = RefineCPFromMask(CP, ppf, [data '/maskleft/'], [data '/maskright/']);
    save([data 'ControlPoints2000_refine.mat'], 'CP', 'ppf');
else
    load([data 'ControlPoints2000_refine.mat']);
end



%% Optimize the paths
stitcher = VideoStitch2([data input_A], [data input_B], pathA, pathB, CP, ppf, UseImageSeq, Smoothness, Cropping, Stitchness);
stitcher.init();
stitcher.optPath(MaxIte);
%stitcher.optPa(1);
%stitcher.optPb(1);
stitcher.render([data OutputPath], 500, 1, 0);
save([data 'stitcher.mat'], 'stitcher');

%% Evaluation Score
padding = 30;
CP_ = stitcher.getStitchedCP(padding + 1, stitcher.nFrames - padding);
save([data 'stitchedCP.mat'], 'CP_');
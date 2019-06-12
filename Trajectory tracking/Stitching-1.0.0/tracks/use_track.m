data = '../case17/';
system(['MyDenseTrajectories.exe ' data 'left_sd.mp4 -W 40 -L 40 -E 499']);
tracksfile = '../case17/left_sd_tracks.txt';

file = fopen(tracksfile, 'r');
nTrack_str = fgetl(file);
nTrack = sscanf(nTrack_str, '%d');
track_list = cell(nTrack, 1);
for trackIndex = 1:nTrack
    track_info = sscanf(fgetl(file), '%d');
    id = track_info(1);
    leng = track_info(2);
    start = track_info(3);
    line = fgetl(file);
    track_list{trackIndex} = Track(id, leng, start, line);
end
tracks = TrackList(nTrack, track_list);
fclose(file);
window_size = 40;
tracks.setWindowSize(window_size);
% [M, M2, ID] = tracks.getM(10);
[f1, f2] = tracks.getF(499);




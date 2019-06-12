video = VideoReader('A.mp4');
k = 0;
mkdir('res');
while (hasFrame(video))
    k = k + 1;
    frame = readFrame(video);
    filename = ['./res/' sprintf('%03d',k) '.png'];
    imwrite(frame, filename);
    %system(['Hsaliency ' filename ' ./SANY0025/']);
end
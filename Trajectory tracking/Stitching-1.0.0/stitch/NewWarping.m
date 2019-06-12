function [ homos ] = NewWarping( pa, pb, H, W, qH, qW, lambda)
% A warpped up wersion of as-similar-as-possible warping, with pre-warp. 
% 
    nP = length(pa);
    if length(pb) ~= nP
        error('Points Numbers met Matching!');
    end
    
%     x1 = [pa'; ones(1,length(pa))];
%     x2 = [pb'; ones(1,length(pb))];
%     [x1, T1] = NormalisePoints(x1);
%     [x2, T2] = NormalisePoints(x2);
%     preH = DLT(x1, x2);
%     preH = T2\preH*T1;
    [preH, ~] = ransacfithomography(pa', pb', 0.001);
    whilecount = 0;
    while isnan(preH(1, 1)) && whilecount < 100
        [preH, ~] = ransacfithomography(pa', pb', 0.001);
        whilecount = whilecount + 1;
    end
    if whilecount == 100
        error('?') ;
    end
%     preH = eye(3);
    pbWarp = preH \ [pb' ; ones(1, nP)];
    pbWarp(1, :) = pbWarp(1, :) ./ pbWarp(3, :);
    pbWarp(2, :) = pbWarp(2, :) ./ pbWarp(3, :);
    pbWarp = pbWarp(1:2, :)';
    
    
%     diff = sum((pa - pbWarp) .* (pa - pbWarp), 2) < 1000;    
%     valid = pbWarp(:, 1) > 0 & pbWarp(:, 1) < W & pbWarp(:, 2) > 0 & pbWarp(:, 2) < H;
%     pa = pa(valid, :);
%     pbWarp = pbWarp(valid, :);
    
    asap = AsSimilarAsPossibleWarping(H, W, qW, qH, lambda);
    asap.SetControlPts(pa, pbWarp);
    asap.Solve();
    e = asap.CalcError();
    
%     grid = asap.Warp(ones(H, W, 3) * 255, 400);
%     imshow(grid);
    
    if e > 10
        disp('?') ;
    end
    homos2 = asap.CalcHomos();
    homos = homos2;
    for row = 1:H/qH
        for col = 1:W/qW
            tempH = preH * squeeze(homos2(row, col, :, :));
            homos(row, col, :, :) = tempH ./ tempH(3, 3);            
        end
    end
    
%     asap = AsSimilarAsPossibleWarping(H, W, W, H, lambda);
%     asap.SetControlPts(pa, pbWarp);
%     asap.Solve();
%     homos2 = asap.CalcHomos();
%     homos = zeros(round(H/qH), round(W/qW), 3, 3);
%     for row = 1:H/qH
%         for col = 1:W/qW
%             homos(row, col, :, :) = preH * squeeze(homos2(1, 1, :, :));
%             homos(row, col, :, :) = homos(1, 1, :, :) ./ homos(1, 1, 3, 3);
%         end
%     end
end


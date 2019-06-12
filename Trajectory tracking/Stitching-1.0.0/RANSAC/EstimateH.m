function H = EstimateH(kp1, kp2)
    M     = 500;  % Number of hypotheses for RANSAC.
    thr   = 0.1;  % RANSAC threshold.
	data_orig = [ kp1(1:2,:) ; ones(1,size(kp1,2)) ; kp2(1:2,:) ; ones(1,size(kp1,2)) ];
    [ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
	[ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
	data_norm = [ dat_norm_img1 ; dat_norm_img2 ];
    
    rng(0);
    [ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);
    con = sum(res<=thr);
    [ ~, maxinx ] = max(con);
    inliers = res(:,maxinx)<=thr;
    
    [ h,~,~,~ ] = feval(fitfn,data_norm(:,inliers));
    H = T2\(reshape(h,3,3)*T1);
    
end
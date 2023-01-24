
%% Computation of HOMOGRAPHY using RANSAC
%% References:
% # Digital Video Processing (http://10.6.4.152/dvp/dvp.html),
% Computer Sc., Indian Institute of Technology, Madras.
% # VLFEAT SIFT Tool box (http://www.vlfeat.org/overview/sift.html)
% # RANSAC algorithm with example of finding homography : Edward Wiggin.
% MATLAB Central 2011.

function H = RANSAC(sift_match_points_f1, sift_match_points_f2)
lmax = size(sift_match_points_f1, 2);
C = 0;

while C < 0.85*lmax 
% Picking the Random points (4 Nos)
    ind = randIndex(lmax, 4);    
    pts1 = sift_match_points_f1(:,ind);
    pts2 = sift_match_points_f2(:,ind);
    
% Compute the HOMOGRAPHY
    H = solveHomo(pts1,pts2);
    
% Compute the residual for rest of points (lmax - 4)
    remain_pts = H*[sift_match_points_f1;ones(1,lmax)];
    norm_remain_pts = [remain_pts(1,:)./remain_pts(3,:);...
                      remain_pts(2,:)./remain_pts(3,:);ones(1,lmax)];
   
    distance_ransac  = norm_remain_pts - [sift_match_points_f2;...
                                                         ones(1,lmax)];
                                                     
    magnitude_distance = sqrt(sum(distance_ransac.^2,1));
    myRansacLogic = magnitude_distance < 10;
    C = size(find(myRansacLogic),2);
    
end


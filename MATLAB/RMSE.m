function [ rmse ] = RMSE( x, pts )
%RMSE Calculate the root-mean-squared error between the mean particle
%estimate and the true values

[d, K] = size(x);

rmse.pos = zeros(1, K);
rmse.vel = zeros(1, K);

pts_x = cell2mat(permute(pts, [3,2,1]));
est_x = mean(pts_x, 3);
error = bsxfun(@minus, est_x, x);

rmse.pos = sqrt(sum(error(1:2,:).^2,1));
rmse.vel = sqrt(sum(error(3:4,:).^2,1));

rmse.mean_pos = sqrt(mean(rmse.pos.^2));
rmse.mean_vel = sqrt(mean(rmse.vel.^2));

end


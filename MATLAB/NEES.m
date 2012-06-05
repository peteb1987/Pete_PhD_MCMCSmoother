function [ nees, nees_over_time ] = NEES( x, pts )
%NEES Calculate the normalised-estimation-error-squared between the mean
% particle estimate and the true values

[d, K] = size(x);

Np = size(pts, 1);

pts_x = cell2mat(permute(pts, [3,2,1]));
est_x = mean(pts_x, 3);
est_P = zeros(d,d,K);
error = zeros(1,K);
for kk = 1:K
%     est_P(:,:,kk) = cov(squeeze(pts_x(:,kk,:))');
    innov = bsxfun(@minus, squeeze(pts_x(:,kk,:)), x(:,kk));
    est_P(:,:,kk) = (1/Np)*innov*(innov');
    if all(all(~diff(innov')))
        error(kk) = 1;
    else
        error(kk) = ((est_x(:,kk)-x(:,kk))'/est_P(:,:,kk))*(est_x(:,kk)-x(:,kk));
    end
end

nees = sum(error)/K;
nees_over_time = error;

end


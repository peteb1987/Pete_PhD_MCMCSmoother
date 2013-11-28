function plot_banana( fig, xy_mean, rng_var, bng_var, c )
%PLOT_BANANA Plot equal probability contour for a range-bearing observation model

% Number of points to plot (must be divisible by 2)
Npts = 2000;

% Convert to polar
[b_mean, r_mean] = cart2pol(xy_mean(1), xy_mean(2));

% Find range of... range
r_min = r_mean-c*sqrt(rng_var);
r_max = r_mean+c*sqrt(rng_var);

% Create an array of points
rb_pts = zeros(2,Npts);

% Grid of range points
rb_pts(1,1:Npts/2) = (r_min:(r_max-r_min)/(Npts/2-1):r_max);
rb_pts(1,Npts/2+1:end) = fliplr(r_min:(r_max-r_min)/(Npts/2-1):r_max);

% Calculate corresponding bearing points
rb_pts(2,1:Npts/2) = b_mean + sqrt(abs(bng_var*(c^2-(rb_pts(1,1:Npts/2)-r_mean).^2/rng_var)));
rb_pts(2,Npts/2+1:end) = b_mean - sqrt(abs(bng_var*(c^2-(rb_pts(1,Npts/2+1:end)-r_mean).^2/rng_var)));

% Convert back to Cartesian
xy_pts = zeros(size(rb_pts));
[xy_pts(1,:), xy_pts(2,:)] = pol2cart(rb_pts(2,:), rb_pts(1,:));

assert(all(all(isreal(xy_pts))));

% Plot
figure(fig)
plot(xy_pts(1,:), xy_pts(2,:), ':r')

end


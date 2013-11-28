function [ y_mn ] = tracking_h( model, x )
%nlng_h Deterministic observation function for the nonlinear non-Gaussian
%benchmark model.

y_mn = zeros(model.do,1);

% Bearing
y_mn(1) = atan2(x(2), x(1));

% Elevation
hoz_range = hypot(x(1), x(2));
y_mn(2) = atan2(x(3), hoz_range);

% Range
range = sqrt( x(1)^2 + x(2)^2 + x(3)^2 );
y_mn(3) = range;

% Range Rate
r = x(1:3);
v = x(4:6);
y_mn(4) = dot(r,v)/range;

end

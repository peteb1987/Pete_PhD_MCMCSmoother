function [ H ] = tracking_obsjacobian( model, x )
%tracking_obsjacobian Calculate observation function jacobian for the
%tracking model.

rng = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
hoz_rng_sq = (x(1)^2+x(2)^2);
hoz_rng = sqrt(hoz_rng_sq);

dh1_dx = [-x(2), x(1), 0, 0, 0, 0]/hoz_rng_sq;
dh2_dx = [-x(1)*x(3)/hoz_rng, -x(2)*x(3)/hoz_rng, hoz_rng, 0, 0, 0]/(rng^2);
dh3_dx = [x(1), x(2), x(3), 0, 0, 0]/rng;
dh4_dr = [x(4)*(x(2)^2+x(3)^2)-x(1)*(x(5)*x(2) + x(6)*x(3)),...
          x(5)*(x(1)^2+x(3)^2)-x(2)*(x(4)*x(1) + x(6)*x(3)),...
          x(6)*(x(1)^2+x(2)^2)-x(3)*(x(4)*x(1) + x(5)*x(2))]/(rng^3);
dh4_dv = [x(1), x(2), x(3)]/rng;

H = [dh1_dx; dh2_dx; dh3_dx; [dh4_dr dh4_dv]];
     
end


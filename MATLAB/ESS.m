function ess = ESS( weight )
%CALCESS Calculates effective sample size from array of (logarithmic) weights
weight = exp(weight);
ess = 1/( sum( weight.^2 ) );
end


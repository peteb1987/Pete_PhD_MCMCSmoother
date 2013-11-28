function [ weight ] = normalise_weights( weight )
%NORMALISE_WEIGHTS Normalise log-weights

weight = weight - max(weight);
weight = weight - logsumexp(weight,2);

end


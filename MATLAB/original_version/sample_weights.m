function idx = sample_weights(w)
%SAMPLE_WEGHTS Sample an index from a column vector of weights

% Normalise weights
sumw = sum(w);
if sumw ~= 1
    w = w./sumw;
end

% Sample
edges = [0; cumsum(w)];
[~, idx] = histc(rand,edges);

end
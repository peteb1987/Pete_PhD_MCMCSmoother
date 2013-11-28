function f = loggausspdf(X, Mu, Sigma)
%MVNPDF Calculates the log of a Gaussian pdf
%
%   Y = MVNPDF(X,MU,SIGMA) returns a vector of Gaussian density
%   values. X is a DxN array of points. Mu is a DxN array of means. Sigma
%   is a DxDxN array of covariance matrices or a single DxD covariance marix.
%
%   Based on mvnpdf from the stats toolbox with all the error checking and
%   special cases stripped out. Also, the vectors are the right way round.

% Constant
log2pi = 1.83787706640935;

% Get size of data
[d,n] = size(X);

if n == 1
    
    X0 = X - Mu;
    R = chol(Sigma);
    logSqrtDetSigma = sum(log(diag(R)));
    xRinv = R' \ X0;
    
else
    
    % Center data
    X0 = bsxfun(@minus, X, Mu);
    
    % Cope with only one covariance matrix
    if size(Sigma, 3) == 1
        sameCov = true;
        R = chol(Sigma);
        logSqrtDetSigma = sum(log(diag(R)));
    else
        sameCov = false;
        % Create arrays
        xRinv = zeros(d,n,superiorfloat(X0,Sigma));
        logSqrtDetSigma = zeros(1,n,class(Sigma));
    end
    
    % Loop through points
    for i = 1:n
        
        if ~sameCov
            R = chol(Sigma(:,:,i));
            logSqrtDetSigma(i) = sum(log(diag(R)));
        end
        
        % Transform to standard normal
        xRinv(:,i) = R' \ X0(:,i);
        
    end
    
end

% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2, 1);

f = (-0.5*quadform - logSqrtDetSigma - d*log2pi/2)';

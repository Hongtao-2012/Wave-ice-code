function [nrmse] = funNRMSE(x, xEst)
%[nrmse] = funNRMSE(x, xrec)
% INPUT:
% x - analytical or simulated signal
% xEst - Estimated signal
% Calculate Normalized root mean square error (NRMSE) using amplitude
% range as denominator
% OUTPUT:
% nrmse - a scalar, normalized root mean squared error

% nrmse = rms( x - xEst ) / (max(x) - min(x));
x = x(:);
xEst = xEst(:);
nrmse = sqrt( sum( (x - xEst).^2 )/numel(x)) / (max(x) - min(x));

end


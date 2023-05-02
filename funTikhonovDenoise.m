function y = funTikhonovDenoise(x, transp_flag, dt, alpha)
% y = funTikhonovDenoise(x, transp_flag, dt, alpha)
%-------------------------------------------
% INPUT:
% x - that should be inversely determined. Need not to be explicitly
% transp_flag - Need not to be explicitly input
% dt - time increment interval
% alpha - regularization parameter. If it is too large, original (or clean
% signal) that should be found by denoising will be greatly attenuated
% Determine alpha by means of L-curve method
% OUTPUT:
% y - has a size of n x 1
%-------------------------------------------
% Example
% alpha = 1.0155e-4; % regularization parameter
% maxiter = 8000; % Maximum number of interations
% y = y(:);
% np = length(y);
% ystack = [y(1:end-2); zeros(np-2,1)];
% Fs = 100; % Sampling frequency in Hz
% dt = 1/Fs; % time increment in seconds
% yTik = lsqr(@(x, trans_flag)funTikhonovDenoise(x, trans_flag, dt, alpha), ystack, [], maxiter);
%-------------------------------------------
% Coded by Hongtao, 19-12-2019
%% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%-------------------------------------------
% References:
% [1] MUELLER, Jennifer L.; SILTANEN, Samuli. Linear and nonlinear inverse problems with practical applications. Siam, 2012.

% see Chapter 5 and section 5.5 in Mueller and Siltanen (2012)
switch transp_flag
    case 'notransp'
    % Perform  Atilde * f multiplication, where Atilde has a size of 2(n-2) x n and 
    % f has a size of n x 1
    n = length(x);
    m = 2*(n-2);
    y = zeros(m,1);
    % y(1:n-2) = eye(n-2,n-2) * x(1:n-2)
    y(1:n-2) = x(1:n-2);
    % y(n-1) = y(n) = 0
    
    % Perform L2 * f, where L2 is the second-order difference operator
    % L2 has a size of (n-2) x n
    y(n-1:m) = sqrt(alpha) /(dt^2) * (x(3:end) - 2* x(2:end-1) + x(1:end-2));
    case 'transp'
        % Here x = Atilde * f, where x has a size of 2(n-2) x 1
        % Atilde has a size of (2n-4) x n, f has a dimension of n x 1
        
        % Perform multiplication of  tranpose of Atilde wtih x
        % x has a dimension of (2n-4), or 2(n-2)
        % Atilde' has a size of n x (2n-4), where first (n-2) columns
        % correspond to identity matrix which ends with two rows of zeros
        % and the last (n-2) columns correspond to transpose of L2
       m = length(x); 
       n = (m+4)/2;
       u = x(1:n-2);
       v = x(n-1:m);
    
       y = zeros(n,1);
       % Peform eye(n-2,n-2) * u(n-2,1)
       y(1:n-2) = u(1:n-2);
       
       % Peform transpose(L2) * v
       z = zeros(n,1);
       z(1) = v(1);
       z(2) = -2 * v(1) + v(2);
       z(3:end-2) = v(1:end-2) - 2 * v(2:end-1) + v(3:end);
       z(end-1) = -2 * v(end) + v(end-1);
       z(end) = v(end);
       z = sqrt(alpha) / (dt^2) * z;
       
       % Sum up
       y = y + z;
    otherwise 
        error('error');
end
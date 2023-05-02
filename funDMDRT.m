function [b, dmdF, dmdAmp, dmdFFTAmp, dmdFR,  dmdAmpR,  xdmd, dampF] = funDMDRT(xN, rr, Fs)
%  [b, dmdF, dmdAmp, xdmd, dampF] = funDMDAmp(xN, Fs)
%  INPUT: 
% xN - scalar time series with size N x 1
% Fs - Sampling rate in Hz
% OUTPUT:
% b - coordinates of time series at first time instant on eigenvector
% basis, with size rr x 1
% dmdF - Frequency (Hz) with size rr x 1 in Amplitude spectrum
% dmdAmp - Amplitude with size rr x 1 in Amplitude spectrum
% dmdFFTAmp  - Amplitude with size rr x 1 in Amplitude spectrum
% dmdFR - Frequency (Hz) with size rOpt x 1 in Amplitude spectrum, rOpt
% is determined by using algorithm proposed by GAVISH and DONOHO (2014)
% dmdAmpR - Amplitude with size rOpt x 1 in Amplitude spectrum, rOpt
% is determined by using algorithm proposed by GAVISH and DONOHO (2014)
% xdmd - Reconstructed signal with components that have largest amplitude,
% xdmd has a size of (N-stack-1) x N, corresponding with dmdF and dmdAmp
% dmdF - damping factor with a size of rr x 1, corresponding with dmdF and
% dmdAMp
% Other parameters determined in this function:
% rr - is determined by TSVD of D1
% rOpt - is determined via denoising routine used by this function
% stack - is determined as ceil(N/2)
%------------------------------------------
% Coded by Hongtao, 17-12-2019
% Updated by Hongao, 13-01-2020
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%------------------------------------------
% References:
% [1] Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016). 
% Dynamic mode decomposition: data-driven modeling of complex systems. 
% Society for Industrial and Applied Mathematics.
% [2] GAVISH, Matan; DONOHO, David L. 
% The optimal hard threshold for singular values is $4/\sqrt {3} $.
% IEEE Transactions on Information Theory, 2014, 60.8: 5040-5053.

% In what follows, if it is not stated explicitly, the pages and equation
% numbers are referred to [1]
% See Kutz et al. (2016) pages 119 - 126
dt = 1/Fs; % Time increment
N = length(xN);

% Dynamic mode decomposition
stack = ceil(N/2); % This is the our innovation
% Construct Hankel matrix
D = zeros(N - stack, stack);
for j = 1:1:stack
    D(:, j) = xN(j:N-stack +j-1);
end
D2 = D(2:end, :);
D1 = D(1:end-1,:);

[U, S, V] = svd(D1,'econ');

% Eq. (1.20)
Atilde = U(:,1:rr)' * D2 * V(:,1:rr) * inv(S(1:rr, 1:rr)); % size of Atilde is rr x rr
% Eq. (1.22)
[W, lamda] = eig(Atilde);
% Eq. (1.23)
% Each column of Phi corresponds to one mode
Phi = D2 * V(:,1:rr) * inv(S(1:rr, 1:rr) ) * W; % size of Phi is (N-stack-1) x rr

% Before Eq. (1.24)
omega = log( diag(lamda) ) * Fs; % size of omega is rr x 1
% Eq. (1.25)
% Each element of b corresponds to one mode
b = Phi\ D1(:,1) ; % size of b is rr x 1. D1(:,1) = Phi * b.

% Page 120, after algorithm 8.1
dmdAmp= abs(b) * 2 / sqrt(stack); % Amplitude spectrum

% Page 123, subsection 8.1.2
xi = size(D1, 2);
ampT = zeros(rr, xi);
for i = 1:1:xi
    bxi = Phi \ D1(:, i);
    ampT(:,i) = abs(bxi) * 2 / sqrt(stack);
end
dmdFFTAmp = mean(ampT, 2);

% Frequency (including frequency for conjugate pairs)
% Postive frequency and negative frequency
dmdTemp = imag(omega) / (2*pi); 

dampF = real(omega);
% Convert all frequencies into positive frequency
dmdF = abs(dmdTemp);

% Select out the components that have the largest amplitude
[~, indM] = max(dmdAmp);
% If dmdF(indM) >0, rowM should have two elements
[rowM, ~] = find( round(dmdF, 12, 'significant') == round(dmdF(indM), 12, 'significant') );

%% DMD reconstruction with the components that have the largest amplitude
% Eq. (1.24)
t = (0:N-1)*dt; % time vector
time_dynamics = zeros( length(rowM), N);
for i = 1:1:N
    time_dynamics(:,i) = b(rowM) .*  exp(omega(rowM) * t(i));
end
% xdmd has a size of (N-stack-1) x N
xdmd =  Phi(:,rowM) * time_dynamics;


% Determine rr by optimal hard thresholding
% Denoising
% Optimal hard thresholding method developed by GAVISH and DONOHO (2014)
beta = min( [size(D1,1) / size(D1, 2), size(D1,2) / size(D1,1)] );
tau = optimal_SVHT_coef(beta, 0) * median(diag(S));
[row,~] = find( diag(S) > tau);
rOpt = max(row);

% Eq. (1.20)
AtildeR = U(:,1:rOpt )' * D2 * V(:,1:rOpt) * inv(S(1:rOpt, 1:rOpt)); % size of Atilde is rr x rr
% Eq. (1.22)
[WR, lamdaR] = eig(AtildeR);

% Before Eq. (1.24)
omegaR = log( diag(lamdaR) ) * Fs; % size of omega is rOpt x 1

dmdTempR = imag(omegaR) / (2*pi); 
% Convert all frequencies into positive frequency
dmdFR = abs(dmdTempR);

PhiR = D2 * V(:,1:rOpt) * inv(S(1:rOpt, 1:rOpt) ) * WR; % size of Phi is (N-stack-1) x rOpt

% Eq. (1.25)
% Each element of b corresponds to one mode
bR = PhiR\ D1(:,1) ; % size of b is rOpt x 1. D1(:,1) = PhiR * b.

% Page 120, after algorithm 8.1
dmdAmpR= abs(bR) * 2 / sqrt(stack); % Amplitude spectrum



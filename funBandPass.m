function yfilt = funBandPass(Fs, y, fLow, fHigh)
% yfilt = funBandPass(Fs, y, fLow, fHigh)
%-------------------------------------
% INPUT:
% Fs - frequency in Hz
% y - scalar time series, column or row vector
% fLow - low cut-off frequency in Hz (scalar)
% fHigh - high cut-off frequency in Hz (scalar)
% OUTPUT:
% yfilt - band-pass filtered scalar time series in column vector
%--------------------------------------
% Example:
% Fs = 200;
% N = 20000;
% t = (0:1:N-1) /Fs;
% t = t(:);
% T = [1:1:10];
% f = 1./T;
% A = [1:-0.1:0.1];
% x = zeros(length(t), length(A));
% rng(1, 'twister');
% phi = rand(length(A),1) * 2 * pi;
% for i = 1:1:length(A)
%     x(:,i) = A(i) * sin( 2 * pi /T(i) * t + phi(i));
% end
% y = sum(x, 2);
% figure; 
% plot(t, y, 'b-'); 
% xlabel('Time (s)');
% ylabel('Signal');
% fLow = 0.25-1e-6; fHigh = 0.25 + 1e-6;
% yfilt = funBandPass(Fs, y, fLow, fHigh);
% figure; 
% plot(x(:,1), 'b-'); hold on; plot(yfilt, 'r-');
%--------------------------------------------
% Coded by Hongtao, 27-12-2019
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%--------------------------------------------

if (size(y,1) >1 && size(y,2) > 1)
    error('Time series must be univariate');
end


y = y(:);
N  = length(y);

% Note: divided by N. This is attributed to the definition of
% fast fourier transform that is implemented in MATLAB
yf = fft(y) / N;

% Wavenumber is stored in this way in MATLAB
% k = 2 * pi / L * [0 : (n/2 -1), (-n/2):1:-1]

% see pages 227 and 234
% in Kutz, J. N. (2013). 
% Data-driven modeling & scientific computation: methods for complex systems & big data. 
% Oxford University Press.

% see page 24 in 
% Trefethen, L. N. (2000). Spectral methods in MATLAB (Vol. 10). SIAM.

% Positive frequency and negative frequency
fHz = Fs/N *  [0 : 1 : (floor(N/2) + mod(N,2) -1),  -floor(N/2) :1:-1] ;

% Perform bandpass filtering
[~, col] = find( abs(fHz) >= fLow & abs(fHz) <= fHigh);

zf = zeros(size(yf));
zf(col) = yf(col);

% Note: multiplied with N. This is attributed to the definition of inverse fast
% fourier transform that is implemented in MATLAB
yfilt = ifft(zf) * N; % Convert back to time series

yfilt  = real( yfilt); % Convert back to real time series






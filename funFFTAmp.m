function [Fx, aFT] = funFFTAmp(a, Fs)
% function [Fx, aFT] = funFFTAmp(a, Fs)
% Input
% a is time series data
% Fs is sampling frequency
% Output
% aFT is amplitude spectrum
% Fx is frequency in Hz (row vector)
% -----------------------------------------
% Coded by Hongtao, 21-11-2018
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%------------------------------------------

aF = abs(fft(a)) / length(a); % Must be divided by length(a), see note: 17/12/2018
aFT = aF(1:1:floor( length(aF) / 2) + 1 );
aFT(2:end) = 2 * aFT(2:end); % Must multiply 2 in order to obtain correct amplitude
Fx = (0 : 1 : length(aFT)-1 ) / length(aFT) * Fs / 2;


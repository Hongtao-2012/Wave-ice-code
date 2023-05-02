function [fkern, fkern1, bf, af] = funBandPassDesign(fLow, fHigh, firOrder, iirOrder, Fs, trans, plotFlag)
% [fkern, fkern1, bf, af] = funBandPassDesign(fLow, fHigh, firOrder, iirOrder, Fs, trans, plotFlag)
% INPUT:
% fLow - low cut-off frequency
% fHigh - high cut-off frequency
% firOrder - finite impulse response filter order (for fir1.m and firls.m)
% iirOrder - infinite impulse response filter order (for butter.m)
% Fs - sampling frequency in Hz
% trans - width of transition zone (set between 0.1 and 0.4)
% plotFlag - 1 for plotting and 0 for no plotting
% OUTPUT: 
% fkern - filter coefficient for firls.m
% fkern1 - filter coefficient for fir1.m
% bf and af - filter coefficients for butter.m
% ------------------------------------------
% Coded by Hongtao, 29-12-2019
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
% -------------------------------------------
% References:
% https://www.udemy.com/course/signal-processing/learn/lecture/11864684#questions
% Section 5: Filtering of course - Signal processing problems, solved in
% MATLAB and in Python

% firls
band = [fLow, fHigh];
shape = [0 0 1 1 0 0];
lowtrans = trans;
uptrans = trans;
bandFreq = [0 (1-lowtrans) * band(1)  band(1)  band(2)  (1+uptrans) * band(2)  Fs/2] / (Fs/2);
fkern = firls(firOrder, bandFreq, shape);

% Power spectrum of fkern
n = length(fkern);
ff = (0:1:floor(n /2)-1) * Fs/n;
fpow = abs( fft(fkern)).^2;
fpow = fpow(1:1:length(ff));

% fir1
bandFreq1 = [0 band(1) band(1)  band(2)  band(2)  Fs/2] / (Fs/2);
fkern1 = fir1(firOrder, band/(Fs/2) );

% Power spectrum of fkern1
n1 = length(fkern1);
ff1 = (0:1:floor(n1 /2)-1) * Fs/n1;
fpow1 = abs( fft(fkern1)).^2;
fpow1 = fpow1(1:1:length(ff1));

% butterworth iir filter
[bf, af] = butter(iirOrder, band/(Fs/2), 'bandpass');
imps = [zeros(500,1); 1; zeros(500,1)];

% Power spectrum of iir filter coefficients
impsFilt = filter(bf, af, imps);
fImp = [0:1:floor( length(imps)/2) - 1] * Fs/(length(imps));
zf = abs( fft(impsFilt) ).^2;
zf = zf(1:1:length(fImp));

if plotFlag == 1

figure('Name',  'function: funBandPassDesign - Fig. 1', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
subplot(2, 2, [1 2]); 
plot(fkern);
xlabel('Time');
ylabel('firls filter coefficients');

subplot(2, 2, 3); 
plot(ff, fpow, 'b-'); 
hold on; plot(bandFreq * Fs/2, shape, 'r-.');
xlabel('Frequency [Hz]');
ylabel('Power spectrum - firls filter coefficients');

subplot(2, 2, 4); 
plot(ff, 10 * log10(fpow), 'b-'); 
hold on; plot(bandFreq * Fs/2, 20 * log10(shape + 1e-7), 'r-.');
xlabel('Frequency [Hz]');
ylabel('Power spectrum - firls filter coefficients [dB]');
legend({'Designed FIRls bandpass filter'; 'Target bandpass filter'});

figure('Name',  'function: funBandPassDesign - Fig. 2', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
subplot(2, 2, [1 2]); 
plot(fkern1);
xlabel('Time');
ylabel('fir1 filter coefficients');

subplot(2, 2, 3); 
plot(ff1, fpow1, 'b-'); 
hold on; plot(bandFreq1 * Fs/2, shape, 'r-.');
xlabel('Frequency [Hz]');
ylabel('Power spectrum - fir1 filter coefficients');

subplot(2, 2, 4); 
plot(ff1, 10 * log10(fpow1), 'b-'); 
hold on; plot(bandFreq1 * Fs/2, 20 * log10(shape + 1e-7), 'r-.');
legend({'Designed FIR1 bandpass filter'; 'Target bandpass filter'});
xlabel('Frequency [Hz]');
ylabel('Power spectrum - fir1 filter coefficients [dB]');

figure('Name',  'function: funBandPassDesign - Fig. 3', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
subplot(1, 2, 1);
plot(fImp, zf, 'b-'); hold on; plot(bandFreq * Fs/2, shape, 'r-.');
xlabel('Frequency [Hz]');
ylabel('Power spectrum - IIR filter coefficients');

subplot(1, 2, 2);
plot(fImp, 10 * log10(zf), 'b-'); hold on; plot(bandFreq * Fs/2, 20 * log10(shape + 1e-20), 'r-.');
xlabel('Frequency [Hz]');
ylabel('Power spectrum - IIR filter coefficients [dB]');
legend({'Designed IIR bandpass filter'; 'Target bandpass filter'});
end
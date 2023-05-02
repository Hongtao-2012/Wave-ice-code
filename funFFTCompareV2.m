function [] = funFFTCompareV2(t, yInfo, Fs)
% [] = funFFTCompareV1(t, yInfo, Fs)
% INPUT:
% t - time vector
% yInfo - a structure
% Fs - sampling frequency
% OUTPUT:
% None
%
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.

y = yInfo.y;
ybp = yInfo.ybp;
yfir1 = yInfo.yfir1;
yTik = yInfo.yTik;
ysm = yInfo.ysm;

figure('Name', 'function: funFFTCompareV1 - Fig. 1', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
plot(t, y, 'b-');  hold on; plot(t, ybp, 'k-.'); 
hold on; plot(t, yfir1, 'g-.'); 
hold on; plot(t, yTik, 'r-.'); hold on; plot(t, ysm, 'm-.');
xlabel('Time [s]');
legend({'Data'; 'bp'; 'fir1'; 'Tikhonov'; 'Smooth'});

[aF, z] = funFFTAmp(y, Fs);
[~, zbp] = funFFTAmp(ybp, Fs);

[~, zfir1] = funFFTAmp(yfir1, Fs);
[~, zTik] = funFFTAmp(yTik, Fs);
[~, zsm] = funFFTAmp(ysm, Fs);

figure('Name', 'function: funFFTCompareV1 - Fig. 2', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
subplot(1, 2, 1);
plot(aF, z, 'b-'); hold on; plot(aF, zbp, 'k-.');
hold on; plot(aF, zfir1, 'g-.');
hold on; plot(aF, zTik, 'r-.'); hold on; plot(aF, zsm, 'm-.');
ylabel('Amplitude spectrum');
xlabel('Frequency (Hz)');
legend({'Data'; 'bp'; 'fir1'; 'Tikhonov'; 'Smooth'});

subplot(1, 2, 2);
semilogx(aF, 20* log10(z), 'b-'); hold on; semilogx(aF, 20* log10(zbp), 'k-.');
hold on; semilogx( aF, 20* log10(zfir1), 'g-.');
hold on; semilogx( aF, 20* log10(zTik), 'r-.'); hold on; semilogx( aF, 20* log10(zsm), 'm-.');
ylabel('Amplitude spectrum (dB)');
xlabel('Frequency (Hz)');
legend({'Data'; 'bp'; 'fir1'; 'Tikhonov'; 'Smooth'});

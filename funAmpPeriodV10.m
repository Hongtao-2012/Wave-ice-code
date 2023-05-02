function [AmpT] = funAmpPeriodV10(Tw, A, t, x, Fs, Fq, tsec, rr, lb, ub, isPlot)
% [AmpT] = funAmpPeriodV4(Tw, A, t, x, Fs, Fq, tsec, lb, ub, isPlot)
%---------------------------------------------------------------
% INPUT:
% Tw - prescribed (target) wave period in seconds
% A -  peak height (determine roughly based on observations made for
% steady wave signal, should be smaller than largest peak height in steady wave signal)
% t - time vector for whole univariate time series x
% x - univariate time series
% Fs - sampling frequency in Hz
% Fq - used upsampling frequency in Hz (best with twice or fourth times of
% Fs)
% rr - modes to be kept for Prony method (use funPronyTSVD.m to determine)
% tsec - a two element vector corresponding with beginning and the end of
% steady wave signal
% lb - 3-elements vector, lower bounds for genetic algorithm. [A, T, phi] - Amplitude, period,
% phase angle
% ub - 3-elements vector, upper bounds for genetic algorithm. [A, T, phi] - Amplitude, period,
% phase angle
% isPlot - 1 for plotting and output to command window, 0 for no plotting
% and output to command window
% OUTPUT:
% AmpT - a structure.
% Field names start with A respresent Amplitude
% Field names begin with f represent frequency in Hz
% Peak analysis of interpolated signal
% AmpT.Am
% AmpT.fm
% Prony analysis of original signal
% AmpT.Ap
% AmpT.fp
% fft of original signal
% AmpT.Af
% AmpT.ff
% fft of interpolated signal
% AmpT.Afm
% AmpT.ffm
% DMD of original signal
% AmpT.Ad
% AmpT.fd
% Hilbert transform of original signal
% AmpT.Ahb
% Frequency is identical with fft of original signal
% Hilbert transform of interpolated signal
% AmpT.Ahbm
% Frequency is identical with fft of interpolated signal
% ----------------------------------------------------------------
% Coded by Hongtao, 30-12-2019
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
% ----------------------------------------------------------------
% References:
% [1] SREE, Dharma KK; LAW, Adrian Wing-Keung;
% SHEN, Hayley H. An experimental study on gravity waves
% through a floating viscoelastic cover.
% Cold Regions Science and Technology, 2018, 155: 289-299.
% [2] https://se.mathworks.com/help/matlab/ref/makima.html, Accessed 03-01-2020
% [3] Kutz, J. N. (2013).
% Data-driven modeling & scientific computation: methods for complex systems & big data.
% Oxford University Press.
% [4] Trefethen, L. N. (2000). Spectral methods in MATLAB (Vol. 10). SIAM.

dt = 1/Fs;
indx = ( t >=tsec(1) & t <= tsec(end));
xN = x(indx);
tsel = t(indx);
tq = tsel(1) : 1/Fq: tsel(end);
% Interpolation
% see [1] and [2]
ym = interp1(tsel, xN, tq, 'makima');
ym = ym(:);

% Peak analysis
[pkp, locp] = findpeaks(ym, 'MinPeakDistance',  0.8 * Tw * Fq, 'MinPeakHeight', A * 0.5, 'Threshold', A * 1e-12);
[pkn, locn] = findpeaks(-ym, 'MinPeakDistance', 0.8 *  Tw * Fq, 'MinPeakHeight', A * 0.5, 'Threshold', A * 1e-12);

if length(pkp) > length(pkn)
    pkp(end) = [];
    locp(end) = [];
end

if length(pkn) > length(pkp)
    pkn(end) = [];
    locn(end) = [];
end

Am  = mean( pkp + pkn) /2;
Tm  = 2 * mean( abs( tq(locp) - tq(locn) ) );
fm = 1/Tm;

% Prony method
[alphac, fpc, Amp_pc, thetac] = funPronyExtract(xN,rr, dt);
[Ap, indM] = max(Amp_pc);
fp = fpc(indM);

% Reconstruct signal using component that has maximum amplitude, which is
% deterimed by using Prony method
tk = tsel - tsel(1);
tk = tk(:);
xprec =  real( Amp_pc(indM) * exp(1i * thetac(indM,1) ) * exp( -alphac(indM,1) * tk + 1i * 2 * pi * fpc(indM,1) * tk) );


% FFT
[Fx, aFT] = funFFTAmp(xN, Fs);

[maFT, ind] = max(aFT);
ff = Fx(ind);
Af = maFT;

% Reconstruct signal by using component that is largest in amplitude
% spectrum obtained by FFT
yf = fft(xN)/length(xN);
xf2 = zeros( size(yf));

% Wavenumber is stored in this way in MATLAB
% k = 2 * pi / L * [0 : (n/2 -1), (-n/2):1:-1]

% see pages 227 and 234
% in Kutz, J. N. (2013).
% Data-driven modeling & scientific computation: methods for complex systems & big data.
% Oxford University Press.

% see page 24 in
% Trefethen, L. N. (2000). Spectral methods in MATLAB (Vol. 10). SIAM.
[pks, locs] = findpeaks(aFT,'SortStr','descend');
indx = locs(1);
% Fourier coefficient corresponding with positive frequency
xf2(indx) = yf(indx);
% Fourier coefficient corresponding with negative frequency
xf2(end-indx+2) = yf(end-indx+2);
xfrec = real(ifft(xf2)) * length(xN);

% FFT of interpolated signal
[Fxm, aFTm] = funFFTAmp(ym, Fq);
[maFT, ind] = max(aFTm);
ffm = Fxm(ind);
Afm = maFT;

% Dynamic mode decomposition
[b, dmdF, dmdAmp, dmdFFTAmp, dmdFR,  dmdAmpR,  xdmd, dampF] = funDMDRT(xN,rr, Fs);
[Ad, indMax] = max(dmdAmp);
fd = dmdF( indMax);

[Adf, indMax] = max(dmdFFTAmp);
fdf = dmdF( indMax);

[AdR, indMax] = max(dmdAmpR);
fdR = dmdFR( indMax);


% Genetic algorithm
isCheck = 0;
[Ag, Tg, alphag] = funGeneProfileV1(tsel, xN, lb, ub, isCheck);
fg = 1/Tg;
% Reconstruct signal by using outputs from Genetic algorithm
xgrec = Ag * cos( 2* pi /Tg * tsel+ alphag);

% Hilbert transform
xhb = hilbert(xN);
Ahb = mean( abs(xhb) );

% Hilbert transform of interpolated signal
xhbm = hilbert(ym);
Ahbm = mean( abs(xhbm) );

% Store prescribed (target) wave frequency
AmpT.Tw = Tw;

% Store frequency and amplitude for components with largest amplitude
% Peak analysis of interpolated signal
AmpT.Am = Am;
AmpT.fm = fm;

% Prony method
AmpT.Ap = Ap;
AmpT.fp = fp;

% FFT of orignal signal
AmpT.Af = Af;
AmpT.ff = ff;

% FFT of interpolated signal
AmpT.Afm = Afm;
AmpT.ffm = ffm;

% Dynamic mode decomposition (rr determined by TSVD)
AmpT.Ad = Ad;
AmpT.fd = fd;

% Dynamic mode decomposition (averged amplitude)
AmpT.Adf = Adf;
AmpT.fdf = fdf;

% Dynamic mode decomposition (rOpt determined by algorithm proposed by Gavish and Donoho (2014))
AmpT.AdR = AdR;
AmpT.fdR = fdR;

% Genetic algorithm
AmpT.Ag = Ag;
AmpT.fg = fg;

% Hilbert transform of original signal
AmpT.Ahb = Ahb;

% Hilbert transform of interpolated signal
AmpT.Ahbm = Ahbm;


% Post-processing
if isPlot == 1
    figure('Name', 'function: funAmpPeriodV4 - Fig. 1', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    plot(t, x, 'k-');
    hold on; plot(tsel, xN, 'r-.');
    hold on; plot(tq, ym, 'b--');
    hold on; plot(tq(locp), pkp, 'ro'); hold on; plot(tq(locn), -pkn, 'bo');
    xlabel('Time [s]');
    ylabel('x');
    legend({'Data'; 'Selected data'; 'Interpolated data'; 'Positive peak'; 'Negative peak'});
    title('Find peaks');
    
    figure('Name',  'function: funAmpPeriodV4 - Fig. 2', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    subplot(3,1,1); plot(tsel, xN); hold on; plot(tsel, real(xhb), 'r-.');
    xlabel('Time [s]');
    ylabel('x');
    legend({'Data'; 'Hilbert'});
    title('Hilbert transform');
    
    subplot(3, 1, 2); plot(tsel, abs(xhb), 'r-o');
    xlabel('Time [s]');
    ylabel('Magnitude');
    subplot(3, 1, 3); plot(tsel,angle(xhb), 'r-o');
    xlabel('Time [s]');
    ylabel('Phase angle');
    
    figure('Name',  'function: funAmpPeriodV4 - Fig. 3', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    subplot(3,1,1); plot(tq, ym); hold on; plot(tq, real(xhbm), 'r-.');
    legend({'Data'; 'Hilbert'});
    xlabel('Time [s]');
    ylabel('x');
    title('Interp - Hilbert transform');
    
    subplot(3, 1, 2); plot(tq, abs(xhbm), 'r-o');
    xlabel('Time [s]');
    ylabel('Magnitude');
    
    subplot(3, 1, 3); plot(tq,angle(xhbm), 'r-o');
    xlabel('Time [s]');
    ylabel('Phase angle');
    
    figure('Name',  'function: funAmpPeriodV4 - Fig. 4', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    plot(tsel, xN, 'b-');
    hold on; plot(tsel, xfrec, 'k-.');
    hold on; plot(tsel, xgrec, 'r-.');
    hold on; plot(tsel, xprec, 'g-.'); hold on; plot(tsel(1:size(xdmd,2)), real(xdmd(1,:)), 'm-.');
    xlabel('Time [s]');
    title('Reconstruction with components that have the largest amplitude');
    legend({'Data'; 'fft'; 'Genetic'; 'Prony'; 'DMD'});
    
    figure('Name',  'function: funAmpPeriodV4 - Fig. 5', 'color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    subplot(1,2,1);
    plot(Fx, aFT, 'b-');
    hold on; plot(Fxm, aFTm, 'r-.');
    hold on; plot(dmdF, dmdAmp, 'ro');
    hold on; plot(fm, Am, 'k>');
    hold on; plot(fpc, Amp_pc, 'm+');
    hold on; plot(fdf, Adf, 'kd');
    hold on; plot(dmdFR, dmdAmpR, 'rs');
    hold on; plot(fg, Ag, 'gs');
    hold on; plot(ff, Ahb, 'b^');
    hold on; plot(ffm, Ahbm, 'ko');
    xlabel('Frequency [Hz]');
    title('Amplitude spectrum');
    legend({'fft'; 'Interp-fft'; 'DMD'; 'Peak'; 'Prony'; 'DMD (averaged)'; 'DMD (rOpt)'; 'Genetic'; 'Hilbert'; 'Interp - Hilbert'});
    
    subplot(1,2,2);
    plot(Fx, 20*log10(aFT), 'b-');
    hold on; plot(Fxm, 20*log10(aFTm), 'r-.');
    hold on; plot(dmdF, 20*log10(dmdAmp), 'ro');
    hold on; plot(fm, 20*log10(Am), 'k>');
    hold on; plot(fpc, 20*log10(Amp_pc), 'm+');
    hold on; plot(fdf, 20 * log10(Adf), 'kd');
    hold on; plot(dmdFR, 20 * log10(dmdAmpR), 'rs');
    hold on; plot(fg, 20*log10(Ag), 'gs');
    hold on; plot(ff, 20*log10(Ahb), 'b^');
    hold on; plot(ffm, 20*log10(Ahbm), 'ko');
    xlabel('Frequency [Hz]');
    title('Amplitude spectrum [dB]');
    legend({'fft'; 'Interp-fft'; 'DMD'; 'Peak'; 'Prony'; 'DMD (averaged)'; 'DMD (rOpt)'; 'Genetic'; 'Hilbert'; 'Interp - Hilbert'});
    
    fprintf('\n---------Function: funAmpPeriodV10.m ------------------- \n');
    fprintf('---------- Maximum amplitude and its corresponding frequency ----------- \n');
    fprintf('Peak analysis of interpolated signal: Amp: %f  Freq:  %f [Hz]  \n\n', Am, fm);
    fprintf('Prony analysis of original signal: Amp: %f  Freq:  %f [Hz]  \n\n', Ap, fp);
    fprintf('fft of original signal: Amp: %f  Freq:  %f [Hz]  \n\n', Af, ff);
    fprintf('fft of interpolated signal: Amp: %f  Freq:  %f [Hz]  \n\n', Afm, ffm);
    fprintf('DMD of original signal: Amp: %f  Freq:  %f [Hz]  \n\n', Ad, fd);
    fprintf('DMD (averaged) of original signal: Amp: %f  Freq:  %f [Hz]  \n\n', Adf, fdf);
    fprintf('DMD (rOpt) of original signal: Amp: %f  Freq:  %f [Hz]  \n\n', AdR, fdR);
    fprintf('Genetic algorithm for original signal: Amp: %f  Freq:  %f [Hz]  \n\n', Ag, fg);
    fprintf('Hilbert transform of original signal: Amp: %f  \n\n', Ahb);
    fprintf('Hilbert transform of interpolated signal: Amp: %f  \n\n', Ahbm);
    fprintf('-----------------------End of Function: funAmpPeriodV10.m ------------------- \n');
end
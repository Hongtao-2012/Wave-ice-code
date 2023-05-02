clear all; close all; clc;

set(0, 'defaultFigureColor', 'w');
set(0, 'defaultLineLineWidth', 1.5);
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');
set(0, 'defaultAxesFontSize', 24);
set(0, 'defaultTextFontSize', 24);
set(0, 'defaultTextInterpreter', 'latex');
set(0,'DefaultLineMarkerSize',7);
set(0, 'defaultAxesTickLabelInterpreter','latex');  
set(0,'defaultAxesLineWidth', 1.5);

%% NOTE:
% If you use any scripts that accompany this example, please cite following
% paper

%% References:
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.


addpath('./utils');

Fs = 200; % Sampling frequency
N = 5601; % Number of data points
t = (0:1:N-1) /Fs; % Time vector
t = t(:);
dt = 1/Fs; % Constant time increment

T = [1 1.1 1.2 4 1.5 1.6 1.8 2 2.2 8]; % Period in seconds. Target wave period is 4 seconds
f = 1./T; % Frequency

rng(50,'twister');
numRand = 0.5 * rand(1,9);
A = [ numRand(1:3) 18 numRand(4:9)]; % Amplitude
x = zeros(length(t), length(A)); 

rng(1, 'twister');
phi = 0.01 * rand(length(A),1) * 2 * pi; % Phase angle

rng(10, 'twister');
dampVec = rand(10, 1) * 10^(-3); % grow/damp rate

rng(20, 'twister');
randsign = 2*(rand(10,1)>0.5) - 1; % sign of grow/damp rate

%Generate synthetic signal
for i = 1:1:length(A)
   switch i
       case 4
        dampFact = randsign(i) * 1e-3;
       otherwise
        dampFact = randsign(i) * dampVec(i) ;
    end
      x(:,i) = A(i) * exp(-dampFact * t) .* sin( 2 * pi /T(i) * t + phi(i));
end
xS = sum(x, 2);

% Add 1% of white noise
xN = xS + 0.01 * std(xS) * randn(size(xS));


figure('Name', 'Select cyclic part with zero periodic boundary');
plot(xN);
indCyclic = [405 5202];
hold on;
plot( [indCyclic(1) indCyclic(1)], get(gca, 'ylim'), 'k-.');
hold on;
plot( [indCyclic(2) indCyclic(2)], get(gca, 'ylim'), 'k-.');
xlabel('Index');
ylabel('y');
title('Select cyclic part with zero periodic boundary');

xCyc = xN(indCyclic(1) : indCyclic(2) );
% Isolate f = 1/T(4) = 0.25 Hz compoent from noise and other components
% by using bp method (frequency domain filtering)
fLow = 0.15; fHigh = 0.4;
ybp= funBandPass(Fs, xCyc , fLow, fHigh);

% FIR and IIR filter design
fLow2 = 0.15; % lower bound in Hz
fHigh2 = 0.4; % upper bound in Hz
firOrder =  800;
iirOrder = 3;
trans = 0.08; % transition zone width for fir1 filter design
plotFlag = 0;
[fkern, fkern1, bf, af] = funBandPassDesign(fLow2, fHigh2, firOrder, iirOrder, Fs, trans, plotFlag);

yfir1= filtfilt( fkern1, 1, xCyc );
yiir= filtfilt(bf, af, xCyc );

% Uncomment following to determine the regularization parameter for Tikhonov
% regularization
% alphaVec = logspace(-9,-4, 400);
%%%  select 500 to 1000 data points 
% alpha= funLCurve(xCyc(500:1050), dt, alphaVec, 'xCyc');

%%% Determined regularization parameter by L-curve method
% alpha = 2.7148e-06;
% alpha = 2.8751e-06;
alpha = 4.1836e-06;

% Tikhonov regularization to denoise signal
np = length(xCyc );
xstack = [xCyc(1:end-2); zeros(np-2,1)];
maxiter = 9000;
yTik = lsqr(@(x, trans_flag)funTikhonovDenoise(x, trans_flag, dt, alpha), xstack, [], maxiter);

ysm = smooth(xCyc, 0.02,  'rloess');

% Point-wise comparison for filtering results obtained by using fir1 and bp
% methods (check this value)
[nrmse] = funNRMSE(yfir1, ybp);
nrmse = nrmse * 100;

% Compare filtering and smoothing results in time and frequency domain
yInfo.nrmse = nrmse;
yInfo.y = xCyc;
yInfo.ybp = ybp;
yInfo.yfir1 = yfir1;
yInfo.yiir = yiir;
yInfo.yTik = yTik;
yInfo.ysm = ysm;
funFFTCompareV2( t(indCyclic(1) : indCyclic(2) ), yInfo, Fs);

lb = [12 2 0]; % [amplitude, period in seconds, phase angle in radians] lower bound used by ga.m to fit a sinusoidal wave
ub = [22 5 2*pi-1e-6]; % [amplitude, period in seconds, phase angle in radians] lower bound used by ga.m to fit a sinusoidal wave


figure('Name', 'Select steady part with zero periodic boundary');
plot(xN);
indSteady = [1195 3609 ];

hold on;
plot( [indCyclic(1) indCyclic(1)], get(gca, 'ylim'), 'k-.');
hold on;
plot( [indCyclic(2) indCyclic(2)], get(gca, 'ylim'), 'k-.');

hold on;
plot( [indSteady(1) indSteady(1)], get(gca, 'ylim'), 'r-.');
hold on;
plot( [indSteady(2) indSteady(2)], get(gca, 'ylim'), 'r-.');

xlabel('Index');
ylabel('y');
title('Select steady part with zero periodic boundary');

% Select number of components to be retained for Prony's method
funPronyTSVD(  xN( indSteady(1) : indSteady(2) )    );
r = 2;

tSteady = t( indSteady(1) : indSteady(2) ); % Steady part of the signal
% Identify wave period and wave amplitude (or oscillation amplitude)
AmpT = funAmpPeriodV10(T(4), 12, t( indCyclic(1) : indCyclic(2) ) , yTik, Fs, 2 * Fs, tSteady ,r, lb, ub, 1);

% save('yInfo_NumericalExample_DMD.mat', 'yInfo');
% save('AmpT_NumericalExample_Tik.mat', 'AmpT');



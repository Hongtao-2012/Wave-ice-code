function [tw, tb] = funSteadyWave(Tt, h, t, x, dw, db, qwTime)
%  [tw, tb] = funSteadyWave(Tt, h, t, x, dw, db, qwTime)
%  -------------------------------------------------------------
%  INPUT:
%  Tt - prescribed/target wave period in seconds
%  h - water depth
%  t - logging time vector
%  x - wave elevation 
%  dw - distance from wavemaker to the wave elevation measurement device
%  db - distance from the wave elevation measurement device to dissipative
%  beach
%  qwTime - the difference between begining of logging and the start of
%  wavemaker
%  OUTPUT:
%  tw - time for waves to reach wave elevation measurement device relative
%  to start of logging time
%  tb - time for waves to travel the round trip between wave elevation measurement device relative
%  to start of logging time
%  -------------------------------------------------------------
%  Example
% Tt = 2; % Target wave period
% h = 2.45; % Water depth
% Fs = 100; % Sampling frequency
% N = 12000; % Number of data points
% t = (0:1:N-1)/Fs; % Logging time
% f = 1/Tt; % 
% x = 2 * sin( 2 * pi * f * t);
% qwTime = 2.5;  % Wavemaker starts 2.5 seconds after logging time
% dw = 25;
% db = 44;
% [tw, tb] = funSteadyWave(Tt, h, t, x, dw, db, qwTime);
%---------------------------------------------------------------
% Coded by Hongtao, 04-01-2020
%% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%---------------------------------------------------------------
% References:
% [1] Cheng, S., Tsarau, A., Evers, K. U., & Shen, H. (2019).
%  Floe size effect on gravity wave propagation through ice covers. 
% Journal of Geophysical Research: Oceans, 124(1), 320-334.
% [2] HUANG, Guoxing; LAW, Adrian Wing-Keung; HUANG, Zhenhua. 
% Wave-induced drift of small floating objects in regular waves. 
% Ocean Engineering, 2011, 38.4: 712-718.

g = 9.81;
 [L, Lw, Lk, LkT, LkS, LkG] = funWLOpen(Tt, h);  % Calculate wavelength
 kT = 2 * pi / Lk; % Wavenumber
Cg = g * Tt/ (4 * pi) * (tanh(kT * h) + kT*h / (cosh(kT *h)^2) );  % Wave group velocity in open water

% Time for wave reaches measurement position
tw = dw / Cg + qwTime ;  

% Time for wave travels round trip between measurement position and damping
% beach
% see Refs. [1] - [2]
tb = (dw + 2 * db) / Cg + qwTime; 
tc = t;

figure('Name', 'function: funSteadyWave', 'Color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
plot(tc, x, 'b-');
hold on;
plot([tw, tw], get(gca,'ylim'), 'r-.');
hold on;
plot([tb tb], get(gca, 'ylim'), 'r-.' );
xlabel('Time [s]');
ylabel('Wave elevation');


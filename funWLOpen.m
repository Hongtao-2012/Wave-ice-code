function [L, Lw, Lk, LkT, LkS, LkG] = funWLOpen(T, h)
% [L, Lw, Lk, LkT, LkS, LkG] = funWLOpen(T, h)
%------------------------------------------
% INPUT:
% T - wave period in second
% h - water depth in meter
% OUTPUT:
% L - wavelength obtained by Bob You's formula
% Lw - wavelength obtained by iteration (wavelengh should be in this range
% 0.5 : 1e-6 : 160)
% Lk - wavelength obtaine by employing w2k algorithm
% LkT - wavelength obtaine by utilizing fminbnd matlab built-in function
% LkS - wavelength obtaine by employing fminsearch optimization matlab built-in function
% LkG - wavelength obtaine by applying ga (genetic algorithm) matlab built-in function
%-------------------------------------------------------------
% Coded by Hongtao Li for TBA4265 2018-09-25
% Updated 2018-11-21
%% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%--------------------------------------------------------------
% References:
% [1] Brodtkorb, P.A., Johannesson,
% P., Lindgren, G., Rychlik, I., Rydén, J. and Sjö, E. (2000). 
% "WAFO - a Matlab toolbox for analysis of random waves and loads",
% Proceedings of the 10th International Offshore and Polar Engineering conference,
% Seattle, Vol III, pp. 343-350.


g = 9.81;  % Gravitational acceleration
% 
% T = 10;    % Wave period
% h = 2000; % Water depth

disp('-----------------Wavelength estimated by funWLOpen.m-----------%');
%% BobYou's formula
L = funBobYou(T, h);
disp('%-----------------%');
fprintf('%s \n', 'wavelength calculated by using BobYou formula');
fprintf(['Wavelength in water depth ', num2str(h), ' m is ', '%6.3f m \n'], L);

%% Iteration
lamda0 = g * T^2 / (2 * pi);
lamda = 0.5 : 1e-6 : 160;
[Lw, minDif1] = funWaveLengthOpen(T,h, lamda);
disp('%-----------------%');
fprintf('%s \n', 'wavelength calculated by using iteration');
fprintf(['Wavelength in water depth ', num2str(h), ' m is ', '%6.3f m \n'], Lw);


%% Employing w2k algorithm
% see [1]
w = 2 * pi / T;  % Wave frequency in rad/s
th = 0;
[k1,k2,ind]=w2k(w,th,h,g);
k = sqrt( k1^2 + k2^2 );
Lk = 2 * pi / k;

disp('%-----------------%');
fprintf('%s \n', 'wavelength calculated by using w2k');
fprintf(['Wavelength in water depth ', num2str(h), ' m is ', '%6.3f m \n'], Lk);



%% Utilizing fminbnd matlab built-in function
options = optimset('MaxFunEvals', 2000, 'MaxIter', 2000, 'TolX', 1e-8);
kT = fminbnd(@(k)abs((2 * pi / T)^2 - 9.81 * k * tanh(k * h) ), 1e-5, 700,  options);
LkT = 2 * pi / kT;

disp('%-----------------%');
fprintf('%s \n', 'wavelength calculated by using fminbnd matlab built-in function');
fprintf(['Wavelength in water depth ', num2str(h), ' m is ', '%6.3f m \n'], LkT);

%% Employing fminsearch optimization matlab built-in function
kS = fminsearch(@(k)abs( (2 * pi / T)^2 - g * k * tanh(k * h) ), 0.01, options);
LkS = 2 * pi / kS;

disp('%-----------------%');
fprintf('%s \n', 'wavelength calculated by using fminsearch optimization matlab built-in function');
fprintf(['Wavelength in water depth ', num2str(h), ' m is ', '%6.3f m \n'], LkS);


%% Applying ga (genetic algorithm) matlab built-in function
goptions = optimoptions('ga'); % Retrieve default options for ga
goptions.MaxGenerations = 600;  % Set option value
goptions.FunctionTolerance = 1e-8; % Set option value

% It is critical to set lower bound of k, otherwise negative k is returned
LowB = 1e-5;  % Lower bound of k
UpB = 6000;   % Upper bound of k

% Read passing extra parameters by means of Anonymous Functions
% Nested Functions and Global Variables
[kG,~,exitflag] = ga(@(k)abs( (2 * pi /T)^2 - g * k * tanh(k * h)), 1,[], [],[],[],LowB,UpB,[],goptions);
LkG = 2 * pi / kG;

disp('%-----------------%');
fprintf('%s \n', 'wavelength calculated by using ga (genetic algorithm) matlab built-in function');
fprintf(['Wavelength in water depth ', num2str(h), ' m is ', '%6.3f m \n \n'], LkG);

fprintf('---------------End of function: funWLOpen.m-------------\n');



function [Amp, T, alpha] = funGeneProfileV1(t, z, lb, ub, isCheck)
% [Amp, T, alpha] = funGeneProfile(t, z, lp, up, isCheck)
% INPUT:
% t - time vector
% z - Vertical oscillation
% lb - 3-elements vector, lower bounds for genetic algorithm. [A, T, phi] - Amplitude, period,
% phase angle
% ub - 3-elements vector, upper bounds for genetic algorithm. [A, T, phi] - Amplitude, period,
% phase angle
% isCheck - check quality
% OUTPUT:
% Amp - amplitude
% T - Period in seconds
% alpha - Phase in radians
% z is fitted by Amp * cos( 2 * pi / T + alpha )
% ---------------------------------------
% Coded by Hongtao, 18-12-2019
%% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
% ---------------------------------------

if isCheck == true
options = optimoptions('ga','PlotFcn', @gaplotbestf, 'MaxGenerations', 9000, 'FunctionTolerance', 1e-9);
[coeff,fval,exitflag,output]  = ga(@(x0)funOptIce(x0, t(:), z(:) ),3, [],[],[],[],lb,ub,[],options);
else
options = optimoptions('ga', 'MaxGenerations', 9000,'FunctionTolerance', 1e-9 );
[coeff,fval,exitflag,output]  = ga(@(x0)funOptIce(x0, t(:), z(:) ),3, [],[],[],[],lb,ub,[],options);
end

if isCheck == true
figure;
plot(t, z(:), 'b-');
hold on;
plot(t, coeff(1)*cos(2*pi/coeff(2) * t + coeff(3)), 'r--');
end

Amp = coeff(1);
T = coeff(2);
alpha = coeff(3);
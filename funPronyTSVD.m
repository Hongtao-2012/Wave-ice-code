function [] = funPronyTSVD(yk)
% [] = funPronyTSVD(yk)
%-------------------
% INPUT:
% yk - univariate time series
% OUTPUT:
% None
%------------------------------------------
% Coded by Hongtao, 03-01-2020
%% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%------------------------------------------
% References:
% [1] HU, Sau-Lon James; YANG, Wen-Long; LI, Hua-Jun. Signal decomposition and reconstruction using complex exponential models. Mechanical Systems and Signal Processing, 2013, 40.2: 421-438.

% See subsection 3.4 Step1 and Fig. 2 in Ref. [1]
xi  = floor( length(yk)/2) -1;
eta = xi;
H0 = hankel( yk(1:xi), yk(xi : xi + eta - 1) );
[U, S, V] = svd(H0,'econ');
figure('Name', 'function: funPronyTSVD', 'Color', 'w', 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
semilogy( diag(S), 'r-o');
xlabel('Subpace dimension');
ylabel('Singular value spectrum in logarithmic scale');

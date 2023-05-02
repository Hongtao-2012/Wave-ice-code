function alphaOpt = funLCurve(xN, dt, alphaVec, figName)
% alphaOpt = funLCurve(xN, dt, alphaVec, figName)
% xN - scalar data vector
% dt - time increment
% alphaVec - test regularization parameter range for Tikhonov
% regularization
% figName - string to give names for figures
% OUTPUT:
% alphaOpt - optimal values of alpha that is determined by L-curve method
% Suggestion: select 300 - 500 data points as xN that are representative
% of the signal
%-------------------------------------------
% Coded by Hongtao, 20-12-2019
%
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%-------------------------------------------
% References:
% [1] MUELLER, Jennifer L.; SILTANEN, Samuli.
% Linear and nonlinear inverse problems with practical applications. SIAM, 2012.

% See pages 73 - 78 in Mueller and Siltanen (2012), subsection 5.4.2
xN = xN(:);
n = length(xN);

A = [ diag(ones(n-2,1)), zeros(n-2,1), zeros(n-2,1)]; % A has a size of (n-2) x n

% Second-order difference operator. L2 has a size of (n-2) x n
L2 = 1/ dt^2 * [ diag(ones(n-2,1)) + diag(-2 * ones(n-3,1), 1) + diag(ones(n-4,1), 2), [zeros(n-4,1); 1; -2], [zeros(n-3,1); 1] ];

% Construct mtilde that has a size of 2(n-2) x 1
mtilde = [xN(1:n-2); zeros(n-2,1)];

% L - curve method
a = zeros( length(alphaVec),1);
b = zeros( length(alphaVec),1);
for i = 1:1:length(alphaVec)
Atilde = [A; sqrt(alphaVec(i)) * L2];
f = Atilde \ mtilde;
a(i) = log(norm(f - xN,2));
b(i) = log(norm( L2 * f, 2) );
end

figure('Name', [figName, 'L-curve'], 'Units', 'Normalized', 'Outerposition', [0 0 1 1], 'Color', 'w');
plot(b, a, 'b-');
xlabel('$\log{ ||Af - m||}$', 'Interpreter', 'Latex');
ylabel('$\log{ ||L_2 f||}$', 'Interpreter', 'Latex');

% Manually locate the optimal regularization parameter on the L-curve plot
ins = input('Input number of significant digits : ', 's');
digit = floor(str2double(ins));
ins = input('Input abscissa of corner in L-curve figure: ', 's');
corner = str2double(ins);

[ind,~] = find( round(b,digit,'significant') == corner);
alphaOpt = alphaVec(ind);

% Use optimal alpha found by L-curve method to inversely determine f
Atilde = [A; sqrt( alphaOpt ) * L2];
f = Atilde \ mtilde;

% Post-processing
figure('Name', [figName], 'Units', 'Normalized', 'Outerposition', [0 0 1 1], 'Color', 'w');
plot( xN, 'k-'); hold on; plot( f, 'r-.');
legend({'Data'; 'Cleaned data'});
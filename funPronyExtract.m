function [alpha_p, f_p, Amp_p, theta_p] = funPronyExtract(yk, r, dt)
%function [alphac, fc, Ampc, thetac] = funPronyExtract(yk, r, dt)
% Input:
% yk - scalar time series
% r - number of subspaces to be kept. Adjacent subspace pairs for
% oscillation components and one subspace for decaying mode without
% oscillation. This number is obtained by checking L-curve of singular
% values plotted in the logrithmic scale (i.e. trungated singular value
% decomposition)
% dt - sampling time interval
% Output:
% alphac - Damping factor column vector.
% f_p - Frequency column vector. Unit: [Hz]
% Amp_p- Amplitude column vector.
% thetac - Phase column vector. Unit: [radians]
% Recontruction:
% The original signal can be reconstructed as:
% yr = zeros(size(yk));
% for i = 1:1:length(f_p)
%     yr = yr + real( Amp_p(i,1) * exp(1i * theta_p(i,1) ) * exp( -alpha_p(i,1) * t + 1i * 2 * pi * f_p(i,1) * t) );
% end
% The t used in reconstruction must start from zero
% ------------------------------------------
% Coded by Hongtao, 12-11-2019
%
% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
%-------------------------------------------
% References:
% [1] HU, Sau-Lon James; YANG, Wen-Long; LI, Hua-Jun. Signal decomposition and reconstruction using complex exponential models. Mechanical Systems and Signal Processing, 2013, 40.2: 421-438.


yk = yk(:);
N = size(yk,1);

% Eq. 26
% See discussions of choosing xi and eta 
% directly after Eq. 34
% and in Step 1 in subsection 3.4 
xi  = floor( length(yk)/2) -1;
eta = xi;
H0 = hankel( yk(1:xi), yk(xi : xi + eta - 1) );
H1 = hankel( yk(2 : xi+1), yk(xi+1 : xi + eta) );

% Eq. 34
[U, S, V] = svd(H0,'econ');

% Sdiag now is a vector
Sdiag = 1./sqrt( diag( S(1:r, 1:r) ) );
% Sdiag now is a matrix
Sdiag = diag(Sdiag);
% Eq. 36
A = Sdiag * U(:, 1:r)'  * H1 * V(:,1:r) * Sdiag;

% After Eq. (36)
[V1, D1] = eig(A);
z = diag(D1);

% Between Eqs. 2 and 4
% After Eq. (36)
alphac = real( log( z) ) / dt;
fc = imag( log( z) ) / (2 * pi * dt);

% Eq. 13
zz = zeros(N, r);
for i = 1:1:r
    zz(:,i) = z(i) .^( (0:1:N-1)');
end

% Eq. 13
% After Eq. 2
% r_n = A_n * exp(1i * theta_n);
gamma = zz \ yk;
Ampc = abs( gamma ) ;
%% Important:
% thetac = atan( imag(gamma) /real(gamma) ); When reconstructing the
% original signal, phase difference is found between recontructed signal
% and the original signal
thetac = atan2( imag(gamma), real(gamma) ); % Use atan2 to obtain Four-Quadrant Inverse Tangent 

%% Based on the non-negative sign requirement of the frequency to 
%% extract physical meaningful results
f_p = fc(fc >= 0); % Frequency
alpha_p = -alphac( fc >= 0); % Damping factor
[c, ia, ib] = unique( round( abs(fc), 12, 'significant'), 'rows', 'stable');

Amp_p = zeros(length(ia),1); % Amplitude
if length(ia) == 1
    Amp_p = sum(Ampc);
end

for i = 1:1:length(ia)-1
    Amp_p(i) = sum( Ampc( ia(i) : ia(i+1)-1) );
    if i== length(ia)-1
        Amp_p(i+1)  = sum( Ampc( ia(i+1)  : end) );
    end
end

theta_p = thetac( fc >= 0); % Phase angle


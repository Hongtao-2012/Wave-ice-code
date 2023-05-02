function E = funOptIce(x0, t, z )
% E = funOptIce(x0, t, z )
% INPUT:
% x0 - Amplitude
% t - time vector
% z - measured signal
% OUTPUT:
% E - L2 - norm difference (not normalized by length 
% of data)
% ------------------------------------------
% Coded by Hongtao, 19-02-2019
%% Please cite this paper
% 
% LI, Hongtao; GEDIKLI, Ersegun Deniz; LUBBAD, Raed.
% Systematic investigation of data analysis methods in wave-ice interaction problemSystematic investigation of data analysis methods in wave-ice interaction problem.
% In: Proceedings of the 25th IAHR International Symposium on Ice. 
% Trondheim, Norway, June 14-18, 2020. International Association for Hydro-Environment Engineering and Research (IAHR), 2020.
% ------------------------------------------
f = x0(1) * cos(2 * pi / x0(2) * t + x0(3));
E = sqrt( sum( (f- z).^2) );
end


function [ M1, M2 ] = bispec_precomp(d)
% For a given signal length, precompute the matrix M for generating
% bispectrum
% Input: 
%        d: signal length
% Output:
%        M1 and M2: Matrix that multiply the Fourier coefficients together.
%        M1 and M2 are sparse matrices, each row only has 3 non-zeros.
%        The idea is for the Fourier amplitude, 
%        log(|B(k1, k2)|) = log(|y(k1)|) + log(|y(k2)|) + log(|y(k2-k1)|
%        and for the phase: log(\tilde{B}(k1, k2)) = i (angle(y(k1)) - angle(y(k2)) + angle(y(k2 - k1)).
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion

[ I , J ] = ind2sub([d, d], 1:d^2);
K = J - I;
K = mod(K, d) + 1;
rows = repmat([1:d^2]', 3, 1);
cols = [I'; J'; K'];
cols = cols(:);
vals = [ ones(d^2, 1); -ones(d^2, 1); ones(d^2, 1) ];
M1 = sparse(rows, cols, ones(3*d^2, 1), d^2, d);
M2 = sparse(rows, cols, vals, d^2, d);

end




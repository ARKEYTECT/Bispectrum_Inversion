function B_mat = get_bispectrum(Y_fft,d,m)
% Y_fft is the column-wise fourier transformed sequence of the observed
% matrix of size d x m
% Given the parameter, this function first computes bipsectrum matrix for
% each column of the observation matrix. (See formula II.3 and II.4)
% The average of the bipsectrum is therefore computed across all m
% observations
%
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion
       B_mat = zeros(d);
       for i = 1:m
            B = (Y_fft(:,i) * Y_fft(:,i)') .* circulant(Y_fft(:,i));
            B_mat = B + B_mat;
       end
       B_mat = B_mat/m;
end
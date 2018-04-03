function [ meanB, B_mat ] = get_bispectrum_v2(Y, M1, M2)
% Y_fft is the column-wise fourier transformed sequence of the observed
% matrix of size d x m
% Given the parameter, this function first computes bipsectrum matrix for
% each column of the observation matrix. (See formula II.3 and II.4)
% The average of the bipsectrum is therefore computed across all m
% observations
% Input: 
%        Y: the Fourier Coefficients of the signal
%        M1 and M2: Matrix that multiply the Fourier coefficients together.
%        M1 and M2 are sparse matrices, each row only has 3 non-zeros.
%        The idea is for the Fourier amplitude, 
%        log(|B(k1, k2)|) = log(|y(k1)|) + log(|y(k2)|) + log(|y(k2-k1)|
%        and for the phase: log(\tilde{B}(k1, k2)) = i (angle(y(k1)) - angle(y(k2)) + angle(y(k2 - k1)).  
% Output:
%        meanB: sample mean of the bispectrum
%        B_mat: all the sample bispectrum
%
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion
tic
d = sqrt(size(M1, 1));
Y1 = log(abs(Y));
phase = atan2(imag(Y), real(Y));
t1 = toc;
tic
B_mat = exp(M1*Y1+sqrt(-1)*M2*phase);
t2 = toc;
tic;
meanB = mean(B_mat, 2);
meanB = reshape(meanB, d, d);
t3 = toc;

end
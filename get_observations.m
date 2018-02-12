function [obs,shifts] = get_observations(x_real, sigma, m, d)
% Given the signal x_real of length d, the function generates a matrix of 
% observed copies of size d x m. Each column of the matrix is a randomly
% and circularly shifted version of x_real with i.i.d Gaussian noise with
% variance equal to sigma^2.
% 
% Second output shifts stores the true shift sequence for every column of obs
% shifts is used to compute known shift oracle 
%
% Jan 2018
% Hua Chen 
% https://github.com/ARKEYTECT/Bispectrum_Inversion

    obs = zeros(d,m);
    
    shifts = randi(d, m, 1) - 1;
    for i = 1 : m
        obs(:, i) = circshift(x_real, shifts(i));
    end
    
    obs = obs +randn(d,m).*sigma; 
    
end

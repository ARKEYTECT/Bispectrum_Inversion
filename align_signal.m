function x_aligned = align_signal(x, x_real)
% Given two column vectors x and xref, returns x_aligned after circularly shifting
% x such that it is optimally aligned with xref in 2-norm.
%
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion

    x_fft = fft(x);
    x_real_fft = fft(x_real);
    
    correlation = real(ifft(conj(x_fft) .* x_real_fft));
    [~, ind] = max(correlation);
    x_aligned = circshift(x, ind-1);
    
end
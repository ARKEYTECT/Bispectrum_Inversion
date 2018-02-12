function X_recon = get_recon(power, phase, X_real)
% Given power spectrum and phase of the estiamted signal, reconstruct the
% signal and align it to X_real
%
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion

        X_power = sqrt(power).*phase; 
        X_fft = real(ifft(X_power));
        X_recon = align_signal(X_fft,X_real);
end

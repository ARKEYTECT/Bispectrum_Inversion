function zout = phases_from_bispectrum_FM_real(B, z0, z1)
% Frequency marching to estimate the phases of the DFT of a signal x given
% an estimator for its bispectrum, B. Entries z0 and z1, if specified,
% provide the phases of the first two entries of the DFT of x, as
% unit-modulus complex numbers (z0 is the phase of the mean, z1 is the
% phase of the first frequency.)
%
% Jan 2018
% Reference to Nocolas Boumal in https://github.com/NicolasBoumal/MRA
% https://github.com/ARKEYTECT/Bispectrum_Inversion

    N = size(B, 1);
    assert(size(B, 2) == N, 'B must be square.');
    
    % For real signals, B is supposed to be Hermitian.
    B = (B+B')/2;

    if ~exist('z0', 'var') || isempty(z0)
        z0 = sign(B(1, 1));
    end
    if ~exist('z1', 'var') || isempty(z1)
        z1 = 1;
    end
    if z0 == 0
        z0 = 1;
    end
    if z1 == 0
        z1 = 1;
    end
    
    Psi = angle(B);

    psi = zeros(floor((N+1)/2), 1);
    psi(1) = angle(z0);
    psi(2) = angle(z1);


    for ii = 3 : (N+1)/2

        psi_est = zeros(ceil(ii/2)-1,1);

        for jj=2:ceil(ii/2)

            psi_est(jj-1) =  Psi( ii, jj ) + psi(jj) + psi(ii-jj+1);

        end

        % Averaging all estimation over SO(2).
        psi(ii) = angle(sum(exp(1j*psi_est)));

    end

    if N/2 == floor(N/2)
        psi_mid = 0;
    else
        psi_mid = [];
    end
    
    % Use symmetries due to the signal being real.
    psi = [angle(z0) ; psi(2:end) ; psi_mid ; flipud(-psi(2:end))];

    zout = exp(1i*psi);

end

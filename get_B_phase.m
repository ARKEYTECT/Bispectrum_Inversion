function B_phase = get_B_phase(B_mat)
% Given the bispectrum matrix, compute the phase of the matrix
% Assign 1 to overly great terms (line 10)
%
% Jan 2018 
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion
 
        B_phase = B_mat./abs(B_mat);
        B_phase(abs(B_mat)<1*10^-12) = 1;
        
end
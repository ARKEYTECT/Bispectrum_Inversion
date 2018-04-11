function [ meanB, toc_bispec ] = get_bispectrum(Y)
% Y is the real signal of size dxn 
% Using 12-fold symmetry of bispoectrum to optimize implement
% Input: 
%        Y: real signals after removing the mean
%         
% Output:
%        meanB: sample mean of the bispectrum
%
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion

% build 3rd moment projection
tic_start = tic;
[ d] = size(Y, 1);
[ I, J ] = ind2sub([d, d], 1:d^2);
I = I.';
J = J.';
K = mod(J - I, d) + 1;
ind1 = find( I >= 0.5*J + 0.5 & I <= J & I <= d - J + 2 );
I1 = I(ind1);
J1 = J(ind1);
K1 = K(ind1);
Btmp = mean(Y(I1, :).*conj(Y(J1, :)).*Y(K1, :), 2);

B = zeros(d^2, 1);
B(ind1) = Btmp;

%Apply symmetries
I1 = I1 - 1;
J1 = J1 - 1;
Itmp = mod(-J1, d) + 1;
Jtmp = mod(-I1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = Btmp;

Itmp = mod(J1 - I1, d) + 1;
Jtmp = mod(-I1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = Btmp;

Itmp = mod(J1 - I1, d) + 1;
Jtmp = mod(J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = Btmp;

Itmp = mod(I1, d) + 1;
Jtmp = mod(I1 - J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = Btmp;

Itmp = mod(-J1, d) + 1;
Jtmp = mod(I1 - J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = Btmp;

Itmp = mod(-I1, d) + 1;
Jtmp = mod(-J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = conj(Btmp);

Itmp = mod(J1, d) + 1;
Jtmp = mod(I1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = conj(Btmp);

Itmp = mod(-J1 + I1, d) + 1;
Jtmp = mod( I1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = conj(Btmp);

Itmp = mod(-J1 + I1, d) + 1;
Jtmp = mod( -J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = conj(Btmp);

Itmp = mod( -I1, d) + 1;
Jtmp = mod( -I1 + J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = conj(Btmp);

Itmp = mod( J1, d) + 1;
Jtmp = mod( - I1 + J1, d) + 1;
ind2 = sub2ind([d, d], Itmp, Jtmp);
B(ind2) = conj(Btmp);

meanB = reshape(B, d, d);

toc_bispec = toc(tic_start);

end
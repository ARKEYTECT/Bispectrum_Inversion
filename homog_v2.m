function [ RelativeError,t1,t2 ] = homog_v2(bias, x,Y, L, sigma, n )
%RUN_HOMOG Generates and solves a homogeneous MRA problem via Jennrich
%   Detailed explanation goes here

% The data. Let's not even bother with random rotations since we only
% symmetrize anyway; it saves a ton of time not generating Y sample by
% sample
%x = randn(L,1)/sqrt(L); % legit version
% x = zeros(L,1); % to debug bias
% x = 2*ones(L,1); % to debug bias
%svar = var(x);
%sigma = sqrt(svar/SNR);
% Y = sigma*randn(L,n) + x * ones(1,n);

% Compute empirical mean, debiased norm, and symmetrized 3rd moment. These
% are rotationally invariant, and we delete everything just after to
% enforce non-cheating

tic;
Emean = mean(mean(Y)); % first moment

%TM = mean(x(:))

% second moment
EM2_asymm = Y*Y'/(L*n);
EM2 = zeros(L,L);
for i = 1:L % circulize
    EM2 = EM2 + circshift(circshift(EM2_asymm,i,1),i,2);
end
M2 = EM2 - sigma^2*eye(L); % debias
%M2 = EM2 - sigma^2*eye; % incorrect but mystically better version
diagM2 = max(0, diag(M2));
M2 = M2 + diag(diagM2 - diag(M2));
normest = sqrt(sum(diagM2));
%normest = sqrt(trace(M2)); % estimate 2-norm


z = randn(L,1)/sqrt(L); % projection direction


% EM3 = zeros(L,L,L);
% if L < 10 % more efficient method for small L: symmetrize the data
%     for sh = 1:L
%         Ysh = circshift(Y,sh);
%         for h = 1:L
%             EM3(:,:,h) = EM3(:,:,h) + (Ysh .* (ones(L,1)*Ysh(h,:))) * Ysh' / (n*L);
%         end
%     end
% else % more efficient for L biggish (but plenty smaller than n): write third moment and then symmetrize it
%     EM3_asymm = zeros(L,L,L);
%     for h = 1:L
%         EM3_asymm(:,:,h) = (Y .* (ones(L,1)*Y(h,:))) * Y' / (n*L);
%     end
%     % symmetrize
%     for h1 = 1:L
%         for h2 = 1:L
%             s = 0;
%             for h3 = 1:L
%                 s = s + EM3_asymm(mod(h1+h3-1,L)+1,mod(h2+h3-1,L)+1,h3);
%             end
%             for h3 = 1:L
%                 EM3(mod(h1+h3-1,L)+1,mod(h2+h3-1,L)+1,h3) = s;
%             end
%         end
%     end
% end
% 
% % unbias third moment
% Id = eye(L);
% T_bias = zeros(L,L,L);
% for h=1:L
%     T_bias(:,:,h) = sigma^2*Emean*(Id + ones(L,1)*Id(:,h)' + Id(:,h)*ones(1,L));
% end
% M3 = EM3 - T_bias;

% project
% EM3pv1 = reshape(reshape(EM3,[L^2 L]) * z, [L L]);
% M3pv1 = reshape(reshape(M3,[L^2 L]) * z, [L L]);



% build 3rd moment projection
EM3p = zeros(L,L);
for sh = 1:L
    delt = (Y .* (ones(L,1)*(circshift(z,-sh)'*Y))) * Y';
    EM3p = EM3p + circshift(circshift(delt,sh,1),sh,2);
end
EM3p = EM3p / (n*L);

% debias
M3p = EM3p - sigma^2*Emean* (sum(z)*eye(L) + ones(L,1)*z' + z*ones(1,L));
% TODO: check M3p bias & error

%{
% for debugging bias
Emean
M2
M3p
sum(z)
%}
t1 = toc;



tic;
% eigs of pseudoinverse
G = M3p * pinv(M2);
%G = M3p * pinv(EM2 + 100*eye(L));
%{
[V,E,U] = svd(G);
E
U
V
x/norm(x)
V = V(:,1);
%}
try
    [V,~] = eigs(G*G,1,'lr');
catch ERR % occasionally the eigenvalue solver throws a convergence error
    'Convergence error'
    RelativeError=1.0;
    return;
end
V = real(V); % proportional to x in the large SNR limit

% {
% estimate the right sign: compute what M3 would be if V were the true
% signal, and compare to actual M3
% M3fake = zeros(L,L,L);
% for h=1:L
%     r = circshift(V,h);
%     for hh=1:L
%         M3fake(:,:,hh) = M3fake(:,:,hh) + r*r'*r(hh) / L;
%     end
% end
% % fix the sign
% V = V * sign( sum(sum(dot(M3,M3fake))) );
% }

% fix the sign via the mean
V = V * sign( Emean * mean(V));

%added by Jane
% V = V - mean(V) + Emean;

% fix the 2-norm
V = V * normest / norm(V);
t2 = toc;
V = V - bias;

minimum_error = Inf;
for h=1:L
   V_new = circshift(V,h);
   V_new_2 = -circshift(V,h);
   error1 = norm(V_new - x);
   error2 = norm(V_new_2-x);
   if error1 > error2
       V_new = V_new_2;
   end
   
   if norm(V_new - x)< minimum_error
      minimum_error = norm(V_new - x);
      V_best = V_new;
   end
end

% V = V * sign( Emean * mean(V));
% 
% fix the 2-norm
% V = V * normest / norm(V);
% 
% minimum_error = Inf;
% for h=1:L
%    V_new = circshift(V,h);
%    if norm(V_new - x)< minimum_error
%       minimum_error = norm(V_new - x);
%       V_best = V_new;
%    end
% end

Estimate = V_best;
GroundTruth = x;
RelativeError = norm(Estimate-GroundTruth)/norm(GroundTruth);
end


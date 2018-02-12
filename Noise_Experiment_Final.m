% Run line 150 - 183 to generate graph in the paper 
% with data given in .mat file

% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion
clc;
clear;
close all;

noise = 0.01*2.^(0:11);
sigma = sqrt(noise);
d = 41; 

r_error_1 = zeros(1,length(noise));
r_error_2 = zeros(1,length(noise));
r_error_3 = zeros(1,length(noise));
% r_error_4 = zeros(1,length(noise));
r_error_5 = zeros(1,length(noise));
r_error_6 = zeros(1,length(noise));
r_error_7 = zeros(1,length(noise));
r_error_8 = zeros(1,length(noise));
runtime = zeros(7,length(noise));

iter = 50;
copy = 10000;
seed = 1234;
rng(seed);

X_real = randn(d,1);
X_real = X_real - mean(X_real);
% compute oracle

for  a = 1:12
    fprintf('Noise Level: %d\n', sigma(a));
    for b = 1:iter  
        fprintf('Iteration: %d\n', b);
        % generate observations 
        [Y,shifts] = get_observations(X_real, sigma(a), copy, d);
        mean_est = mean(Y);
        mean_est = repmat(mean_est,d,1);
        Y = Y - mean_est;
        
        % compute DFT for every signal   
        Y_hat = fft(Y,[],1);    
        
        % compute Bispectrum matrix
        B_mat = get_bispectrum(Y_hat,d,copy);
        % Power estimate and phase estimate 
        Y_power = mean(abs(Y_hat).^2,2)-d*sigma(a)^2;
        Y_power = max(0, Y_power);
        B_phase = get_B_phase(B_mat);
        
        %%% Phase Recovery & Reconstruction
        
        % phase with greatest gap
        tic
        Est_phase_1 = get_phase_from_bispectrum_gap(B_phase,d); 
        runtime(1,a) = runtime(1,a) + toc; 
        
        % phase with greatest absolute eigenvector
%         tic
%         Est_phase_4 = get_phase_from_bispectrum_absval(B_phase,d); 
%         runtime(2,a) = runtime(2,a) + toc; 
        
        % phase from iterative phase sync
        tic
        Est_phase_2 = phases_from_bispectrum_APS_real(B_mat); 
        runtime(3,a) = runtime(3,a) + toc; 
        
        % Optim. on phase manifold
        tic
        [Est_phase_3, problem] = phases_from_bispectrum_real(B_mat); 
        runtime(4,a) = runtime(4,a) + toc; 
        
        % Freq. marching
        y = fft(X_real);
        tic
        Est_phase_5 = phases_from_bispectrum_FM_real(B_mat,sign(y(1)), sign(y(2))); 
        runtime(5,a) = runtime(5,a) + toc; 
        
        % Phase Unwrapping
        y = fft(X_real);
        tic
        Est_phase_6 = phases_from_bispectrum_LLL_real(B_mat,sign(y(1)), sign(y(2))); 
        runtime(6,a) = runtime(6,a) + toc; 
        
        % SDP
        y = fft(X_real);
        tic
        Est_phase_7 = phases_from_bispectrum_SDP_real(B_mat,sign(y(1)), sign(y(2))); 
        runtime(7,a) = runtime(7,a) + toc; 
        
        %%% Compute X_oracle
        X_oracle = zeros(d, copy);
        for m = 1 : copy
            X_oracle(:, m) = circshift(Y(:, m), -shifts(m));
        end
        X_oracle = mean(X_oracle, 2);
        %%% 
        
        X_Recon_1 = get_recon(Y_power,Est_phase_1,X_real);
        X_Recon_2 = get_recon(Y_power,Est_phase_2,X_real);
        X_Recon_3 = get_recon(Y_power,Est_phase_3,X_real);
%         X_Recon_4 = get_recon(Y_power,Est_phase_4,X_real);
        X_Recon_5 = get_recon(Y_power,Est_phase_5,X_real);
        X_Recon_6 = get_recon(Y_power,Est_phase_6,X_real);
        X_Recon_7 = get_recon(Y_power,Est_phase_7,X_real);
        
        %%%%% Validation 
        relative_error = norm(X_Recon_1-X_real,2)/norm(X_real,2);
        r_error_1(a) = r_error_1(a) + relative_error;
        
        relative_error = norm(X_Recon_2-X_real,2)/norm(X_real,2);
        r_error_2(a) = r_error_2(a) + relative_error;
        
        relative_error = norm(X_Recon_3-X_real,2)/norm(X_real,2);
        r_error_3(a) = r_error_3(a) + relative_error;
        
%         relative_error = norm(X_Recon_4-X_real,2)/norm(X_real,2);
%         r_error_4(a) = r_error_4(a) + relative_error;
        
        relative_error = norm(X_Recon_5-X_real,2)/norm(X_real,2);
        r_error_5(a) = r_error_5(a) + relative_error;
        
        relative_error = norm(X_Recon_6-X_real,2)/norm(X_real,2);
        r_error_6(a) = r_error_6(a) + relative_error;
        
        relative_error = norm(X_Recon_7-X_real,2)/norm(X_real,2);
        r_error_7(a) = r_error_7(a) + relative_error;
        
        relative_error = norm(X_oracle-X_real,2)/norm(X_real,2);
        r_error_8(a) = r_error_8(a) + relative_error;
        
        
        
    end
    r_error_1(a) = r_error_1(a)/iter; 
    r_error_2(a) = r_error_2(a)/iter;
    r_error_3(a) = r_error_3(a)/iter;
%     r_error_4(a) = r_error_4(a)/iter;
    r_error_5(a) = r_error_5(a)/iter;
    r_error_6(a) = r_error_6(a)/iter;
    r_error_7(a) = r_error_7(a)/iter;
    r_error_8(a) = r_error_8(a)/iter;
    runtime(:,a) = runtime(:,a)./iter;
end



figure;
loglog(noise,r_error_1,'k');
hold all;
% loglog(noise,r_error_4,'k');
loglog(noise,r_error_2,'g');
loglog(noise,r_error_3,'Color',[1.0 0.5 0.0]);
loglog(noise,r_error_5,'b--');
loglog(noise,r_error_6,'c--');
loglog(noise,r_error_7,'--','Color',[0.9 0.4 0.4]);
loglog(noise,r_error_8,'m');
axis([0.01 40 0.0008 10]);
xlabel('Noise level \sigma');
ylabel('Relative Error (up to circular shift)');
legend('Spectral M. largest spectral gap','Iterative phase sync.','Optim. on phase manifold', ...,
'FM','Phase unwrapping','SDP','Known-shifts oracle','Location','best');
set(gca, 'FontSize', 12)
legend('boxoff')

figure;
loglog(noise,runtime(1,:),'k');
hold all;
% loglog(noise,runtime(2,:),'k');
loglog(noise,runtime(3,:),'g');
loglog(noise,runtime(4,:),'Color',[1.0 0.5 0.0]);
loglog(noise,runtime(5,:),'b--');
loglog(noise,runtime(6,:),'c--');
loglog(noise,runtime(7,:),'--','Color',[0.9 0.4 0.4]);
xlabel('Noise level \sigma');
ylabel('Computation time(s)');
axis([0.01 25 0.0001 100000]);
legend('Spectral M. largest spectral gap','Iterative phase sync.','Optim. on phase manifold', ...,
'FM','Phase unwrapping','SDP','Location','best');
set(gca, 'FontSize', 12)
legend('boxoff')
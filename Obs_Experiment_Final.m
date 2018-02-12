% Run line 139-172 to generate graph in the paper 
% with data given in .mat file

% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion 
clc;
clear;
close all;

d = 41; 
obs = [10 100 1000 10000 100000];

r_error_1 = zeros(1,length(obs));
r_error_2 = zeros(1,length(obs));
r_error_3 = zeros(1,length(obs));
% r_error_4 = zeros(1,length(obs));
r_error_5 = zeros(1,length(obs));
r_error_6 = zeros(1,length(obs));
r_error_7 = zeros(1,length(obs));
r_error_8 = zeros(1,length(obs));
runtime = zeros(7,length(obs));

iter = 50;
seed = 1234;
rng(seed);
X_real = randn(d,1);
X_real = X_real - mean(X_real);


for  a = 1:length(obs) 
    copy = obs(a);
    fprintf('Observation: %d\n', obs(a));
    for b = 1:iter       
        fprintf('Iteration: %d\n', b);
        % generate observations 
        [Y,shifts] = get_observations(X_real, 1, copy, d);
        mean_est = mean(Y);
        mean_est = repmat(mean_est,d,1);
        Y = Y - mean_est;
        
        % compute DFT for every signal   
        Y_hat = fft(Y,[],1);  
        
        % compute Bispectrum matrix
        B_mat = get_bispectrum(Y_hat,d,copy);
        % Power estimate and phase estimate 
        Y_power = mean(abs(Y_hat).^2,2)-d;
        Y_power = max(0, Y_power);
        B_phase = get_B_phase(B_mat);
        
        %%% Phase Recovery & Reconstruction
        
        % phase with greatest gap
        tic
        Est_phase_1 = get_phase_from_bispectrum_gap(B_phase,d); 
        runtime(1,a) = runtime(1,a) + toc; 
        
%         % phase with greatest abs value 
%         tic
%         Est_phase_4 = get_phase_from_bispectrum_absval(B_phase,d); 
%         runtime(2,a) = runtime(2,a) + toc; 
        
        % phase from iterative phase sync
        tic
        Est_phase_2 = phases_from_bispectrum_APS_real(B_mat); 
        runtime(3,a) = runtime(3,a) + toc; 
        tic
        [Est_phase_3, problem] = phases_from_bispectrum_real(B_mat); % Optim. on phase manifold
        runtime(4,a) = runtime(4,a) + toc; 
        y = fft(X_real);
        tic
        Est_phase_5 = phases_from_bispectrum_FM_real(B_mat,sign(y(1)), sign(y(2))); % Freq. marching
        runtime(5,a) = runtime(5,a) + toc; 
        y = fft(X_real);
        tic
        Est_phase_6 = phases_from_bispectrum_LLL_real(B_mat,sign(y(1)), sign(y(2))); % Unwrapping
        runtime(6,a) = runtime(6,a) + toc; 
        y = fft(X_real);
        tic
        Est_phase_7 = phases_from_bispectrum_SDP_real(B_mat,sign(y(1)), sign(y(2))); % SDP
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
loglog(obs,r_error_1,'r'); 
hold all;
% loglog(obs,r_error_4,'k');
loglog(obs,r_error_2,'g');
loglog(obs,r_error_3,'Color',[1.0 0.5 0.0]);
loglog(obs,r_error_5,'b--');
loglog(obs,r_error_6,'c--');
loglog(obs,r_error_7,'--','Color',[0.9 0.4 0.4]);
loglog(obs,r_error_8,'m');
xlabel('# Observations M');
ylabel('Relative Error (up to circular shift)');
legend('Spectral M. largest spectral gap','Spectral M. greatest absolute eig-value','Iterative phase sync.','Optim. on phase manifold', ...,
'FM','Phase unwrapping','SDP','Known-shifts oracle','Location','best');
set(gca, 'FontSize', 12)
legend('boxoff');

figure;
loglog(obs,runtime(1,:),'r');
hold all;
loglog(obs,runtime(2,:),'k');
loglog(obs,runtime(3,:),'g');
% loglog(obs,runtime(4,:),'Color',[1.0 0.5 0.0]);
loglog(obs,runtime(5,:),'b--');
loglog(obs,runtime(6,:),'c--');
loglog(obs,runtime(7,:),'--','Color',[0.9 0.4 0.4]);
xlabel('# Observations M');
ylabel('Computation time(s)');

axis([10 100000 0.0001 100000]);
legend('Spectral M. largest spectral gap','Spectral M. greatest absolute eig-value','Iterative phase sync.','Optim. on phase manifold', ...,
'FM','Phase unwrapping','SDP','Location','best');
set(gca, 'FontSize', 12)
legend('boxoff');

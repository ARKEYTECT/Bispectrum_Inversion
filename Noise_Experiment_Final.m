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
iter = 50;

r_error_1 = zeros(1,length(noise));
r_error_2 = zeros(1,length(noise));
r_error_3 = zeros(1,length(noise));
r_error_4 = zeros(1,length(noise));
r_error_5 = zeros(1,length(noise));
r_error_6 = zeros(1,length(noise));
r_error_7 = zeros(1,length(noise));
r_error_8 = zeros(1,length(noise));
runtime1 = zeros(iter+1,length(noise));
runtime2 = zeros(iter+1,length(noise));
runtime3 = zeros(iter+1,length(noise));
runtime4 = zeros(iter+1,length(noise));
runtime5 = zeros(iter+1,length(noise));
runtime6 = zeros(iter+1,length(noise));
runtime7 = zeros(iter+1,length(noise));
% gap_runtime = zeros(2,length(sigma));
% jenn_runtime = zeros(2,length(sigma));


copy = 10000;
seed = 1234;
rng(seed);

X_real = randn(d,1);
X_real = X_real - mean(X_real);

bias = 0.1;
x_jenn = X_real + bias;

[ M1, M2 ] = bispec_precomp(d);
 
for  a = 1:length(noise)
    fprintf('Noise Level: %d\n', sigma(a));
    for b = 1:iter  
        fprintf('Iteration: %d\n', b);
        % jennrich algorithm 
        [Y,~] = get_observations(x_jenn, sigma(a), copy, d);
        tic
        [ relerr, t1, t2 ] = homog_v2(bias,X_real,Y, d, sigma(a), copy );
        runtime2(b+1,a) = t2;
        r_error_4(1,a) = r_error_4(1,a) + relerr;
        
        % Spectral gap method
        
        [Y,shifts] = get_observations(X_real, sigma(a), copy, d); 
        
        tic;
        Emean = mean(mean(Y));
        mean_est = mean(Y);
        mean_est = repmat(mean_est,d,1);
        Y = Y - mean_est;
        Y_hat = fft(Y,[],1);   
        Y_power = mean(abs(Y_hat).^2,2)-d*sigma(a)^2;
        Y_power = max(0, Y_power);
        t_power = toc;

        [ meanB, ~ ] = get_bispectrum_v2(Y_hat, M1, M2);
        tic;
        B_phase = get_B_phase(meanB);
        t_bphase = toc;
        
        %%% Phase Recovery & Reconstruction
        
        % phase with greatest gap
        tic
        Est_phase_1 = get_phase_from_bispectrum_gap(B_phase,d); 
        X_Recon_1 = get_recon(Emean,Y_power,Est_phase_1);
        runtime1(b+1,a) = t_power+t_bphase+toc; 
        X_Recon_1 = align_signal(X_Recon_1,X_real);
        
        % phase from iterative phase sync
        tic
        Est_phase_2 = phases_from_bispectrum_APS_real(meanB); 
        X_Recon_2 = get_recon(Emean,Y_power,Est_phase_2);
        runtime3(b+1,a) = t_power+t_bphase+toc; 
        X_Recon_2 = align_signal(X_Recon_2,X_real);
        
        % Optim. on phase manifold
        tic
        [Est_phase_3, problem] = phases_from_bispectrum_real(meanB,1); 
        X_Recon_3 = get_recon(Emean,Y_power,Est_phase_3);
        runtime4(b+1,a) = t_power+t_bphase+toc; 
        X_Recon_3 = align_signal(X_Recon_3,X_real);
        
        % Freq. marching
        y = fft(X_real);
        tic
        Est_phase_5 = phases_from_bispectrum_FM_real(meanB,sign(y(1)), sign(y(2))); 
        X_Recon_5 = get_recon(Emean,Y_power,Est_phase_5);
        runtime5(b+1,a) = t_power+t_bphase+toc; 
        X_Recon_5 = align_signal(X_Recon_5,X_real);
        
        % Phase Unwrapping
        y = fft(X_real);
        tic
        Est_phase_6 = phases_from_bispectrum_LLL_real(meanB,sign(y(1)), sign(y(2))); 
        X_Recon_6 = get_recon(Emean,Y_power,Est_phase_6);
        runtime6(b+1,a) = t_power+t_bphase+toc; 
        X_Recon_6 = align_signal(X_Recon_6,X_real);
        
        % SDP
        y = fft(X_real);
        tic
        Est_phase_7 = phases_from_bispectrum_SDP_real(meanB,sign(y(1)), sign(y(2))); 
        X_Recon_7 = get_recon(Emean,Y_power,Est_phase_7);
        runtime7(b+1,a) = t_power+t_bphase+toc; 
        X_Recon_7 = align_signal(X_Recon_7,X_real);
        
        %%% Compute X_oracle
        X_oracle = zeros(d, copy);
        for m = 1 : copy
            X_oracle(:, m) = circshift(Y(:, m), -shifts(m));
        end
        X_oracle = mean(X_oracle, 2);
        %%% 
        
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
    r_error_4(a) = r_error_4(a)/iter;
    r_error_5(a) = r_error_5(a)/iter;
    r_error_6(a) = r_error_6(a)/iter;
    r_error_7(a) = r_error_7(a)/iter;
    r_error_8(a) = r_error_8(a)/iter;
    runtime1(1,a) = median(runtime1(2:end,a));
    runtime2(1,a) = median(runtime2(2:end,a));
    runtime3(1,a) = median(runtime3(2:end,a));
    runtime4(1,a) = median(runtime4(2:end,a));
    runtime5(1,a) = median(runtime5(2:end,a));
    runtime6(1,a) = median(runtime6(2:end,a));
    runtime7(1,a) = median(runtime7(2:end,a));
    
end



figure;
loglog(noise,r_error_1,'k');
hold all;
% loglog(noise,r_error_4,'k');
loglog(noise,r_error_4,'r');
loglog(noise,r_error_2,'g');
loglog(noise,r_error_3,'Color',[1.0 0.5 0.0]);
loglog(noise,r_error_5,'b--');
loglog(noise,r_error_6,'c--');
loglog(noise,r_error_7,'--','Color',[0.9 0.4 0.4]);
loglog(noise,r_error_8,'m');
axis([0.01 30 0.0008 10]);
xlabel('Noise level \sigma^2');
ylabel('Relative Error (up to circular shift)');
legend('Spectral M. largest spectral gap','Jennrich Algorithm','Iterative phase sync.','Optim. on phase manifold', ...,
'FM','Phase unwrapping','SDP','Known-shifts oracle','Location','best');
set(gca, 'FontSize', 12)
legend('boxoff')

figure;
loglog(noise,runtime1(1,:),'k','linewidth',1);
hold all;
loglog(noise,runtime2(1,:),'r','linewidth',1);
% loglog(noise,runtime(2,:),'k');
% loglog(noise,gap_runtime(2,:),'m');
% loglog(noise,jenn_runtime(2,:),'Color',[0.5,0.5,0.6]);
loglog(noise,runtime3(1,:),'g','linewidth',1);
loglog(noise,runtime4(1,:),'Color',[1.0 0.5 0.0],'linewidth',1);
loglog(noise,runtime5(1,:),'b--','linewidth',1);
loglog(noise,runtime6(1,:),'c--','linewidth',1);
loglog(noise,runtime7(1,:),'--','Color',[0.9 0.4 0.4],'linewidth',1);
xlabel('Noise level \sigma^2');
ylabel('Total computation time(s)');
axis([0.01 25 0.001 100000]);
lgd = legend('Spectral method','Jennrich Algorithm',...,
       'Iterative phase sync.','Optim. on phase manifold',..., 
       'FM','Phase unwrapping','SDP','Location','best');
set(gca, 'FontSize', 12);
legend('boxoff');
% lgd.FontSize = 8;
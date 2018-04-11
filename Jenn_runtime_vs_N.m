%%runtime vs n
clc; clear;

%{comments 
% 1. Noise Varaince 
% check zero frequency 
% 2. runtime 
% start from line 40 
%}
% d = 41;
d = 21:20:501;
iter = 40;
obs = [10 40 70 100 150 200 600 1000 5000 10000 100000];

rng('default');

% x = x - mean(x);
% bias = 0.5;
% x_jenn = x + bias;
% % bias = 0;
% noise = 0.01*2.^(0:11);
% sigma = sqrt(noise);

gap_runtime = zeros(iter+1,length(d));
jenn_runtime = zeros(iter+1,length(d));

r_error_gap = zeros(1,length(d));


for  i = 1:length(d) 
%      copy = obs(i);
     copy = 100;
    x = randn(d(i), 1);
    fprintf('Round: %d\n', i);
    for b = 1:iter       
    [Y,~] = get_observations(x, 0.1, copy, d(i)); 
    [ relerr, t1, t2 ] = homog_v2(0,x,Y, d(i), 0.1, copy );
    jenn_runtime(b+1,i) = t1;
    
        tic_spectral = tic;
        Emean = mean(mean(Y));
        mean_est = mean(Y);
        mean_est = repmat(mean_est,d(i),1);
        Y = Y - mean_est;
        Y_hat = fft(Y,[],1);   
        Y_power = mean(abs(Y_hat).^2,2)-d(i)*1^2;
        Y_power = max(0, Y_power); 
        meanB = get_bispectrum(Y_hat);
        gap_runtime(b+1,i) = toc(tic_spectral);
        
%         B_phase = get_B_phase(meanB);
%         Est_phase_1 = get_phase_from_bispectrum_gap(B_phase,d(i)); 
%         X_Recon_1 = get_recon(Emean,Y_power,Est_phase_1);
%         X_Recon_1 = align_to_reference(X_Recon_1,x);
%         
%         relative_error_gap = norm(X_Recon_1-x)/norm(x);
%         r_error_gap(i) = r_error_gap(i) + relative_error_gap;
        
    end
%     r_error_gap(i) = r_error_gap(i)/iter;

    jenn_runtime(1,i) = median(jenn_runtime(2:end,i));
    gap_runtime(1,i) = median(gap_runtime(2:end,i));
    
end

figure;
plot(obs,gap_runtime(1,:),'r-','linewidth',1);
hold all;
plot(obs,jenn_runtime(1,:),'k--','linewidth',1);
xlabel('Length of signal d');
ylabel('Computation time (s)');
% axis([20 120 0 1.4]);
legend('bispectral invariant feature','homoJenn invariant feature','Location','northwest');
set(gca, 'FontSize', 12)
legend('boxoff');

figure;
plot(d,gap_runtime(1,:),'r-','linewidth',1);
hold all;
plot(d,jenn_runtime(1,:),'k--','linewidth',1);
xlabel('Length of signal d');
ylabel('Computation time (s)');
% axis([20 120 0 1.4]);
legend('bispectral invariant feature','homoJenn invariant feature','Location','northwest');
set(gca, 'FontSize', 12)
legend('boxoff');
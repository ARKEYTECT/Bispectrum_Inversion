%%%%% 1.1 
clc;
clear;
close all;

noise = 0.01*2.^(0:11);
sigma = sqrt(noise);
d = 41; 
obs = [10 100 1000 10000 100000];
%copy = 100000;
r_error_1 = 0;
r_error_2 = 0;
iter = 50;
copy = 10000;
X_real = randn(d,1);
X_real = X_real - mean(X_real);
% compute oracle

    %copy = obs(a);
%         result_1 = zeros(length(noise),1);
%         result_2 = zeros(length(noise),1);
        result_3 = zeros(length(noise),iter+1);
        result_4 = zeros(length(noise),iter+1);
        
        rand_runtime = zeros(length(noise),iter+1);
        spectral_runtime = zeros(length(noise),iter+1);
for a = 1:length(noise)
        for q = 1:iter
        % generate observations 
        
        Y = get_observations(X_real, sigma(a), copy, 41);
        

        mean_est = mean(Y);
        mean_est = repmat(mean_est,41,1);
        Y = Y - mean_est;
        % compute DFT for every signal 
        Y_hat = fft(Y,[],1);

        
        % compute Bispectrum matrix
        B_mat = get_bispectrum(Y_hat,41,copy);
        
        % Power estimate and phase estimate 
        
        Y_power = mean(abs(Y_hat).^2,2)-d*sigma(a)^2;
        Y_power = max(0, Y_power);
        B_phase = get_B_phase(B_mat);
        
        %%% Phase Recovery & Reconstruction
        tic;
        guess = get_phase_from_bispectrum_gap(B_phase,41); %% phase
        t2 = toc;
%         [Est_phase_1,cost_1,info_1] = phases_from_bispectrum_APS_real(B_mat);
%         [Est_phase_2,cost_2,info_2] = phases_from_bispectrum_APS_real(B_mat,mean(guess));
        tic;
        [Est_phase_3,cost_3,info_3] = phases_from_bispectrum_real(B_mat);
        rand_runtime(a,q+1) = rand_runtime(a,q+1) + toc;
        tic;
        [Est_phase_4,cost_4,info_4] = phases_from_bispectrum_real(B_mat,1,guess);
        spectral_runtime(a,q+1) = spectral_runtime(a,q+1) + t2+toc;
        
%         run_1 = struct2array(info_1);
%         run_2 = struct2array(info_2);
%         
%         run_1 = struct2array(info_1);
%         d = size(info_1,2);
%         l = length(run_1)/d;
%         result_1(a,1) = result_1(a,1) + run_1(d*l-l+1)+1;
%         
%         run_2 = struct2array(info_2);
%         d = size(info_2,2);
%         l = length(run_2)/d;
%         result_2(a,1) = result_2(a,1) + run_2(d*l-l+1)+1;


        run_3 = struct2array(info_3);
        d = size(info_3,2);
        l = length(run_3)/d;
        result_3(a,q+1) = result_3(a,q+1) + run_3(d*l-l+1)+1;
        
        run_4 = struct2array(info_4);
        d = size(info_4,2);
        l = length(run_4)/d;
        result_4(a,q+1) = result_4(a,q+1) + run_4(d*l-l+1)+1;
        end
%         result_1(a,1) = result_1(a,1)/50;
%         result_2(a,1) = result_2(a,1)/50;
        result_3(a,1) = median(result_3(a,2:end));
        result_4(a,1) = median(result_4(a,2:end));
        rand_runtime(a,1) = median(rand_runtime(a,2:end));
        spectral_runtime(a,1) = median(spectral_runtime(a,2:end));
end     

%         figure(1);
%         loglog(noise,result_1,'r-');
%         hold on;
%         loglog(noise,result_2,'b--');
%         xlabel('Noise level \sigma');
%         ylabel('Iterations to covnerge');
%         legend('Random Init.','Spectral Init.');
        
        
        cfigure = figure;
        yyaxis left;
        line1 = semilogx(noise,rand_runtime(:,1),'r-');
        hold on;
        line2 = semilogx(noise,spectral_runtime(:,1),'k-');
        xlabel('Noise level \sigma^2');
        ylabel('Total computation time(s)')
        axis([0.01 25 0 0.025]);
        set(gca, 'FontSize', 12);
        
        yyaxis right;
        line3 = semilogx(noise,result_3(:,1),'r--');
        hold on;
        line4 = semilogx(noise,result_4(:,1),'k--');
        ylabel('Iterations to ocnverge');
        axis([0.01 25 0 45]);
        lgd = legend([line1;line2;line3;line4],{'Random Init. runtime','Spectral Init. runtime','Random Init.iterations to converge','Spectral Init. iterations to converge'});
        set(gca, 'FontSize', 12);
        legend('boxoff');




% ========================================================================
% Ruslan Masharipov, June, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================
%
% The impact of amplitude differences on regression coefficients asymmetry.
% EXAMPLE #2: Adding correlated signals
%
% Step 1: Generate correlated time series for 30 ROIs 
% Step 2: Add correlated time series to the first 15 ROIs (Module #1)
% Step 3: Calculate correlation and regression coefficients

clear
close all

N_ROIs = 30;
dur = 400;

mu = ones(1,N_ROIs).*100;
sigma = ones(N_ROIs,N_ROIs).*0.1;
sigma(find(eye(size(sigma)))) = 0.5;

mu2 = ones(1,N_ROIs/2);
sigma2 = ones(N_ROIs/2,N_ROIs/2).*0.15;
sigma2(find(eye(size(sigma2)))) = 0.5;

for sub = 1:100

    % Generate correlated time series (ts) for 30 ROIs without amplitude difference (amp_diff)
    ts_without_amp_diff = mvnrnd(mu,sigma,dur);
    % Add correlated time series to the first 15 ROIs (Module #1)
    correlated_ts = [mvnrnd(mu2,sigma2,dur), zeros(dur,N_ROIs/2)];
    ts_with_amp_diff = ts_without_amp_diff + correlated_ts;
    
    %Calculate correlation and regression coefficients
    for i = 1:N_ROIs
        for j = 1:N_ROIs

            FCcorr1(:,:,sub) = corr(ts_without_amp_diff);
            FCcorr2(:,:,sub) = corr(ts_with_amp_diff);
            
            b1 = regress(ts_without_amp_diff(:,i),[ts_without_amp_diff(:,j) ones(dur,1)]);
            FCreg1(i,j,sub) = b1(1);
    
            b2 = regress(ts_with_amp_diff(:,i),[ts_with_amp_diff(:,j) ones(dur,1)]);
            FCreg2(i,j,sub) = b2(1);
    
            amp_diff1(i,j,sub) = mean(ts_without_amp_diff(:,i)) - mean(ts_without_amp_diff(:,j));
            amp_diff2(i,j,sub) = mean(ts_with_amp_diff(:,i))  - mean(ts_with_amp_diff(:,j));
    
            clear b1 b2
        end
    end
    
    mean_amplitude1(sub) = mean(max(ts_without_amp_diff) - min(ts_without_amp_diff));
    mean_amplitude2(sub) = mean(max(ts_with_amp_diff) - min(ts_with_amp_diff));
       
    symm1(sub) = check_symmetry(FCreg1(:,:,sub));
    symm2(sub) = check_symmetry(FCreg2(:,:,sub));
        
    amp_diff_up1(:,sub) = upper_triangle(amp_diff1(:,:,sub));
    amp_diff_up2(:,sub) = upper_triangle(amp_diff2(:,:,sub));

    FCreg_diff1 = upper_triangle(FCreg1(:,:,sub)) - lower_triangle(FCreg1(:,:,sub));
    FCreg_diff2 = upper_triangle(FCreg2(:,:,sub)) - lower_triangle(FCreg2(:,:,sub));
       
    amp_diff_and_FCreg_diff1(sub) = corr(abs(amp_diff_up1(:,sub)),abs(FCreg_diff1)');
    amp_diff_and_FCreg_diff2(sub) = corr(abs(amp_diff_up2(:,sub)),abs(FCreg_diff2)');
end

group_mean_amplitude1 = mean(mean_amplitude1);
group_mean_amplitude2 = mean(mean_amplitude2);
group_mean_amp_diff1 = mean(mean(amp_diff_up1));
group_mean_amp_diff2 = mean(mean(amp_diff_up2));

%% Plot FC matrices and amplitude differences
max_ax = 0.3;
figure('Position', [400 250 800 250])
subplot(131); imagesc(mean(FCcorr1,3));   title('Correlation'); axis square; caxis([-max_ax max_ax]); set(gca,'XTick',[1 15 30], 'YTick', [1 15 30]); colorbar
subplot(132); imagesc(mean(FCreg1,3));    title('Regression');  axis square; caxis([-max_ax max_ax]); set(gca,'XTick',[1 15 30], 'YTick', [1 15 30]); colorbar
subplot(133); imagesc(mean(amp_diff1,3)); title('Amplitude differences'); axis square; caxis([-max_ax max_ax]); set(gca,'XTick',[1 15 30], 'YTick', [1 15 30]); colorbar
try
    sgtitle('Without amplitude differences')
catch
    suptitle('Without amplitude differences')
end
try
    colormap('turbo')
catch
    colormap('jet')
end
movegui('center')

figure('Position', [400 250 800 250])
subplot(131); imagesc(mean(FCcorr2,3)); title('Correlation'); axis square; caxis([-max_ax max_ax]); set(gca,'XTick',[1 15 30], 'YTick', [1 15 30]); colorbar
subplot(132); imagesc(mean(FCreg2,3)); title(['Regression']); axis square; caxis([-max_ax max_ax]); set(gca,'XTick',[1 15 30], 'YTick', [1 15 30]); colorbar
subplot(133); imagesc(mean(amp_diff2,3)); title('Amplitude differences'); axis square; caxis([-max_ax max_ax]); set(gca,'XTick',[1 15 30], 'YTick', [1 15 30]); colorbar
try
    sgtitle('Without amplitude differences')
catch
    suptitle('Without amplitude differences')
end
try
    colormap('turbo')
catch
    colormap('jet')
end
movegui('center')

group_symm1 = check_symmetry(mean(FCreg1,3));
group_symm2 = check_symmetry(mean(FCreg2,3));
fprintf(['Correlation between upper and lower diagonal elements WITHOUT amplitude differences: r = ' num2str(group_symm1,'%4.2f') ' \n']);
fprintf(['Correlation between upper and lower diagonal elements WITH amplitude differences: r = ' num2str(group_symm2,'%4.2f') ' \n']);

%% Plot histograms 
figure
subplot(121); hist(amp_diff_and_FCreg_diff1,7); set(gca,'XTick',[-1 1]); title(['Without ampl. diff. Mean r = ' num2str(mean(amp_diff_and_FCreg_diff1),'%4.2f')]); axis square; 
subplot(122); hist(amp_diff_and_FCreg_diff2,7); set(gca,'XTick',[-1 1]); title(['With ampl. diff. Mean r = ' num2str(mean(amp_diff_and_FCreg_diff2),'%4.2f')]); axis square; 
try
    sgtitle({'Correlations between amplitude differences','and regression coefficient differences across ROIs'})
catch
    suptitle({'Correlations between amplitude differences','and regression coefficient differences across ROIs'})
end
movegui('center')



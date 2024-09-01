function corrdiff_FIR(stat_path,sots_path,exp_folder,N,TR,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
load(sots_path);
load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(1,'%.3d') filesep  'VOI_ALL.mat']);

length_N = size(VOI,1);
imagei = 1:length_N;
imaget = (imagei)*TR;

% TaskA
for i = 1:length(onsets{1,1})
    TaskA_block_time(i,:) = (imaget>onsets{1,1}(i)+6 & imaget<=onsets{1,1}(i)+durations{1,1});
end
TaskA_block_time = sum(TaskA_block_time,1);

% TaskB
for i = 1:length(onsets{1,2})
    TaskB_block_time(i,:) = (imaget>onsets{1,2}(i)+6 & imaget<=onsets{1,2}(i)+durations{1,2});
end
TaskB_block_time = sum(TaskB_block_time,1);

% Rest
all_onsets = sort([onsets{1,1}; onsets{1,2}]);
dur = durations{1,1};
for i = 1:length(all_onsets)
    Rest_block_time(i,:) = (imaget>all_onsets(i)+dur+6 & imaget<=all_onsets(i)+2*dur);
end
Rest_block_time = sum(Rest_block_time,1);

for subji = 1:N
    load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);

    TS_TaskA = detrend(VOI(TaskA_block_time==1,:),'constant');
    TS_TaskB = detrend(VOI(TaskB_block_time==1,:),'constant');
    TS_Rest  = detrend(VOI(Rest_block_time==1,:),'constant');
    corr_TaskA = atanh(corr(TS_TaskA));
    corr_TaskB = atanh(corr(TS_TaskB));
    corr_Rest  = atanh(corr(TS_Rest));
    corr_diff = atanh(corr(TS_TaskA)) - atanh(corr(TS_TaskB));

    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'CorrDiff_FIR.mat'],'corr_*');       
    clear VOI TS_TaskA TS_TaskB TS_Rest corr_*
end
time=toc;
fprintf(['Correlation difference :: Done in: ' num2str(time) 's \n']);

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'CorrDiff_FIR.mat']);
    corrdiff_TaskA(:,:,subji) = corr_TaskA;
    corrdiff_TaskB(:,:,subji) = corr_TaskB;
    corrdiff_Rest(:,:,subji) = corr_Rest;
    corrdiff_TaskA_vs_TaskB(:,:,subji) = corr_diff;   
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'CorrDiff_FIR.mat'],'corrdiff*');
time = toc;
fprintf(['Merge Direct correlation difference matrices :: Done in: ' num2str(time) 's \n']); 

%% Results
[corrdiff_Rest_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_Rest,q_level);
[corrdiff_TaskA_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_TaskA,q_level);
[corrdiff_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_TaskB,q_level);
[corrdiff_TaskA_vs_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_TaskA_vs_TaskB,q_level);
gm_corrdiff_TaskA_vs_TaskB = mean(corrdiff_TaskA_vs_TaskB,3);
gm_corrdiff_TaskA_vs_TaskB(1:1+size(gm_corrdiff_TaskA_vs_TaskB,1):end) = 0;
corrdiff_TaskA_vs_TaskB_FDR(1:1+size(corrdiff_TaskA_vs_TaskB_FDR,1):end) = 0;

figure
try
    sgtitle('Direct correlation difference (After FIR task regression)')
catch
    suptitle('Direct correlation difference (After FIR task regression)')
end
subplot(241); imagesc(mean(corrdiff_Rest,3)); title('Rest');  axis square; caxis(max_ax(mean(corrdiff_Rest,3),1));
subplot(242); imagesc(mean(corrdiff_TaskA,3)); title('Task A');  axis square; caxis(max_ax(mean(corrdiff_TaskA,3),1));
subplot(243); imagesc(mean(corrdiff_TaskB,3)); title('Task B');  axis square; caxis(max_ax(mean(corrdiff_TaskB,3),1));
subplot(244); imagesc(gm_corrdiff_TaskA_vs_TaskB); title('Task AvsB');  axis square; caxis(max_ax(gm_corrdiff_TaskA_vs_TaskB,1));
subplot(245); imagesc(corrdiff_Rest_FDR); title('Rest FDR');  axis square; 
subplot(246); imagesc(corrdiff_TaskA_FDR); title('Task A FDR');  axis square; 
subplot(247); imagesc(corrdiff_TaskB_FDR); title('Task B FDR');  axis square; 
subplot(248); imagesc(corrdiff_TaskA_vs_TaskB_FDR); title('Task AvsB FDR');  axis square;

set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 
colormap(subplot(245),'parula') 
colormap(subplot(246),'parula')
colormap(subplot(247),'parula') 
colormap(subplot(248),'parula') 

%% Sensitivity and Specificity
[TPR_corrdiff, TNR_corrdiff] = TPR_TNR(corrdiff_TaskA_vs_TaskB_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'CorrDiff_FIR_RESULTS.mat'],'TPR*','TNR*');


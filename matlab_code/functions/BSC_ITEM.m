function BSC_ITEM(stat_path,sots_path,exp_folder,N,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, May, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic

% Stimulus onset times (SOTs)
load(sots_path);

% Number of events
E = length(onsets{1,1})+length(onsets{1,2}); 

for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'SPM.mat']);    
    if SPM.Sess(1).U(1).ons(1) < 0
        SPM.Sess(1).U(1).ons(1) = 0;
    elseif SPM.Sess(1).U(2).ons(1) < 0
        SPM.Sess(1).U(2).ons(1) = 0;
    end
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_get_defaults('cmdline',true);
    spm_get_defaults('stats.resmem',1);
    spm_get_defaults('stats.maxmem',2^34);
    spm_get_defaults('stats.fmri.ufp',1);
    ITEM_est_1st_lvl(SPM)
    clear SPM
    
    % Save betas
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'ITEM_est_1st_lvl' filesep 'GLM1.mat'])
    for j = 1:E
        beta_hdr = spm_vol([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep ...
            'ITEM_est_1st_lvl' filesep 'gamma_' num2str(j,'%.4d') '.nii']);
        betas_ITEM(GLM1.Sess.is(j),:) = spm_read_vols(beta_hdr)';
        clear beta_hdr
    end
    
    % Remove dir (to save drive space)
    rmdir([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'ITEM_est_1st_lvl'],'s');
    
    % Sort betas
    betas_TaskA = betas_ITEM(1:length(onsets{1,1}),:);
    betas_TaskB = betas_ITEM(length(onsets{1,1})+1:length(onsets{1,1})+length(onsets{1,2}),:);
    
    % Beta-series correlations
    BSC_ITEM_TaskA = atanh(corr(betas_TaskA));
    BSC_ITEM_TaskB = atanh(corr(betas_TaskB));
    BSC_ITEM_TaskA_vs_TaskB = BSC_ITEM_TaskA - BSC_ITEM_TaskB;
    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_ITEM.mat'],'BSC*','betas_ITEM');
    clear BSC* betas_ITEM

    fprintf(['BSC-ITEM :: Subject #' num2str(subji) ' :: Done in: ' num2str(toc) 's \n']);
end

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_ITEM.mat']);
    BSC_ITEM_TaskA_group(:,:,subji)  = BSC_ITEM_TaskA;
    BSC_ITEM_TaskB_group(:,:,subji)  = BSC_ITEM_TaskB;
    BSC_ITEM_TaskA_vs_TaskB_group(:,:,subji)  = BSC_ITEM_TaskA_vs_TaskB;
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_ITEM.mat'],'BSC*');
time = toc;
fprintf(['Merge BSC-ITEM matrices :: Done in: ' num2str(time) 's \n']); 

%% Results: BSC ITEM
[BSC_ITEM_TaskA_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_ITEM_TaskA_group,q_level);
[BSC_ITEM_TaskB_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_ITEM_TaskB_group,q_level);
[BSC_ITEM_TaskA_vs_TaskB_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_ITEM_TaskA_vs_TaskB_group,q_level);
gm_BSC_ITEM_TaskA_vs_TaskB = mean(BSC_ITEM_TaskA_vs_TaskB_group,3);
gm_BSC_ITEM_TaskA_vs_TaskB(1:1+size(gm_BSC_ITEM_TaskA_vs_TaskB,1):end) = 0;
BSC_ITEM_TaskA_vs_TaskB_group_FDR(1:1+size(BSC_ITEM_TaskA_vs_TaskB_group_FDR,1):end) = 0;

figure
subplot(231); imagesc(mean(BSC_ITEM_TaskA_group,3)); title('ITEM Task A');      axis square; caxis(max_ax(mean(BSC_ITEM_TaskA_group,3),1));
subplot(232); imagesc(mean(BSC_ITEM_TaskB_group,3)); title('ITEM Task B');      axis square; caxis(max_ax(mean(BSC_ITEM_TaskB_group,3),1));
subplot(233); imagesc(gm_BSC_ITEM_TaskA_vs_TaskB);   title('ITEM Task AvsB');   axis square; caxis(max_ax(gm_BSC_ITEM_TaskA_vs_TaskB,1));
subplot(234); imagesc(BSC_ITEM_TaskA_group_FDR);     title('ITEM Task A FDR');  axis square; 
subplot(235); imagesc(BSC_ITEM_TaskB_group_FDR);     title('ITEM Task B FDR');  axis square; 
subplot(236); imagesc(BSC_ITEM_TaskA_vs_TaskB_group_FDR); title('ITEM Task AvsB FDR');  axis square; 

sgtitle('Beta-series correlations: ITEM approach')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(234),'parula') 
colormap(subplot(235),'parula') 
colormap(subplot(236),'parula') 

%% Sensitivity and Specificity
[TPR_BSC_ITEM, TNR_BSC_ITEM] = TPR_TNR(BSC_ITEM_TaskA_vs_TaskB_group_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_ITEM_RESULTS.mat'],'TPR*','TNR*');
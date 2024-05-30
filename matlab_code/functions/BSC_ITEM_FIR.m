function BSC_ITEM_FIR(stat_path,sots_path,exp_folder,N,TR,model,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, May, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic

% Number of *.nii images per subject
dur = length(dir([stat_path filesep exp_folder filesep 'funct_images' filesep 'Sub_001*']));

% Stimulus onset times (SOTs)
load(sots_path);

% Number of events
E = length(onsets{1,1})+length(onsets{1,2}); 

spm('defaults','fmri');
spm_jobman('initcfg');

for subji = 1:N
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = ...
        {[stat_path filesep exp_folder filesep 'GLMs_after_FIR' filesep 'Sub_' num2str(subji,'%.3d')]};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    for image = 1:dur
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{image,1} = ...
            [stat_path filesep exp_folder filesep 'funct_images_after_FIR_task_regression' filesep 'Sub_' num2str(subji,'%.3d') '_Image_' num2str(image,'%.4d') '.nii,1'];
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {[sots_path]};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = model;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    batch{subji} = matlabbatch;
    clear matlabbatch
end

parfor subji = 1:N
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_get_defaults('cmdline',true);
    spm_get_defaults('stats.resmem',1);
    spm_get_defaults('stats.maxmem',2^34);
    spm_get_defaults('stats.fmri.ufp',1);
    spm_jobman('run',batch{subji});
end

for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs_after_FIR' filesep 'Sub_' num2str(subji,'%.3d') filesep 'SPM.mat']);    
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
    load([stat_path filesep exp_folder filesep 'GLMs_after_FIR' filesep 'Sub_' num2str(subji,'%.3d') filesep 'ITEM_est_1st_lvl' filesep 'GLM1.mat'])
    for j = 1:E
        beta_hdr = spm_vol([stat_path filesep exp_folder filesep 'GLMs_after_FIR' filesep 'Sub_' num2str(subji,'%.3d') filesep ...
            'ITEM_est_1st_lvl' filesep 'gamma_' num2str(j,'%.4d') '.nii']);
        betas_ITEM(GLM1.Sess.is(j),:) = spm_read_vols(beta_hdr)';
        clear beta_hdr
    end
       
    % Sort betas
    betas_TaskA = betas_ITEM(1:length(onsets{1,1}),:);
    betas_TaskB = betas_ITEM(length(onsets{1,1})+1:length(onsets{1,1})+length(onsets{1,2}),:);
    
    % Beta-series correlations
    BSC_ITEM_TaskA = atanh(corr(betas_TaskA));
    BSC_ITEM_TaskB = atanh(corr(betas_TaskB));
    BSC_ITEM_TaskA_vs_TaskB = BSC_ITEM_TaskA - BSC_ITEM_TaskB;
    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_ITEM_FIR.mat'],'BSC*','betas_ITEM');
    clear BSC* betas_ITEM

    fprintf(['BSC-ITEM (After FIR task regression) :: Subject #' num2str(subji) ' :: Done in: ' num2str(toc) 's \n']);
end

% Remove dir (to save drive space)
rmdir([stat_path filesep exp_folder filesep 'GLMs_after_FIR'],'s');

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_ITEM_FIR.mat']);
    BSC_ITEM_TaskA_group(:,:,subji)  = BSC_ITEM_TaskA;
    BSC_ITEM_TaskB_group(:,:,subji)  = BSC_ITEM_TaskB;
    BSC_ITEM_TaskA_vs_TaskB_group(:,:,subji)  = BSC_ITEM_TaskA_vs_TaskB;
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_ITEM_FIR.mat'],'BSC*');
time = toc;
fprintf(['Merge BSC-ITEM matrices (After FIR task regression) :: Done in: ' num2str(time) 's \n']); 

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

sgtitle('Beta-series correlations: ITEM approach (After FIR task regression)')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(234),'parula') 
colormap(subplot(235),'parula') 
colormap(subplot(236),'parula') 

%% Sensitivity and Specificity
[TPR_BSC_ITEM, TNR_BSC_ITEM] = TPR_TNR(BSC_ITEM_TaskA_vs_TaskB_group_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_ITEM_FIR_RESULTS.mat'],'TPR*','TNR*');
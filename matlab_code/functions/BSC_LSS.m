function BSC_LSS(stat_path,sots_path,exp_folder,N,TR,model,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic

% Number of *.nii images per subject
dur = (length(dir([stat_path filesep exp_folder filesep 'funct_images'])) - 2)/N;
% Stimulus onset times (SOTs)
load(sots_path);
% Number of events
E = length(onsets{1,1})+length(onsets{1,2}); 

for subji = 1:N
    for j = 1:E
        %% Estimate
        matlabbatch{1}.spm.stats.fmri_spec.dir = ...
            {[stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_LSS']};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        for image = 1:dur
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{image,1} = ...
                [stat_path filesep exp_folder filesep 'funct_images' filesep 'Sub_' num2str(subji,'%.3d') '_Image_' num2str(image,'%.4d') '.nii,1'];
        end        

        all_onsets = [onsets{1,1}; onsets{1,2}];
        onset_event = all_onsets(j);
        all_events = (1:E)';
        other_events = all_events(find(all_events~=j));
        onset_other = all_onsets(other_events);
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Event';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = onset_event;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = durations{1,1};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Other';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = onset_other;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = durations{1,2};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
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
        spm_jobman('run',matlabbatch);
        clear matlabbatch
        %% Save betas
        beta_hdr = spm_vol([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_LSS' filesep 'beta_0001.nii']);
        betas_LSS(j,:) = spm_read_vols(beta_hdr)';
        clear beta_hdr
        %% Remove dir (to save drive space)
        rmdir([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_LSS'],'s');
    end
    %% Sort betas
    betas_TaskA = betas_LSS(1:length(onsets{1,1}),:);
    betas_TaskB = betas_LSS(length(onsets{1,1})+1:length(onsets{1,1})+length(onsets{1,2}),:);
    %% Beta-series correlations
    BSC_LSS_TaskA = atanh(corr(betas_TaskA));
    BSC_LSS_TaskB = atanh(corr(betas_TaskB));
    BSC_LSS_TaskA_vs_TaskB = BSC_LSS_TaskA - BSC_LSS_TaskB;
    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_LSS.mat'],'BSC*','betas_LSS');
    clear BSC* betas_LSS

    fprintf(['BSC-LSS :: Subject #' num2str(subji) ' :: Done in: ' num2str(toc) 's \n']);
end

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_LSS.mat']);
    BSC_LSS_TaskA_group(:,:,subji)  = BSC_LSS_TaskA;
    BSC_LSS_TaskB_group(:,:,subji)  = BSC_LSS_TaskB;
    BSC_LSS_TaskA_vs_TaskB_group(:,:,subji)  = BSC_LSS_TaskA_vs_TaskB;
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_LSS.mat'],'BSC*');
time = toc;
fprintf(['Merge BSC-LSS matrices :: Done in: ' num2str(time) 's \n']); 

%% Results: BSC LSS
[BSC_LSS_TaskA_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_LSS_TaskA_group,q_level);
[BSC_LSS_TaskB_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_LSS_TaskB_group,q_level);
[BSC_LSS_TaskA_vs_TaskB_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_LSS_TaskA_vs_TaskB_group,q_level);

figure(2)
subplot(231); imagesc(mean(BSC_LSS_TaskA_group,3)); title('LSS TaskA');  axis square; caxis(max_ax(mean(BSC_LSS_TaskA_group,3),1));
subplot(232); imagesc(mean(BSC_LSS_TaskB_group,3)); title('LSS TaskB');  axis square; caxis(max_ax(mean(BSC_LSS_TaskB_group,3),1));
subplot(233); imagesc(mean(BSC_LSS_TaskA_vs_TaskB_group,3)); title('LSS TaskAvsB');  axis square; caxis(max_ax(mean(BSC_LSS_TaskA_vs_TaskB_group,3),1));
subplot(234); imagesc(BSC_LSS_TaskA_group_FDR); title('LSS Task A FDR');  axis square; 
subplot(235); imagesc(BSC_LSS_TaskB_group_FDR); title('LSS Task B FDR');  axis square; 
subplot(236); imagesc(BSC_LSS_TaskA_vs_TaskB_group_FDR); title('LSS Task AvsB FDR');  axis square; 

sgtitle('Beta-series correlations: LSS approach')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')

%% Sensitivity and Specificity
[TPR_BSC_LSS, TNR_BSC_LSS] = TPR_TNR(BSC_LSS_TaskA_vs_TaskB_group_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_LSS_RESULTS.mat'],'TPR*','TNR*');
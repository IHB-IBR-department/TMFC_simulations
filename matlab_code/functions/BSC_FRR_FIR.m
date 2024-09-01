function BSC_FRR_FIR(stat_path,sots_path,exp_folder,N,TR,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, July, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic

% Number of *.nii images per subject
dur = length(dir([stat_path filesep exp_folder filesep 'funct_images' filesep 'Sub_001*']));

% Define unconvolved design matrix without upsampling
load(sots_path);
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.stats.fmri_design.dir = {[stat_path filesep exp_folder filesep 'GLM_no_upsampling' filesep 'Sub_001']};
matlabbatch{1}.spm.stats.fmri_design.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_design.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_design.timing.fmri_t = 1;
matlabbatch{1}.spm.stats.fmri_design.timing.fmri_t0 = 1;
matlabbatch{1}.spm.stats.fmri_design.sess.nscan = dur;
matlabbatch{1}.spm.stats.fmri_design.sess.cond(1).name = 'Task_A';
matlabbatch{1}.spm.stats.fmri_design.sess.cond(1).onset = onsets{1,1};
matlabbatch{1}.spm.stats.fmri_design.sess.cond(1).duration = TR/4;
matlabbatch{1}.spm.stats.fmri_design.sess.cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_design.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_design.sess.cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_design.sess.cond(2).name = 'Task_B';
matlabbatch{1}.spm.stats.fmri_design.sess.cond(2).onset = onsets{1,2};
matlabbatch{1}.spm.stats.fmri_design.sess.cond(2).duration = TR/4;
matlabbatch{1}.spm.stats.fmri_design.sess.cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_design.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_design.sess.cond(2).orth = 1;
matlabbatch{1}.spm.stats.fmri_design.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_design.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_design.sess.multi_reg = {''};
matlabbatch{1}.spm.stats.fmri_design.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_design.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_design.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_design.volt = 1;
matlabbatch{1}.spm.stats.fmri_design.global = 'None';
matlabbatch{1}.spm.stats.fmri_design.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_design.cvi = 'AR(1)';

spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);
spm_get_defaults('stats.resmem',1);
spm_get_defaults('stats.maxmem',2^33);
spm_jobman('run',matlabbatch);
clear matlabbatch

load([stat_path filesep exp_folder filesep 'GLM_no_upsampling' filesep 'Sub_001' filesep 'SPM.mat']);
task_design = [full(SPM.Sess.U(1).u(33:end)), full(SPM.Sess.U(2).u(33:end))]; 
rmdir([stat_path filesep exp_folder filesep 'GLM_no_upsampling'],'s');

half = ceil(dur/2);
last_8sec = ceil(8/TR);
first_half = task_design(1:half,:);
first_half(half-last_8sec:half,:) = 0;
second_half = task_design(half+1:dur,:);
second_half(length(second_half)-last_8sec:end,:) = 0;

design{1,1} = first_half;
design{1,2} = second_half;

task_design_no_last_8s = [first_half;second_half];

stimorder = task_design_no_last_8s(:,1) + 2*task_design_no_last_8s(:,2);
stimorder = nonzeros(stimorder)';

stimdur = durations{1,1};

%% Run GLMsingle

opt.wantglmdenoise = 0;  % Do not perform GLMdenoise
opt.wantlibrary = 0;     % Assume canonical HRF

for subji = 1:N
    VOI = load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);
    
    outputdir = [stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'GLMsingle_cHRF'];
    
    data = {VOI.VOI(1:half,:)',VOI.VOI(half+1:dur,:)'};
    
    out = 0;
    while(~out)
        try
            GLMestimatesingletrial(design,data,stimdur,TR,outputdir,opt);
            out = 1;
        end
    end

    load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'GLMsingle_cHRF' filesep 'TYPED_FITHRF_GLMDENOISE_RR.mat'], 'modelmd');
    betas = squeeze(modelmd);
    betas_TaskA = betas(:,stimorder == 1)';
    betas_TaskB = betas(:,stimorder == 2)';
    % Beta-series correlations
    BSC_FRR_TaskA = atanh(corr(betas_TaskA));
    BSC_FRR_TaskB = atanh(corr(betas_TaskB));
    BSC_FRR_TaskA_vs_TaskB = BSC_FRR_TaskA - BSC_FRR_TaskB;
    save([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_FRR.mat'],'BSC*');
    clear BSC* betas*
    fprintf(['BSC-FRR after FIR :: Subject #' num2str(subji) ' :: Done in: ' num2str(toc) 's \n']);
end

%% Merge matrices & Clear FRR folders
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'BSC_FRR.mat']);
    BSC_FRR_TaskA_group(:,:,subji)  = BSC_FRR_TaskA;
    BSC_FRR_TaskB_group(:,:,subji)  = BSC_FRR_TaskB;
    BSC_FRR_TaskA_vs_TaskB_group(:,:,subji)  = BSC_FRR_TaskA_vs_TaskB;
    rmdir([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'GLMsingle_cHRF'],'s');
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_FRR_FIR.mat'],'BSC*');
time = toc;
fprintf(['Merge BSC-FRR matrices (after FIR) :: Done in: ' num2str(time) 's \n']); 

%% Results: BSC FRR (fractional ridge regression)
[BSC_FRR_TaskA_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_FRR_TaskA_group,q_level);
[BSC_FRR_TaskB_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_FRR_TaskB_group,q_level);
[BSC_FRR_TaskA_vs_TaskB_group_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(BSC_FRR_TaskA_vs_TaskB_group,q_level);
gm_BSC_FRR_TaskA_vs_TaskB = mean(BSC_FRR_TaskA_vs_TaskB_group,3);
gm_BSC_FRR_TaskA_vs_TaskB(1:1+size(gm_BSC_FRR_TaskA_vs_TaskB,1):end) = 0;
BSC_FRR_TaskA_vs_TaskB_group_FDR(1:1+size(BSC_FRR_TaskA_vs_TaskB_group_FDR,1):end) = 0;

figure
try
    sgtitle('Beta-series correlations: FRR approach (after FIR task regression)')
catch
    suptitle('Beta-series correlations: FRR approach (after FIR task regression)')
end
subplot(231); imagesc(mean(BSC_FRR_TaskA_group,3)); title('FRR TaskA');       axis square; caxis(max_ax(mean(BSC_FRR_TaskA_group,3),1));
subplot(232); imagesc(mean(BSC_FRR_TaskB_group,3)); title('FRR TaskB');       axis square; caxis(max_ax(mean(BSC_FRR_TaskB_group,3),1));
subplot(233); imagesc(gm_BSC_FRR_TaskA_vs_TaskB);   title('FRR TaskAvsB');    axis square; caxis(max_ax(gm_BSC_FRR_TaskA_vs_TaskB,1));
subplot(234); imagesc(BSC_FRR_TaskA_group_FDR);     title('FRR Task A FDR');  axis square; 
subplot(235); imagesc(BSC_FRR_TaskB_group_FDR);     title('FRR Task B FDR');  axis square; 
subplot(236); imagesc(BSC_FRR_TaskA_vs_TaskB_group_FDR); title('FRR Task AvsB FDR');  axis square; 
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(234),'parula') 
colormap(subplot(235),'parula') 
colormap(subplot(236),'parula') 

%% Sensitivity and Specificity
[TPR_BSC_FRR, TNR_BSC_FRR] = TPR_TNR(BSC_FRR_TaskA_vs_TaskB_group_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'BSC_FRR_FIR_RESULTS.mat'],'TPR*','TNR*');




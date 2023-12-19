function parallel_extract_VOI(stat_path,exp_folder,N,N_ROIs)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

for subji = 1:N
    for j = 1:N_ROIs
        matlabbatch{1}.spm.util.voi.spmmat = ...
            {[stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'SPM.mat']};
        matlabbatch{1}.spm.util.voi.adjust = 1;
        matlabbatch{1}.spm.util.voi.session = 1;
        matlabbatch{1}.spm.util.voi.name = ['ROI_' num2str(j,'%03.f')];
        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = ...
            {[stat_path filesep exp_folder filesep 'ROI_masks' filesep 'ROI_' num2str(j,'%03.f') '.nii,1']};
        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.1;
        matlabbatch{1}.spm.util.voi.expression = 'i1';    
        
        batch{j} = matlabbatch;
        clear matlabbatch
    end
    
    parfor j = 1:N_ROIs
        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_get_defaults('cmdline',true);
        spm_get_defaults('stats.resmem',1);
        spm_get_defaults('stats.maxmem',2^34);
        spm_jobman('run',batch{j});
    end

    clear batch
end

%% Merge VOIs
for subji = 1:N
    for j = 1:N_ROIs
        load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'VOI_ROI_' num2str(j,'%03.f') '_1.mat']);
        VOI(:,j) = Y;
        clear Y
    end    
    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat'], 'VOI');
         
    clear VOI
end
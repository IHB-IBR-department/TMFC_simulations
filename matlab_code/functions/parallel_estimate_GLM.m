function parallel_estimate_GLM(stat_path,sots_path,exp_folder,N,TR,model)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

% Number of *.nii images per subject
dur = (length(dir([stat_path filesep exp_folder filesep 'funct_images'])) - 2)/N;

spm('defaults','fmri');
spm_jobman('initcfg');

for subji = 1:N    
    matlabbatch{1}.spm.stats.fmri_spec.dir = ...
        {[stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d')]};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    for image = 1:dur
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{image,1} = ...
            [stat_path filesep exp_folder filesep 'funct_images' filesep 'Sub_' num2str(subji,'%03.f') '_Image_' num2str(image,'%04.f') '.nii,1'];
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
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Effects of interest';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(2);
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;
    
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
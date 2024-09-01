function generate_funct_images_FIR(stat_path,sim_path,exp_folder,SF,SNR,N,N_ROIs,dummy)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);
    BOLD(:,:,subji) = VOI;
end

mkdir([stat_path filesep exp_folder filesep 'funct_images_after_FIR_task_regression'])
BOLD = BOLD + 100.*ones(size(BOLD));

% Generate .nii images for SPM 
for i = 1:N
    for j = 1:length(BOLD)
        nii_image = make_nii(BOLD(j,:,i)');
        save_nii(nii_image,[stat_path filesep exp_folder filesep 'funct_images_after_FIR_task_regression' filesep 'Sub_' num2str(i,'%03.f') '_Image_' num2str(j,'%04.f') '.nii']); 
        clear nii_image
    end
end

time = toc;
fprintf(['Generate images :: ' exp_folder ' :: Done in ' num2str(time)  '\n']);
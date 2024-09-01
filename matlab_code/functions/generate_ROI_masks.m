function generate_ROI_masks(stat_path,exp_folder,N_ROIs)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
mask = eye(N_ROIs);

% Make folder for ROI masks
mkdir([stat_path filesep exp_folder filesep 'ROI_masks']);

for i = 1:N_ROIs
        nii_image = make_nii(mask(:,i));
        save_nii(nii_image,[stat_path filesep exp_folder filesep 'ROI_masks' filesep ...
            'ROI_' num2str(i,'%03.f') '.nii']); 
        clear nii_image
end

time = toc;
fprintf(['Generate ROI masks :: ' exp_folder ' :: Done in ' num2str(time)  '\n']);

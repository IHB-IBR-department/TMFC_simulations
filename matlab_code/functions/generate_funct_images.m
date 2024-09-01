function generate_funct_images(stat_path,sim_path,exp_folder,SF,SNR,N,N_ROIs,dummy)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
load(sim_path)

% Remove first dummy dynamics
oscill(1:dummy,:,:) = [];
coact(1:dummy,:) = [];

% Simulation duration
dur = size(oscill,1);

% Normalize oscillations
oscill = zscore(oscill) + 100.*ones(size(oscill));

% Scale coactivations
coact = zscore(coact)./SF;
coact(isinf(coact)) = 0;
coact = repmat(coact,1,1,N);

% Additive white gaussian noise (AWGN)
noise = randn(dur,N_ROIs,N)./SNR;
noise(isinf(noise)) = 0;

% BOLD = WC-oscillations + coactivations + noise
BOLD = oscill(1:dur,1:N_ROIs,1:N) + coact(1:dur,1:N_ROIs,:) + noise;

% Make folder for .nii images
mkdir([stat_path filesep exp_folder filesep 'funct_images']);

% Generate .nii images for SPM 
for i = 1:N
    for j = 1:dur
        nii_image = make_nii(BOLD(j,:,i)');
        save_nii(nii_image,[stat_path filesep exp_folder filesep 'funct_images' filesep 'Sub_' num2str(i,'%03.f') '_Image_' num2str(j,'%04.f') '.nii']); 
        clear nii_image
    end
end

time = toc;
fprintf(['Generate images :: ' exp_folder ' :: Done in ' num2str(time)  '\n']);
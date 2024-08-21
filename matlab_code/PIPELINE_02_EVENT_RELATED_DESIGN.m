% ========================================================================
% TMFC analysis for event-related design
% Requires SPM12 (v7771)
%
% Beta-Series Correlations based on Fractional Ridge Regression (BSC-FRR)
% requires GLMsingle toolbox (https://github.com/cvnlab/GLMsingle)
%
% ========================================================================
% Ruslan Masharipov, August 7, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

%% Define parameters
close all
clear

% Set path for stat folder 
stat_path = 'C:\TMFC_simulations\experiments\02_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS]';

% Set path for simulated BOLD time series *.mat file
sim_path = 'C:\TMFC_simulations\simulated_BOLD_time_series\SIM_BOLD_02_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].mat';

% Set path for task design *.mat file (stimulus onset times, SOTs)
% Simular to the multiple condition *.mat file used in SPM 12
sots_path = 'C:\TMFC_simulations\task_designs\02_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].mat';

% Symmetric ground truth matrix
load('C:\TMFC_simulations\matlab_code\ground_truth_symm_matrix.mat');

% BOLD-signal = BOLD(Wilson-Cowan oscillations) + BOLD(Co-activations) + White Gaussian Noise

% Co-activations
% Scaling Factor (SF):
% SF = SD_oscill/SD_coact
% Set SF = 0 for no co-activations
SF = 1;

% Additive white gaussian noise
% Signal-to-noise ratio (SNR):
% SNR = SD_signal/SD_noise
% Set SNR = 0 for no noise
SNR = 0.4;

% Short-term synaptic plasticity (STP) delay
% Delay between stimulus onset and change in synaptic weights, [s]
STP_delay = 0.2;

% Sample size
N = 100;

% Number of ROIs
N_ROIs = 100;

% Remove first X dummy scans
dummy = 3;

% Repetition time (TR), [s]
TR = 2;

% FIR parameters
FIR_window = 24; % length of FIR window [s]
FIR_bins = 24;   % number of time bins

% Autocorrelation modeling
% 'AR(1)' - default model in SPM
% 'FAST' - model for faster sampling rates (shorter TRs)
model = 'AR(1)';

% FDR correction
q_level = 0.001/2; % Corrected for two-sided comparison

% Experiment folder
exp_folder = ['SF_[' num2str(SF,'%.2f') ']_SNR_[' num2str(SNR,'%.2f') ']_STP_[' num2str(STP_delay,'%.2f') ']_'  model];

% Add STP delay
load(sots_path);
onsets{1,1} = onsets{1,1} - STP_delay;
onsets{1,2} = onsets{1,2} - STP_delay;
sots_path = sots_path(1:end-4);
sots_path = [sots_path '_[' num2str(data.STP_delay,'%.2f') 's_STP].mat'];
save(sots_path,'activations','onsets','durations','names','rest_matrix','task_matrices');

%% Generate .nii functional images for SPM
generate_funct_images(stat_path,sim_path,exp_folder,SF,SNR,N,N_ROIs,dummy)

%% Generate .nii ROI binary masks for SPM VOI extraction 
generate_ROI_masks(stat_path,exp_folder,N_ROIs)

%% Estimate GLM
tic
parallel_estimate_GLM(stat_path,sots_path,exp_folder,N,TR,model)
fprintf(['Estimate GLM :: Done in: ' num2str(toc) 's \n']);

%% Estimate FIR GLM
tic
parallel_estimate_FIR_GLM(stat_path,sots_path,exp_folder,N,TR,model,FIR_window,FIR_bins);
fprintf(['Estimate FIR GLM :: Done in: ' num2str(toc) 's \n']);

%% VOI time-series extraction of adjusted data
tic
parallel_extract_VOI(stat_path,exp_folder,N,N_ROIs)
fprintf(['Extract VOI :: Done in: ' num2str(toc) 's \n']);

%% VOI time-series extraction of adjusted data (after FIR)
tic
parallel_extract_FIR_VOI(stat_path,exp_folder,N,N_ROIs)
fprintf(['Extract FIR VOI :: Done in: ' num2str(toc) 's \n']);

%% Calculate PPIs
tic
parallel_calculate_PPIs(stat_path,exp_folder,N,N_ROIs)
fprintf(['Calculate FIR PPI :: Done in: ' num2str(toc) 's \n']);

%% Calculate PPIs (after FIR)
tic
parallel_calculate_FIR_PPIs(stat_path,exp_folder,N,N_ROIs)
fprintf(['Calculate FIR PPI :: Done in: ' num2str(toc) 's \n']);

%% sPPI and gPPI with deconvolution
sPPI_and_gPPI_with_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% sPPI and gPPI with deconvolution (after FIR)
sPPI_and_gPPI_with_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% sPPI and gPPI without deconvolution
sPPI_and_gPPI_without_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% sPPI and gPPI without deconvolution (after FIR)
sPPI_and_gPPI_without_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% cPPI with deconvolution
cPPI_with_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% cPPI with deconvolution (after FIR)
cPPI_with_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% cPPI without deconvolution
cPPI_without_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% cPPI without deconvolution (after FIR)
cPPI_without_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% Generate .nii residual images (after FIR)
generate_funct_images_FIR(stat_path,sim_path,exp_folder,SF,SNR,N,N_ROIs,dummy)

%% Beta-Series Correlations: Least Squares All
BSC_LSA(stat_path,sots_path,exp_folder,N,TR,model,q_level,ground_truth)

%% Beta-Series Correlations: Least Squares All (after FIR)
BSC_LSA_FIR(stat_path,sots_path,exp_folder,N,TR,model,q_level,ground_truth)

%% Beta-Series Correlations: Least Squares Separate
parallel_BSC_LSS(stat_path,sots_path,exp_folder,N,TR,model,q_level,ground_truth);

%% Beta-Series Correlations: Least Squares Separate (after FIR)
parallel_BSC_LSS_FIR(stat_path,sots_path,exp_folder,N,TR,model,q_level,ground_truth);

%% Beta-Series Correlations: Fractional Ridge Regression 
BSC_FRR(stat_path,sots_path,exp_folder,N,TR,q_level,ground_truth)

%% Beta-Series Correlations: Fractional Ridge Regression (after FIR)
BSC_FRR_FIR(stat_path,sots_path,exp_folder,N,TR,q_level,ground_truth)

%% Task-state FC and Background FC
TSFC_BGFC(stat_path,exp_folder,N,N_ROIs)

%% TMFC results 
% (Plot all figures for a given q-level)
close all
q_level = 0.001/2;
event_related_results(stat_path,exp_folder,q_level,ground_truth)

%% TMFC results after FIR task regression 
% (Plot all figures for a given q-level)
q_level = 0.05/2;
event_related_results_FIR(stat_path,exp_folder,q_level,ground_truth)






%% ===============[No centering before PPI calculation]====================
% Open spm_peb_ppi and comment out line 439: PSY = spm_detrend(PSY);
% Run PPI estimation
parallel_calculate_FIR_PPIs_no_centering(stat_path,exp_folder,N,N_ROIs)
% Uncomment PSY = spm_detrend(PSY);

%% sPPI and gPPI with deconvolution + no cetering (after FIR)
sPPI_and_gPPI_with_deconv_no_centering_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

%% sPPI and gPPI without deconvolution + no cetering (after FIR)
sPPI_and_gPPI_without_deconv_no_centering_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)


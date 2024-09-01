% ========================================================================
% TMFC (gPPI) analysis for event-related design with 
% asymmetric synaptic weight matrices
% Requires SPM12 (v7771)
% ========================================================================
% Ruslan Masharipov, July, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

%% Define parameters
close all
clear

% Set path for stat folder 
stat_path = 'C:\TMFC_simulations\experiments\27_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[200_TRIALS]_ASYMMETRIC';

% Set path for Wilson-Cowan (WC) simulations .mat file
sim_path = 'C:\TMFC_simulations\simulated_BOLD_time_series\files\SIM_BOLD_27_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[200_TRIALS]_ASYMMETRIC.mat';

% Set path for task design *.mat file (stimulus onset times, SOTs)
% Simular to the multiple condition *.mat file used in SPM 12
sots_path = 'C:\TMFC_simulations\task_designs\files\27_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[200_TRIALS]_ASYMMETRIC.mat';

% Asymmetric ground truth matrix
load('C:\TMFC_simulations\matlab_code\ground_truth_asymm_matrix.mat');

% BOLD-signal = BOLD(Wilson-Cowan oscillations) + BOLD(Co-activations) + White Gaussian Noise

% Co-activations
% Scaling Factor (SF):
% SF = SD_oscill/SD_coact
% Set SF = 0 for no co-activations
SF = 0;

% Additive white gaussian noise
% Signal-to-noise ratio (SNR):
% SNR = SD_signal/SD_noise
% Set SNR = 0 for no noise
SNR = 5;

% Short-term synaptic plasticity (STP) delay
% Delay between stimulus onset and change in synaptic weights, [s]
STP_delay = 0.2;

% Sample size
N = 100;

% Number of ROIs
N_ROIs = 100;

% Remove first X dummy scanns
dummy = 3;

% Repetition Time (TR), [s]
TR = 2;

% FIR parameters
FIR_window = 24; % length of FIR window [s]
FIR_bins = 24;   % number of time bins

% Autocorrelation modeling
% 'AR(1)' - default model in SPM
% 'FAST' - model for faster sampling rates (shorter TRs)
model = 'AR(1)';

% FDR correction
q_level = 0.001/2;

% Experiment folder
exp_folder = ['SF_[' num2str(SF,'%.2f') ']_SNR_[' num2str(SNR,'%.2f') ']_STP_[' num2str(STP_delay,'%.2f') ']_'  model];

% STP delay
load(sots_path);
onsets{1,1} = onsets{1,1} - STP_delay;
onsets{1,2} = onsets{1,2} - STP_delay;
sots_path = sots_path(1:end-4);
sots_path = [sots_path '_[' num2str(STP_delay,'%.2f') 's_STP].mat'];
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
parallel_calculate_PPIs (stat_path,exp_folder,N,N_ROIs)
fprintf(['Calculate PPI :: Done in: ' num2str(toc) 's \n']);

%% Calculate PPIs (after FIR)
tic
parallel_calculate_FIR_PPIs(stat_path,exp_folder,N,N_ROIs)
fprintf(['Calculate FIR PPI :: Done in: ' num2str(toc) 's \n']);

%% sPPI and gPPI with deconvolution
sPPI_and_gPPI_with_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% sPPI and gPPI without deconvolution
sPPI_and_gPPI_without_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% sPPI and gPPI with deconvolution (after FIR)
sPPI_and_gPPI_with_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% sPPI and gPPI without deconvolution (after FIR)
sPPI_and_gPPI_without_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% Correlations between ground truth and gPPI
% load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_with_Deconv_FIR.mat'])    % If scaling factor, SF =/= 0 
% load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv_FIR.mat']) % If scaling factor, SF =/= 0 
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_with_Deconv.mat'])          % If scaling factor, SF = 0 (no co-activations)
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv.mat'])       % If scaling factor, SF = 0 (no co-activations)


mean_gPPI_WD  = mean(gPPI_WD_TaskA_vs_TaskB_asymm,3);
mean_gPPI_WoD  = mean(gPPI_WoD_TaskA_vs_TaskB_asymm,3);
mean_gPPI_WD(1:1+size(mean_gPPI_WD,1):end) = 0;
mean_gPPI_WoD(1:1+size(mean_gPPI_WoD,1):end) = 0;

figure
subplot(131); imagesc(ground_truth_asymm);  axis square; caxis(max_ax(ground_truth_asymm,1)); 
subplot(132); imagesc(mean_gPPI_WoD);       axis square; caxis(max_ax(mean_gPPI_WoD,1)); 
subplot(133); imagesc(mean_gPPI_WD);        axis square; caxis(max_ax(mean_gPPI_WD,1)); 

set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 

truth = [lower_triangle(ground_truth_asymm), upper_triangle(ground_truth_asymm)]; 
gPPI_WD =  [lower_triangle(mean_gPPI_WD), upper_triangle(mean_gPPI_WD)]; 
gPPI_WoD =  [lower_triangle(mean_gPPI_WoD), upper_triangle(mean_gPPI_WoD)]; 

truth_vs_gPPI_WD = corr(truth',gPPI_WD');
truth_vs_gPPI_WoD = corr(truth',gPPI_WoD');

% Correct Sign Rate (CSR)
gPPI_WD_CSR = sum((gPPI_WD.*truth)>0)/nnz(truth)*100;
gPPI_WoD_CSR = sum((gPPI_WoD.*truth)>0)/nnz(truth)*100;


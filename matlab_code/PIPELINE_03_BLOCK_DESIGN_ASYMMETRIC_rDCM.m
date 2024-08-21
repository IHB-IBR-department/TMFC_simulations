% ========================================================================
% TMFC (gPPI) and TMEC (rDCM) analyses for block design with asymmetric
% synaptic weight matrices
% Requires SPM12 (v7771) and TAPAS toolbox (rDCM)
% ========================================================================
% Ruslan Masharipov, July, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

%% Define parameters
close all
clear

% Set path for stat folder 
stat_path = 'C:\TMFC_simulations\experiments\25_BLOCK_[2s_TR]_[20s_DUR]_[20_BLOCKS]_ASYMMETRIC';

% Set path for Wilson-Cowan (WC) simulations .mat file
sim_path = 'C:\TMFC_simulations\simulated_BOLD_time_series\SIM_BOLD_25_BLOCK_[2s_TR]_[20s_DUR]_[20_BLOCKS]_ASYMMETRIC.mat';

% Set path for task design *.mat file (stimulus onset times, SOTs)
% Simular to the multiple condition *.mat file used in SPM 12
sots_path = 'C:\TMFC_simulations\task_designs\25_BLOCK_[2s_TR]_[20s_DUR]_[20_BLOCKS]_ASYMMETRIC.mat';

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
FIR_window = 32; % length of FIR window [s]
FIR_bins = 32;   % number of time bins

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
sots_path = [sots_path '_[' num2str(data.STP_delay,'%.2f') 's_STP].mat'];
save(sots_path,'activations','onsets','durations','names','rest_matrix','task_matrices');

%% Generate .nii functional images for SPM estimate
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
parallel_calculate_FIR_PPIs (stat_path,exp_folder,N,N_ROIs)
fprintf(['Calculate FIR PPI :: Done in: ' num2str(toc) 's \n']);

%% sPPI and gPPI with deconvolution
sPPI_and_gPPI_with_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% sPPI and gPPI without deconvolution
sPPI_and_gPPI_without_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% sPPI and gPPI with deconvolution (after FIR)
sPPI_and_gPPI_with_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% sPPI and gPPI without deconvolution (after FIR)
sPPI_and_gPPI_without_deconv_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth_asymm)

%% Prepare block time series for rDCM
load(sots_path);
%load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' %num2str(1,'%.3d') filesep  'VOI_ALL.mat']);   % If scaling factor, SF =/= 0 
load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(1,'%.3d') filesep  'VOI_ALL.mat']);         % If scaling factor, SF = 0 (no co-activations)

length_N = size(VOI,1);
imagei = 1:length_N;
imaget = (imagei)*TR;

% TaskA
for i = 1:length(onsets{1,1})
    TaskA_block_time(i,:) = (imaget>onsets{1,1}(i)+6 & imaget<=onsets{1,1}(i)+durations{1,1});
end
TaskA_block_time = sum(TaskA_block_time,1);

% TaskB
for i = 1:length(onsets{1,2})
    TaskB_block_time(i,:) = (imaget>onsets{1,2}(i)+6 & imaget<=onsets{1,2}(i)+durations{1,2});
end
TaskB_block_time = sum(TaskB_block_time,1);

% Rest
all_onsets = sort([onsets{1,1}; onsets{1,2}]);
dur = durations{1,1};
for i = 1:length(all_onsets)
    Rest_block_time(i,:) = (imaget>all_onsets(i)+dur+6 & imaget<=all_onsets(i)+2*dur);
end
Rest_block_time = sum(Rest_block_time,1);

%% rDCM
tic
for subji = 1:N
    %load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' %num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);    % If scaling factor, SF =/= 0 
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);          % If scaling factor, SF = 0 (no co-activations)

    TS_TaskA.y = VOI(TaskA_block_time==1,:);
    TS_TaskB.y = VOI(TaskB_block_time==1,:);
    TS_Rest.y  = VOI(Rest_block_time==1,:);
    TS_TaskA.dt = TR;
    TS_TaskB.dt = TR;
    TS_Rest.dt  = TR;

    DCM_TaskA = tapas_rdcm_model_specification(TS_TaskA,[],[]); 
    DCM_TaskB = tapas_rdcm_model_specification(TS_TaskB,[],[]); 
    DCM_Rest  = tapas_rdcm_model_specification(TS_Rest,[],[]); 
    [output_TaskA options] = tapas_rdcm_estimate(DCM_TaskA,'r',[],1);
    rDCM_TaskA(:,:,subji) = output_TaskA.Ep.A;
    clear output_TaskA
    [output_TaskB options] = tapas_rdcm_estimate(DCM_TaskB,'r',[],1);
    rDCM_TaskB(:,:,subji) = output_TaskB.Ep.A;
    clear output_TaskB 
    [output_Rest options] = tapas_rdcm_estimate(DCM_Rest,'r',[],1);
    rDCM_Rest(:,:,subji) = output_Rest.Ep.A;
    clear output_TaskB 

    rDCM_TaskA_vs_TaskB(:,:,subji) = rDCM_TaskA(:,:,subji) - rDCM_TaskB(:,:,subji);
    
    clear VOI TS_TaskA TS_TaskB TS_Rest DCM_TaskA DCM_TaskB DCM_Rest
    
    fprintf(['rDCM :: Subject #' num2str(subji) ': Done :: Time: ' num2str(toc) 's \n'])
end

mkdir([stat_path filesep exp_folder filesep 'group_stat']);
%save([stat_path filesep exp_folder filesep 'group_stat' filesep 'rDCM_FIR.mat'],'rDCM_TaskA','rDCM_TaskB','rDCM_Rest','rDCM_TaskA_vs_TaskB');     % If scaling factor, SF =/= 0 
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'rDCM.mat'],'rDCM_TaskA','rDCM_TaskB','rDCM_Rest','rDCM_TaskA_vs_TaskB');          % If scaling factor, SF = 0 (no co-activations)

%% rDCM matrices
%load([stat_path filesep exp_folder filesep 'group_stat' filesep 'rDCM_FIR.mat']) % If scaling factor, SF =/= 0 
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'rDCM.mat'])      % If scaling factor, SF = 0 (no co-activations)
Rest  = mean(rDCM_Rest,3);
TaskA = mean(rDCM_TaskA,3);
TaskB = mean(rDCM_TaskB,3);
TaskA_vs_TaskB = mean(rDCM_TaskA_vs_TaskB,3);

Rest(1:1+size(Rest,1):end) = 0;
TaskA(1:1+size(TaskA,1):end) = 0;
TaskB(1:1+size(TaskB,1):end) = 0;
TaskA_vs_TaskB(1:1+size(TaskA_vs_TaskB,1):end) = 0;

[rDCM_Rest_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(rDCM_Rest,q_level);
[rDCM_TaskA_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(rDCM_TaskA,q_level);
[rDCM_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(rDCM_TaskB,q_level);
[rDCM_TaskA_vs_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(rDCM_TaskA_vs_TaskB,q_level);
rDCM_Rest_FDR(1:1+size(rDCM_Rest_FDR,1):end) = 0;
rDCM_TaskA_FDR(1:1+size(rDCM_TaskA_FDR,1):end) = 0;
rDCM_TaskB_FDR(1:1+size(rDCM_TaskB_FDR,1):end) = 0;
rDCM_TaskA_vs_TaskB_FDR(1:1+size(rDCM_TaskA_vs_TaskB_FDR,1):end) = 0;

figure
subplot(241); imagesc(Rest); title('Rest');  axis square; 
subplot(242); imagesc(TaskA); title('Task A');  axis square; 
subplot(243); imagesc(TaskB); title('Task B');  axis square; 
subplot(244); imagesc(TaskA_vs_TaskB); title('Task AvsB');  axis square; 
subplot(245); imagesc(rDCM_Rest_FDR); title('Rest FDR');  axis square; 
subplot(246); imagesc(rDCM_TaskA_FDR); title('Task A FDR');  axis square; 
subplot(247); imagesc(rDCM_TaskB_FDR); title('Task B FDR');  axis square; 
subplot(248); imagesc(rDCM_TaskA_vs_TaskB_FDR); title('Task AvsB FDR');  axis square; 

sgtitle('rDCM')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 
colormap(subplot(245),'parula')
colormap(subplot(246),'parula')
colormap(subplot(247),'parula')
colormap(subplot(248),'parula')

%% Correlations between ground truth, rDCM and gPPI
% load([stat_path filesep exp_folder filesep 'group_stat' filesep 'rDCM_FIR.mat'])
% load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_with_Deconv_FIR.mat'])
% load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv_FIR.mat'])
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'rDCM.mat'])
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_with_Deconv.mat'])
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv.mat'])

mean_rDCM = mean(rDCM_TaskA_vs_TaskB,3);
mean_gPPI_WD  = mean(gPPI_WD_TaskA_vs_TaskB_asymm,3);
mean_gPPI_WoD  = mean(gPPI_WoD_TaskA_vs_TaskB_asymm,3);
mean_gPPI_WD(1:1+size(mean_gPPI_WD,1):end) = 0;
mean_gPPI_WoD(1:1+size(mean_gPPI_WoD,1):end) = 0;

figure
subplot(221); imagesc(ground_truth_asymm);  axis square; caxis(max_ax(ground_truth_asymm,1)); 
subplot(222); imagesc(mean_rDCM);           axis square; caxis(max_ax(mean_rDCM,1)); 
subplot(223); imagesc(mean_gPPI_WoD);       axis square; caxis(max_ax(mean_gPPI_WoD,1)); 
subplot(224); imagesc(mean_gPPI_WD);        axis square; caxis(max_ax(mean_gPPI_WD,1)); 

set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 

truth = [lower_triangle(ground_truth_asymm), upper_triangle(ground_truth_asymm)]; 
rDCM =  [lower_triangle(mean_rDCM), upper_triangle(mean_rDCM)]; 
gPPI_WD =  [lower_triangle(mean_gPPI_WD), upper_triangle(mean_gPPI_WD)]; 
gPPI_WoD =  [lower_triangle(mean_gPPI_WoD), upper_triangle(mean_gPPI_WoD)]; 

truth_vs_rDCM = corr(truth',rDCM');
truth_vs_gPPI_WD = corr(truth',gPPI_WD');
truth_vs_gPPI_WoD = corr(truth',gPPI_WoD');
rDCM_vs_gPPI_WD = corr(rDCM',gPPI_WD');
rDCM_vs_gPPI_WoD = corr(rDCM',gPPI_WoD'); 

% Correct Sign Rate (CSR)
rDCM_CSR = sum((rDCM.*truth)>0)/nnz(truth)*100;
gPPI_WD_CSR = sum((gPPI_WD.*truth)>0)/nnz(truth)*100;
gPPI_WoD_CSR = sum((gPPI_WoD.*truth)>0)/nnz(truth)*100;



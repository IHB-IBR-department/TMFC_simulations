function sPPI_and_gPPI_without_deconv_no_centering_FIR(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, July, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'No_centering_PPI_FIR_ALL.mat']);
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'SPM.mat']);
    Psy_TaskA = SPM.xX.X(:,1);
    Psy_TaskB = SPM.xX.X(:,2);
    Psy_TaskA_vs_TaskB = (SPM.xX.X(:,1) - SPM.xX.X(:,2));
    PPI_TaskA_vs_TaskB = [];
    PPI_TaskA = [];
    PPI_TaskB = [];

    for j = 1:N_ROIs    
        PPI_TaskA_vs_TaskB = VOI(:,j).*Psy_TaskA_vs_TaskB;
        PPI_TaskA =          VOI(:,j).*Psy_TaskA;
        PPI_TaskB =          VOI(:,j).*Psy_TaskB; 
        
        sPPI_design = [PPI_TaskA_vs_TaskB VOI(:,j) Psy_TaskA_vs_TaskB ones(size(Psy_TaskA_vs_TaskB))];
        gPPI_design = [PPI_TaskA PPI_TaskB VOI(:,j) Psy_TaskA Psy_TaskB ones(size(Psy_TaskA))];
       
       
        for k = 1:N_ROIs
            if k==j
                beta_coeff_sPPI(:,j,k) = NaN(4,1);
                con_TaskA_vs_TaskB_sPPI(j,k) = NaN;
                task_indep_sPPI(j,k) = NaN;
                
                beta_coeff_gPPI(:,j,k) = NaN(6,1); 
                con_TaskA(j,k) = NaN;
                con_TaskB(j,k) = NaN;
                con_TaskA_vs_TaskB_gPPI(j,k) = NaN;
                task_indep_gPPI(j,k) = NaN;
            else
                [b_sppi,bint,res] = regress(zscore(VOI(:,k)),sPPI_design);
                [b_gppi,bint,res] = regress(zscore(VOI(:,k)),gPPI_design);
                
                beta_coeff_sPPI(:,j,k) = b_sppi;
                con_TaskA_vs_TaskB_sPPI(j,k) = [1 0 0 0]*b_sppi;
                task_indep_sPPI(j,k)   = [0 1 0 0]*b_sppi;
                
                beta_coeff_gPPI(:,j,k) = b_gppi;
                con_TaskA(j,k) = [1 0 0 0 0 0]*b_gppi;
                con_TaskB(j,k)   = [0 1 0 0 0 0]*b_gppi;
                con_TaskA_vs_TaskB_gPPI(j,k) = [1 -1 0 0 0 0]*b_gppi;
                task_indep_gPPI(j,k)   = [0 0 1 0 0 0]*b_gppi;
            end
        end
    end
   
    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'sPPI_and_gPPI_without_Deconv_No_centering_FIR.mat'],...
          'beta_*', 'con_*', 'task_indep_*');
    time(subji)=toc;
    fprintf(['sPPI and gPPI without deconvolution (No centering) :: Subject #' num2str(subji) ' :: Done in: ' num2str(time(subji)) 's \n']);     
    
    clearvars -except stat_path exp_folder N N_ROIs q_level ground_truth
end   

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'sPPI_and_gPPI_without_Deconv_No_centering_FIR.mat']);
    sPPI_WoD_TaskA_vs_TaskB_asymm(:,:,subji) = con_TaskA_vs_TaskB_sPPI;
    sPPI_WoD_TaskA_vs_TaskB_symm(:,:,subji) = (con_TaskA_vs_TaskB_sPPI + con_TaskA_vs_TaskB_sPPI')/2;
    sPPI_WoD_TaskIndep_asymm(:,:,subji) = task_indep_sPPI;
    sPPI_WoD_TaskIndep_symm(:,:,subji) = (task_indep_sPPI + task_indep_sPPI')/2;
    
    gPPI_WoD_TaskA_asymm(:,:,subji) = con_TaskA;
    gPPI_WoD_TaskA_symm(:,:,subji) = (con_TaskA + con_TaskA')/2;
    gPPI_WoD_TaskB_asymm(:,:,subji) = con_TaskB;
    gPPI_WoD_TaskB_symm(:,:,subji) = (con_TaskB + con_TaskB')/2;
    gPPI_WoD_TaskA_vs_TaskB_asymm(:,:,subji) = con_TaskA_vs_TaskB_gPPI;
    gPPI_WoD_TaskA_vs_TaskB_symm(:,:,subji) = (con_TaskA_vs_TaskB_gPPI + con_TaskA_vs_TaskB_gPPI')/2;
    gPPI_WoD_TaskIndep_asymm(:,:,subji) = task_indep_gPPI;
    gPPI_WoD_TaskIndep_symm(:,:,subji) = (task_indep_gPPI + task_indep_gPPI')/2;
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv_No_centering_FIR.mat'],'sPPI*','gPPI*');
time = toc;
fprintf(['Merge sPPI_and_gPPI_without_Deconv (No centering) matrices :: Done in: ' num2str(time) 's \n']); 

%% Results: sPPI without Deconv
[sPPI_WoD_TaskA_vs_TaskB_asymm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(sPPI_WoD_TaskA_vs_TaskB_asymm,q_level);
[sPPI_WoD_TaskA_vs_TaskB_symm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(sPPI_WoD_TaskA_vs_TaskB_symm,q_level);
gm_sPPI_WoD_TaskIndep_asymm = mean(sPPI_WoD_TaskIndep_asymm,3);
gm_sPPI_WoD_TaskIndep_symm  = mean(sPPI_WoD_TaskIndep_symm,3);
gm_sPPI_WoD_TaskA_vs_TaskB_asymm = mean(sPPI_WoD_TaskA_vs_TaskB_asymm,3);
gm_sPPI_WoD_TaskA_vs_TaskB_symm = mean(sPPI_WoD_TaskA_vs_TaskB_symm,3);
gm_sPPI_WoD_TaskIndep_asymm(1:1+size(sPPI_WoD_TaskIndep_asymm,1):end) = 1;
gm_sPPI_WoD_TaskIndep_symm(1:1+size(sPPI_WoD_TaskIndep_symm,1):end) = 1;
gm_sPPI_WoD_TaskA_vs_TaskB_asymm(1:1+size(gm_sPPI_WoD_TaskA_vs_TaskB_asymm,1):end) = 0;
gm_sPPI_WoD_TaskA_vs_TaskB_symm(1:1+size(gm_sPPI_WoD_TaskA_vs_TaskB_symm,1):end) = 0;
sPPI_WoD_TaskA_vs_TaskB_asymm_FDR(1:1+size(sPPI_WoD_TaskA_vs_TaskB_asymm_FDR,1):end) = 0;
sPPI_WoD_TaskA_vs_TaskB_symm_FDR(1:1+size(sPPI_WoD_TaskA_vs_TaskB_symm_FDR,1):end) = 0;

figure
try
    sgtitle('sPPI without Deconvolution No centering (After FIR task regression)')
catch
    suptitle('sPPI without Deconvolution No centering (After FIR task regression)')
end
subplot(231); imagesc(gm_sPPI_WoD_TaskIndep_asymm);       title('TaskIndep asymm');      axis square; caxis(max_ax(gm_sPPI_WoD_TaskIndep_asymm,1));          
subplot(232); imagesc(gm_sPPI_WoD_TaskA_vs_TaskB_asymm);  title('Task AvsB asymm');      axis square; caxis(max_ax(gm_sPPI_WoD_TaskA_vs_TaskB_asymm,1)); 
subplot(233); imagesc(sPPI_WoD_TaskA_vs_TaskB_asymm_FDR); title('Task AvsB asymm FDR');  axis square; 
subplot(234); imagesc(gm_sPPI_WoD_TaskIndep_symm);        title('TaskIndep symm');       axis square; caxis(max_ax(gm_sPPI_WoD_TaskIndep_symm,1));           
subplot(235); imagesc(gm_sPPI_WoD_TaskA_vs_TaskB_symm);   title('Task AvsB symm');       axis square; caxis(max_ax(gm_sPPI_WoD_TaskA_vs_TaskB_symm,1));  
subplot(236); imagesc(sPPI_WoD_TaskA_vs_TaskB_symm_FDR);  title('Task AvsB symm FDR');   axis square; 
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(233),'parula') 
colormap(subplot(236),'parula') 

sPPI_WoD_symmetry = check_symmetry(mean(sPPI_WoD_TaskA_vs_TaskB_asymm,3));
fprintf(['sPPI without deconvolution symmetry :: r = ' num2str(sPPI_WoD_symmetry) ' \n']);

%% Results: gPPI without Deconv
[gPPI_WoD_TaskA_vs_TaskB_asymm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(gPPI_WoD_TaskA_vs_TaskB_asymm,q_level);
[gPPI_WoD_TaskA_vs_TaskB_symm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(gPPI_WoD_TaskA_vs_TaskB_symm,q_level);
gm_gPPI_WoD_TaskIndep_asymm = mean(gPPI_WoD_TaskIndep_asymm,3);
gm_gPPI_WoD_TaskIndep_symm  = mean(gPPI_WoD_TaskIndep_symm,3);
gm_gPPI_WoD_TaskA_vs_TaskB_asymm = mean(gPPI_WoD_TaskA_vs_TaskB_asymm,3);
gm_gPPI_WoD_TaskA_vs_TaskB_symm = mean(gPPI_WoD_TaskA_vs_TaskB_symm,3);
gm_gPPI_WoD_TaskIndep_asymm(1:1+size(gPPI_WoD_TaskIndep_asymm,1):end) = 1;
gm_gPPI_WoD_TaskIndep_symm(1:1+size(gPPI_WoD_TaskIndep_symm,1):end) = 1;
gm_gPPI_WoD_TaskA_vs_TaskB_asymm(1:1+size(gm_gPPI_WoD_TaskA_vs_TaskB_asymm,1):end) = 0;
gm_gPPI_WoD_TaskA_vs_TaskB_symm(1:1+size(gm_gPPI_WoD_TaskA_vs_TaskB_symm,1):end) = 0;
gPPI_WoD_TaskA_vs_TaskB_asymm_FDR(1:1+size(gPPI_WoD_TaskA_vs_TaskB_asymm_FDR,1):end) = 0;
gPPI_WoD_TaskA_vs_TaskB_symm_FDR(1:1+size(gPPI_WoD_TaskA_vs_TaskB_symm_FDR,1):end) = 0;

figure
try
    sgtitle('gPPI without Deconvolution No centering (After FIR task regression)')
catch
    suptitle('gPPI without Deconvolution No centering (After FIR task regression)')
end
subplot(251); imagesc(gm_gPPI_WoD_TaskIndep_asymm);         title('TaskIndep asymm');      axis square; caxis(max_ax(gm_gPPI_WoD_TaskIndep_asymm,1));
subplot(252); imagesc(mean(gPPI_WoD_TaskA_asymm,3));        title('Task A asymm');         axis square; caxis(max_ax(mean(gPPI_WoD_TaskA_asymm,3),1));
subplot(253); imagesc(mean(gPPI_WoD_TaskB_asymm,3));        title('Task B asymm');         axis square; caxis(max_ax(mean(gPPI_WoD_TaskB_asymm,3),1));
subplot(254); imagesc(gm_gPPI_WoD_TaskA_vs_TaskB_asymm);    title('Task AvsB asymm');      axis square; caxis(max_ax(gm_gPPI_WoD_TaskA_vs_TaskB_asymm,1));
subplot(255); imagesc(gPPI_WoD_TaskA_vs_TaskB_asymm_FDR);   title('Task AvsB asymm FDR');  axis square;
subplot(256); imagesc(gm_gPPI_WoD_TaskIndep_symm);          title('TaskIndep symm');       axis square; caxis(max_ax(gm_gPPI_WoD_TaskIndep_symm,1));
subplot(257); imagesc(mean(gPPI_WoD_TaskA_symm,3));         title('Task A symm');          axis square; caxis(max_ax(mean(gPPI_WoD_TaskA_symm,3),1));
subplot(258); imagesc(mean(gPPI_WoD_TaskB_symm,3));         title('Task B symm');          axis square; caxis(max_ax(mean(gPPI_WoD_TaskB_symm,3),1));
subplot(259); imagesc(gm_gPPI_WoD_TaskA_vs_TaskB_symm);     title('Task AvsB symm');       axis square; caxis(max_ax(gm_gPPI_WoD_TaskA_vs_TaskB_symm,1));
subplot(2,5,10); imagesc(gPPI_WoD_TaskA_vs_TaskB_symm_FDR); title('Task AvsB symm FDR');   axis square;
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(255),'parula') 
colormap(subplot(2,5,10),'parula') 

gPPI_WoD_symmetry = check_symmetry(mean(gPPI_WoD_TaskA_vs_TaskB_asymm,3));
fprintf(['gPPI without deconvolution symmetry :: r = ' num2str(gPPI_WoD_symmetry) ' \n']);

%% Sensitivity and Specificity
[TPR_sPPI_WoD_asymm, TNR_sPPI_WoD_asymm] = TPR_TNR(sPPI_WoD_TaskA_vs_TaskB_asymm_FDR,ground_truth);
[TPR_sPPI_WoD_symm,  TNR_sPPI_WoD_symm]  = TPR_TNR(sPPI_WoD_TaskA_vs_TaskB_symm_FDR,ground_truth);
[TPR_gPPI_WoD_asymm, TNR_gPPI_WoD_asymm] = TPR_TNR(gPPI_WoD_TaskA_vs_TaskB_asymm_FDR,ground_truth);
[TPR_gPPI_WoD_symm,  TNR_gPPI_WoD_symm]  = TPR_TNR(gPPI_WoD_TaskA_vs_TaskB_symm_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv_No_centering_FIR_RESULTS.mat'],'TPR*','TNR*','*symmetry');
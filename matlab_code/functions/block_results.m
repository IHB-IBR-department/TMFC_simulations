function block_results(stat_path,exp_folder,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_with_Deconv.mat']);
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'sPPI_and_gPPI_without_Deconv.mat']);
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'cPPI_with_Deconv.mat']);
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'cPPI_without_Deconv.mat']);
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'CorrDiff.mat']);
load([stat_path filesep exp_folder filesep 'group_stat' filesep 'TSFC_BGFC.mat']);

%% sPPI with Deconv
[sPPI_WD_TaskA_vs_TaskB_asymm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(sPPI_WD_TaskA_vs_TaskB_asymm,q_level);
[sPPI_WD_TaskA_vs_TaskB_symm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(sPPI_WD_TaskA_vs_TaskB_symm,q_level);
gm_sPPI_WD_TaskIndep_asymm = mean(sPPI_WD_TaskIndep_asymm,3);
gm_sPPI_WD_TaskIndep_symm  = mean(sPPI_WD_TaskIndep_symm,3);
gm_sPPI_WD_TaskA_vs_TaskB_asymm = mean(sPPI_WD_TaskA_vs_TaskB_asymm,3);
gm_sPPI_WD_TaskA_vs_TaskB_symm = mean(sPPI_WD_TaskA_vs_TaskB_symm,3);
gm_sPPI_WD_TaskIndep_asymm(1:1+size(sPPI_WD_TaskIndep_asymm,1):end) = 1;
gm_sPPI_WD_TaskIndep_symm(1:1+size(sPPI_WD_TaskIndep_symm,1):end) = 1;
gm_sPPI_WD_TaskA_vs_TaskB_asymm(1:1+size(gm_sPPI_WD_TaskA_vs_TaskB_asymm,1):end) = 0;
gm_sPPI_WD_TaskA_vs_TaskB_symm(1:1+size(gm_sPPI_WD_TaskA_vs_TaskB_symm,1):end) = 0;
sPPI_WD_TaskA_vs_TaskB_asymm_FDR(1:1+size(sPPI_WD_TaskA_vs_TaskB_asymm_FDR,1):end) = 0;
sPPI_WD_TaskA_vs_TaskB_symm_FDR(1:1+size(sPPI_WD_TaskA_vs_TaskB_symm_FDR,1):end) = 0;

figure
try
    sgtitle('sPPI with Deconvolution')
catch
    suptitle('sPPI with Deconvolution')
end
subplot(231); imagesc(gm_sPPI_WD_TaskIndep_asymm);       title('TaskIndep asymm');      axis square; caxis(max_ax(gm_sPPI_WD_TaskIndep_asymm,1));          
subplot(232); imagesc(gm_sPPI_WD_TaskA_vs_TaskB_asymm);  title('Task AvsB asymm');      axis square; caxis(max_ax(gm_sPPI_WD_TaskA_vs_TaskB_asymm,1)); 
subplot(233); imagesc(sPPI_WD_TaskA_vs_TaskB_asymm_FDR); title('Task AvsB asymm FDR');  axis square; 
subplot(234); imagesc(gm_sPPI_WD_TaskIndep_symm);        title('TaskIndep symm');       axis square; caxis(max_ax(gm_sPPI_WD_TaskIndep_symm,1));           
subplot(235); imagesc(gm_sPPI_WD_TaskA_vs_TaskB_symm);   title('Task AvsB symm');       axis square; caxis(max_ax(gm_sPPI_WD_TaskA_vs_TaskB_symm,1));  
subplot(236); imagesc(sPPI_WD_TaskA_vs_TaskB_symm_FDR);  title('Task AvsB symm FDR');   axis square; 
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(233),'parula') 
colormap(subplot(236),'parula') 

fprintf(['sPPI with deconvolution assymetry :: r = ' num2str(check_symmetry(mean(sPPI_WD_TaskA_vs_TaskB_asymm,3))) ' \n']);

%% gPPI with Deconv
[gPPI_WD_TaskA_vs_TaskB_asymm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(gPPI_WD_TaskA_vs_TaskB_asymm,q_level);
[gPPI_WD_TaskA_vs_TaskB_symm_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(gPPI_WD_TaskA_vs_TaskB_symm,q_level);
gm_gPPI_WD_TaskIndep_asymm = mean(gPPI_WD_TaskIndep_asymm,3);
gm_gPPI_WD_TaskIndep_symm  = mean(gPPI_WD_TaskIndep_symm,3);
gm_gPPI_WD_TaskA_vs_TaskB_asymm = mean(gPPI_WD_TaskA_vs_TaskB_asymm,3);
gm_gPPI_WD_TaskA_vs_TaskB_symm = mean(gPPI_WD_TaskA_vs_TaskB_symm,3);
gm_gPPI_WD_TaskIndep_asymm(1:1+size(gPPI_WD_TaskIndep_asymm,1):end) = 1;
gm_gPPI_WD_TaskIndep_symm(1:1+size(gPPI_WD_TaskIndep_symm,1):end) = 1;
gm_gPPI_WD_TaskA_vs_TaskB_asymm(1:1+size(gm_gPPI_WD_TaskA_vs_TaskB_asymm,1):end) = 0;
gm_gPPI_WD_TaskA_vs_TaskB_symm(1:1+size(gm_gPPI_WD_TaskA_vs_TaskB_symm,1):end) = 0;
gPPI_WD_TaskA_vs_TaskB_asymm_FDR(1:1+size(gPPI_WD_TaskA_vs_TaskB_asymm_FDR,1):end) = 0;
gPPI_WD_TaskA_vs_TaskB_symm_FDR(1:1+size(gPPI_WD_TaskA_vs_TaskB_symm_FDR,1):end) = 0;

figure
try
    sgtitle('gPPI with Deconvolution')
catch
    suptitle('gPPI with Deconvolution')
end
subplot(251); imagesc(gm_gPPI_WD_TaskIndep_asymm);         title('TaskIndep asymm');      axis square; caxis(max_ax(gm_gPPI_WD_TaskIndep_asymm,1));
subplot(252); imagesc(mean(gPPI_WD_TaskA_asymm,3));        title('Task A asymm');         axis square; caxis(max_ax(mean(gPPI_WD_TaskA_asymm,3),1));
subplot(253); imagesc(mean(gPPI_WD_TaskB_asymm,3));        title('Task B asymm');         axis square; caxis(max_ax(mean(gPPI_WD_TaskB_asymm,3),1));
subplot(254); imagesc(gm_gPPI_WD_TaskA_vs_TaskB_asymm);    title('Task AvsB asymm');      axis square; caxis(max_ax(gm_gPPI_WD_TaskA_vs_TaskB_asymm,1));
subplot(255); imagesc(gPPI_WD_TaskA_vs_TaskB_asymm_FDR);   title('Task AvsB asymm FDR');  axis square;
subplot(256); imagesc(gm_gPPI_WD_TaskIndep_symm);          title('TaskIndep symm');       axis square; caxis(max_ax(gm_gPPI_WD_TaskIndep_symm,1));
subplot(257); imagesc(mean(gPPI_WD_TaskA_symm,3));         title('Task A symm');          axis square; caxis(max_ax(mean(gPPI_WD_TaskA_symm,3),1));
subplot(258); imagesc(mean(gPPI_WD_TaskB_symm,3));         title('Task B symm');          axis square; caxis(max_ax(mean(gPPI_WD_TaskB_symm,3),1));
subplot(259); imagesc(gm_gPPI_WD_TaskA_vs_TaskB_symm);     title('Task AvsB symm');       axis square; caxis(max_ax(gm_gPPI_WD_TaskA_vs_TaskB_symm,1));
subplot(2,5,10); imagesc(gPPI_WD_TaskA_vs_TaskB_symm_FDR); title('Task AvsB symm FDR');   axis square;
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(255),'parula') 
colormap(subplot(2,5,10),'parula') 

fprintf(['gPPI with deconvolution assymetry :: r = ' num2str(check_symmetry(gm_gPPI_WD_TaskA_vs_TaskB_asymm)) ' \n']);

%% sPPI without Deconv
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
    sgtitle('sPPI without Deconvolution')
catch
    suptitle('sPPI without Deconvolution')
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

fprintf(['sPPI without deconvolution assymetry :: r = ' num2str(check_symmetry(mean(sPPI_WoD_TaskA_vs_TaskB_asymm,3))) ' \n']);

%% gPPI without Deconv
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
    sgtitle('gPPI without Deconvolution')
catch
    suptitle('gPPI without Deconvolution')
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

fprintf(['gPPI without deconvolution assymetry :: r = ' num2str(check_symmetry(gm_gPPI_WoD_TaskA_vs_TaskB_asymm)) ' \n']);

%% cPPI with deconv and without deconv
[cPPI_WD_TaskA_vs_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(cPPI_WD_TaskA_vs_TaskB,q_level);
[cPPI_WoD_TaskA_vs_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(cPPI_WoD_TaskA_vs_TaskB,q_level);
gm_cPPI_WD_TaskA_vs_TaskB = mean(cPPI_WD_TaskA_vs_TaskB,3);
gm_cPPI_WoD_TaskA_vs_TaskB = mean(cPPI_WoD_TaskA_vs_TaskB,3);
gm_cPPI_WD_TaskA_vs_TaskB(1:1+size(gm_cPPI_WD_TaskA_vs_TaskB,1):end) = 1;
gm_cPPI_WoD_TaskA_vs_TaskB(1:1+size(gm_cPPI_WoD_TaskA_vs_TaskB,1):end) = 1;

figure
try
    sgtitle('cPPI')
catch
    suptitle('cPPI')
end
subplot(221); imagesc(gm_cPPI_WD_TaskA_vs_TaskB);   title('cPPI WD Task AvsB');       axis square; caxis(max_ax(gm_cPPI_WD_TaskA_vs_TaskB,1));
subplot(223); imagesc(gm_cPPI_WoD_TaskA_vs_TaskB);  title('cPPI WoD Task AvsB');      axis square; caxis(max_ax(gm_cPPI_WoD_TaskA_vs_TaskB,1));
subplot(222); imagesc(cPPI_WD_TaskA_vs_TaskB_FDR);  title('cPPI WD Task AvsB FDR');   axis square; 
subplot(224); imagesc(cPPI_WoD_TaskA_vs_TaskB_FDR); title('cPPI WoD Task AvsB FDR');  axis square; 
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 
colormap(subplot(222),'parula') 
colormap(subplot(224),'parula') 

%% Direct correlation difference
[corrdiff_Rest_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_Rest,q_level);
[corrdiff_TaskA_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_TaskA,q_level);
[corrdiff_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_TaskB,q_level);
[corrdiff_TaskA_vs_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(corrdiff_TaskA_vs_TaskB,q_level);
gm_corrdiff_TaskA_vs_TaskB = mean(corrdiff_TaskA_vs_TaskB,3);
gm_corrdiff_TaskA_vs_TaskB(1:1+size(gm_corrdiff_TaskA_vs_TaskB,1):end) = 0;
corrdiff_TaskA_vs_TaskB_FDR(1:1+size(corrdiff_TaskA_vs_TaskB_FDR,1):end) = 0;

figure
try
    sgtitle('Direct correlation difference')
catch
    suptitle('Direct correlation difference')
end
subplot(241); imagesc(mean(corrdiff_Rest,3)); title('Rest');  axis square; caxis(max_ax(mean(corrdiff_Rest,3),1));
subplot(242); imagesc(mean(corrdiff_TaskA,3)); title('Task A');  axis square; caxis(max_ax(mean(corrdiff_TaskA,3),1));
subplot(243); imagesc(mean(corrdiff_TaskB,3)); title('Task B');  axis square; caxis(max_ax(mean(corrdiff_TaskB,3),1));
subplot(244); imagesc(gm_corrdiff_TaskA_vs_TaskB); title('Task AvsB');  axis square; caxis(max_ax(gm_corrdiff_TaskA_vs_TaskB,1));
subplot(245); imagesc(corrdiff_Rest_FDR); title('Rest FDR');  axis square; 
subplot(246); imagesc(corrdiff_TaskA_FDR); title('Task A FDR');  axis square; 
subplot(247); imagesc(corrdiff_TaskB_FDR); title('Task B FDR');  axis square; 
subplot(248); imagesc(corrdiff_TaskA_vs_TaskB_FDR); title('Task AvsB FDR');  axis square;
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 
colormap(subplot(245),'parula') 
colormap(subplot(246),'parula')
colormap(subplot(247),'parula') 
colormap(subplot(248),'parula') 

%% TSFC and BGFC
figure
subplot(121); imagesc(mean(TSFC_group,3)); title('TSFC');  axis square; caxis(max_ax(mean(TSFC_group,3),1));
subplot(122); imagesc(mean(BGFC_group,3)); title('BGFC');  axis square; caxis(max_ax(mean(BGFC_group,3),1));
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue') 

%% Sensitivity and Specificity
T = ground_truth;
TMFC_FDR_results(:,:,1) = sPPI_WD_TaskA_vs_TaskB_symm_FDR;
TMFC_FDR_results(:,:,2) = sPPI_WD_TaskA_vs_TaskB_asymm_FDR;
TMFC_FDR_results(:,:,3) = gPPI_WD_TaskA_vs_TaskB_symm_FDR;
TMFC_FDR_results(:,:,4) = gPPI_WD_TaskA_vs_TaskB_asymm_FDR;
TMFC_FDR_results(:,:,5) = sPPI_WoD_TaskA_vs_TaskB_symm_FDR;
TMFC_FDR_results(:,:,6) = sPPI_WoD_TaskA_vs_TaskB_asymm_FDR;
TMFC_FDR_results(:,:,7) = gPPI_WoD_TaskA_vs_TaskB_symm_FDR;
TMFC_FDR_results(:,:,8) = gPPI_WoD_TaskA_vs_TaskB_asymm_FDR;
TMFC_FDR_results(:,:,9) = cPPI_WD_TaskA_vs_TaskB_FDR;
TMFC_FDR_results(:,:,10) = cPPI_WoD_TaskA_vs_TaskB_FDR;
TMFC_FDR_results(:,:,11) = corrdiff_TaskA_vs_TaskB_FDR;

for i=1:11
    P = TMFC_FDR_results(:,:,i);
    P(1:1+size(P,1):end) = 0;
    FP = sum(sum((1-T).*P));
    FN = sum(sum(T.*(1-P)));
    TP = sum(sum(T.*P));
    TN = sum(sum((1-T).*(1-P)));
    TPR(i) = TP./(TP+FN);
    TNR(i) = TN./(TN+FP);
    PPV(i) = TP./(TP+FP);
    NPV(i) = TN./(TN+FN);
    FDR(i) = 1-PPV(i);
    clear P FP FN TP TN
end

save([stat_path filesep exp_folder filesep 'group_stat' filesep 'all_results_[q_level_' num2str(q_level) '].mat'],'TPR','TNR','PPV','NPV','FDR')



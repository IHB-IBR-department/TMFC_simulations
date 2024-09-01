function cPPI_with_deconv(stat_path,exp_folder,N,N_ROIs,q_level,ground_truth)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'PPI_ALL.mat']);
    partcorr_TaskA_vs_TaskB_cPPI = NaN(N_ROIs,N_ROIs);  

    for j = 1:N_ROIs
        for k = j+1:N_ROIs
            % Time courses to be correlated
            X_PPI_TaskA_vs_TaskB = [];
            X_PPI_TaskA_vs_TaskB = [PPI_TaskA_vs_TaskB(:,j),PPI_TaskA_vs_TaskB(:,k)];

            % Confound matrix
            Y_cPPI = [];
            Y_cPPI = [VOI(:,j), VOI(:,k), Psy_TaskA_vs_TaskB];

            % Compute partial correlations
            part_corr = [];
            part_corr = atanh(partialcorr(X_PPI_TaskA_vs_TaskB,Y_cPPI));  
            partcorr_TaskA_vs_TaskB_cPPI(j,k) = part_corr(2,1);
            partcorr_TaskA_vs_TaskB_cPPI(k,j) = part_corr(2,1);

        end
    end

    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'cPPI_with_Deconv.mat'],'partcorr_TaskA_vs_TaskB_cPPI');
    time(subji)=toc;
    fprintf(['cPPI with deconvolution :: Subject #' num2str(subji) ' :: Done in: ' num2str(time(subji)) 's \n']);    

    clearvars -except stat_path exp_folder N N_ROIs q_level ground_truth
end

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'cPPI_with_Deconv.mat']);
    cPPI_WD_TaskA_vs_TaskB(:,:,subji) = partcorr_TaskA_vs_TaskB_cPPI; 
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'cPPI_with_Deconv.mat'],'cPPI*');
time = toc;
fprintf(['Merge cPPI_with_Deconv matrices :: Done in: ' num2str(time) 's \n']); 

%% Results: cPPI with Deconv
[cPPI_WD_TaskA_vs_TaskB_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(cPPI_WD_TaskA_vs_TaskB,q_level);
gm_cPPI_WD_TaskA_vs_TaskB = mean(cPPI_WD_TaskA_vs_TaskB,3);
gm_cPPI_WD_TaskA_vs_TaskB(1:1+size(gm_cPPI_WD_TaskA_vs_TaskB,1):end) = 1;

figure
try
    sgtitle('cPPI with Deconvolution')
catch
    suptitle('cPPI with Deconvolution')
end
subplot(121); imagesc(gm_cPPI_WD_TaskA_vs_TaskB);  title('cPPI Task AvsB');      axis square; caxis(max_ax(gm_cPPI_WD_TaskA_vs_TaskB,1));
subplot(122); imagesc(cPPI_WD_TaskA_vs_TaskB_FDR); title('cPPI Task AvsB FDR');  axis square; 

set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')
colormap(subplot(122),'parula') 

%% Sensitivity and Specificity
[TPR_cPPI_WD, TNR_cPPI_WD] = TPR_TNR(cPPI_WD_TaskA_vs_TaskB_FDR,ground_truth);

%% Save results
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'cPPI_with_Deconv_RESULTS.mat'],'TPR*','TNR*');


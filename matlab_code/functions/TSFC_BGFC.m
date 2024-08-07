function TSFC_BGFC(stat_path,exp_folder,N)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);
    TSFC = atanh(corr(VOI));  
    clear VOI

    load([stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'VOI_ALL.mat']);
    BGFC = atanh(corr(VOI));  
    clear VOI

    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'TSFC_BGFC.mat'],'TSFC','BGFC');
    clear TSFC BGFC
end

time=toc;
fprintf(['TSFC and BGFC :: Done in: ' num2str(time) 's \n']); 

%% Merge matrices
tic
mkdir([stat_path filesep exp_folder filesep 'group_stat']);
for subji = 1:N
    load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'TSFC_BGFC.mat']);
    TSFC_group(:,:,subji) = TSFC;
    BGFC_group(:,:,subji) = BGFC; 
end
save([stat_path filesep exp_folder filesep 'group_stat' filesep 'TSFC_BGFC.mat'],'TSFC*','BGFC*');
time = toc;
fprintf(['Merge TSFC and BGFC matrices :: Done in: ' num2str(time) 's \n']); 

%% Results
figure
subplot(121); imagesc(mean(TSFC_group,3)); title('TSFC');  axis square; caxis(max_ax(mean(TSFC_group,3),1));
subplot(122); imagesc(mean(BGFC_group,3)); title('BGFC');  axis square; caxis(max_ax(mean(BGFC_group,3),1));

set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap('redblue')


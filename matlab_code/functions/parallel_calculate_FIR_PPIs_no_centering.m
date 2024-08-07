function parallel_calculate_FIR_PPIs_no_centering(stat_path,exp_folder,N,N_ROIs)

% ========================================================================
% Ruslan Masharipov, July, 2024
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

tic
for subji = 1:N
    for j = 1:N_ROIs
        %% PPI_Task_A for gPPI 
        matlabbatch{1}.spm.stats.ppi.spmmat = ...
            {[stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.ppi.type.ppi.voi = ...
            {[stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'VOI_ROI_' num2str(j,'%03.f') '_1.mat']};
        matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 1 1];
        matlabbatch{1}.spm.stats.ppi.name = ['No_centering_FIR_TaskA_' num2str(j,'%.3d')];
        matlabbatch{1}.spm.stats.ppi.disp = 0;
        batch_A{j} = matlabbatch;
        clear matlabbatch     
       
        %% PPI_Task_B for gPPI 
        matlabbatch{1}.spm.stats.ppi.spmmat = ...
            {[stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.ppi.type.ppi.voi = ...
            {[stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'VOI_ROI_' num2str(j,'%03.f') '_1.mat']};
        matlabbatch{1}.spm.stats.ppi.type.ppi.u = [2 1 1];
        matlabbatch{1}.spm.stats.ppi.name = ['No_centering_FIR_TaskB_' num2str(j,'%.3d')];
        matlabbatch{1}.spm.stats.ppi.disp = 0;
        batch_B{j} = matlabbatch;
        clear matlabbatch
        
        %% PPI_TaskA_vs_TaskB for sPPI 
        matlabbatch{1}.spm.stats.ppi.spmmat = ...
            {[stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.ppi.type.ppi.voi = ...
            {[stat_path filesep exp_folder filesep 'FIR_GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'VOI_ROI_' num2str(j,'%03.f') '_1.mat']};
        matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 1 1; 2 1 -1];
        matlabbatch{1}.spm.stats.ppi.name = ['No_centering_FIR_TaskA_vs_TaskB_' num2str(j,'%.3d')];
        matlabbatch{1}.spm.stats.ppi.disp = 0;
        batch_AB{j} = matlabbatch;
        clear matlabbatch

    end

    parfor j = 1:N_ROIs
        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_get_defaults('cmdline',true);
        spm_get_defaults('stats.resmem',1);
        spm_get_defaults('stats.maxmem',2^33);
        spm_jobman('run',batch_A{j});
        spm_jobman('run',batch_B{j});
        spm_jobman('run',batch_AB{j});
    end

    clear batch_A batch_B batch_AB
    time(subji)=toc;
    fprintf(['Calculate PPIs :: Subject #' num2str(subji) ' :: Done in ' num2str(time(subji)) 's \n']);
end

%% Merge PPIs
tic
for subji = 1:N
    for j = 1:N_ROIs
        TaskA = load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'PPI_No_centering_FIR_TaskA_' num2str(j,'%03.f') '.mat']);
        TaskB = load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'PPI_No_centering_FIR_TaskB_' num2str(j,'%03.f') '.mat']);
        TaskA_vs_TaskB = load([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep 'PPI_No_centering_FIR_TaskA_vs_TaskB_' num2str(j,'%03.f') '.mat']);
        
        VOI(:,j) = TaskA.PPI.Y;
                
        Psy_TaskA = TaskA.PPI.P;
        Psy_TaskB = TaskB.PPI.P;
        Psy_TaskA_vs_TaskB = TaskA_vs_TaskB.PPI.P;
                
        PPI_TaskA(:,j) = TaskA.PPI.ppi;
        PPI_TaskB(:,j) = TaskB.PPI.ppi;
        PPI_TaskA_vs_TaskB(:,j) = TaskA_vs_TaskB.PPI.ppi;
    end    
          
    save([stat_path filesep exp_folder filesep 'GLMs' filesep 'Sub_' num2str(subji,'%.3d') filesep  'No_centering_PPI_FIR_ALL.mat'],...
        'VOI','Psy_TaskA','Psy_TaskB','Psy_TaskA_vs_TaskB',...
              'PPI_TaskA','PPI_TaskB','PPI_TaskA_vs_TaskB');
   
    time(subji)=toc;
    fprintf(['Merge PPIs :: Subject #' num2str(subji) ' :: Done in: ' num2str(time(subji)) 's \n']);       
    
    clearvars -except stat_path exp_folder N N_ROIs
end

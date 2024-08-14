# Simulated BOLD time series 

Task design *.mat file is used as an input for TMFC simulations using Python code. 

The output of the Python code is simulated BOLD time series *.mat file. 

Task design and simulated BOLD time series *.mat files are used as input for TMFC analysis using MATLAB code. 

## "Simulated BOLD time series" data structures

* oscill - oscillatory-related BOLD signal, [m x n x o] matrix, where m - time point, n - brain region, o - subject   
* coact - coactivation-related BOLD signal, [m x n] matrix, where m - time point, n - brain region

## Experiments:

* [SIM_BOLD 01_BLOCK [2s_TR] [20s_DUR] [10_BLOCKS]](/simulated_BOLD_time_series/files/SIM_BOLD_01_BLOCK_[2s_TR]_[20s_DUR]_[10_BLOCKS].mat)
   * Block design
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 10 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 13.3 min

* [SIM_BOLD 01_BLOCK [2s_TR] [20s_DUR] [10_BLOCKS] [**NO_TMFC**]](/simulated_BOLD_time_series/files/SIM_BOLD_01_BLOCK_[2s_TR]_[20s_DUR]_[10_BLOCKS]_[NO_TMFC].mat)
   * Block design
   * **No task-modulated functional connectivity** (synaptic weight matrices are the same for Rest, Cond A and Cond B)
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 10 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 13.3 min

* [SIM_BOLD 01_BLOCK [2s_TR] [20s_DUR] [10_BLOCKS] [**Variable_HRF**]](/simulated_BOLD_time_series/files/SIM_BOLD_01_BLOCK_[2s_TR]_[20s_DUR]_[10_BLOCKS]_[Variable_HRF].mat)
   * Block design
   * **Variable parameters of the Balloon-Windkessel hemodynamic model (time-to-peak = 3-7 s)**
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 10 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 13.3 min

* [SIM_BOLD 02_EVENT [2s_TR] [1s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_02_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Default event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 23.6 min

* [SIM_BOLD 02_EVENT [2s_TR] [1s_DUR] [6s_ISI] [100_TRIALS] [**Variable_HRF**]](/simulated_BOLD_time_series/files/SIM_BOLD_02_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS]_[Variable_HRF].rar)
   * Default event-related design
   * **Variable parameters of the Balloon-Windkessel hemodynamic model (time-to-peak = 3-7 s)**
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 23.6 min

* [SIM_BOLD 03_EVENT [500ms_TR] [1s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_03_EVENT_[500ms_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].zip)
   * Event-related design
   * TR = 500 ms
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 14 time points (7 s)
   * Total scan time = 23.6 min
     
* [SIM_BOLD 04_EVENT [700ms_TR] [1s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_04_EVENT_[700ms_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].zip)
   * Event-related design
   * TR = 700 ms
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 10 time points (7 s)
   * Total scan time = 23.6 min

* [SIM_BOLD 05_EVENT [1s_TR] [1s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_05_EVENT_[1s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 1 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 6 time points (6 s)
   * Total scan time = 23.6 min

* [SIM_BOLD 06_EVENT [3s_TR] [1s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_06_EVENT_[3s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 3 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (9 s)
   * Total scan time = 23.6 min
 
* [SIM_BOLD 07_EVENT [2s_TR] [100ms_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_07_EVENT_[2s_TR]_[100ms_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 100 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 20.6 min
 
* [SIM_BOLD 08_EVENT [2s_TR] [250ms_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_08_EVENT_[2s_TR]_[250ms_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 250 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 21.1 min
 
* [SIM_BOLD 09_EVENT [2s_TR] [500ms_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_09_EVENT_[2s_TR]_[500ms_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 500 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 22 min
 
* [SIM_BOLD 10_EVENT [2s_TR] [2s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_10_EVENT_[2s_TR]_[2s_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 2 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 27 min
 
* [SIM_BOLD 11_EVENT [2s_TR] [4s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_11_EVENT_[2s_TR]_[4s_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 4 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 33.8 min
 
* [SIM_BOLD 12_EVENT [2s_TR] [1s_DUR] [2s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_12_EVENT_[2s_TR]_[1s_DUR]_[2s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 1-3 s (mean ISI = 2 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 10.4 min
 
* [SIM_BOLD 13_EVENT [2s_TR] [1s_DUR] [4s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_13_EVENT_[2s_TR]_[1s_DUR]_[4s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 2-6 s (mean ISI = 4 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 17 min
 
* [SIM_BOLD 14_EVENT [2s_TR] [1s_DUR] [8s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_14_EVENT_[2s_TR]_[1s_DUR]_[8s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 6-10 s (mean ISI = 8 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 30.4 min
 
* [SIM_BOLD 15_EVENT [2s_TR] [1s_DUR] [12s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_15_EVENT_[2s_TR]_[1s_DUR]_[12s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 10-14 s (mean ISI = 12 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 44 min

* [SIM_BOLD 16_EVENT [2s_TR] [1s_DUR] [6s_ISI] [20_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_16_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[20_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 20 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 5 min
 
* [SIM_BOLD 17_EVENT [2s_TR] [1s_DUR] [6s_ISI] [40_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_17_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[40_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 40 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 9.7 min

* [SIM_BOLD 18_EVENT [2s_TR] [1s_DUR] [6s_ISI] [60_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_18_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[60_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 60 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 14 min
 
* [SIM_BOLD 19_EVENT [2s_TR] [1s_DUR] [6s_ISI] [80_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_19_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[80_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 80 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 19 min
 
* [SIM_BOLD 20_EVENT [2s_TR] [100ms_DUR] [6s_ISI] [200_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_20_EVENT_[2s_TR]_[100ms_DUR]_[6s_ISI]_[200_TRIALS].mat)
   * Event-related design
   * TR = 2 s
   * Event duration = 100 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 200 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 40.9 min
 
* [SIM_BOLD 21_EVENT [500ms_TR] [100ms_DUR] [6s_ISI] [50_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_21_EVENT_[500ms_TR]_[100ms_DUR]_[6s_ISI]_[50_TRIALS].mat)
   * Event-related design
   * TR = 500 ms
   * Event duration = 100 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 50 events per condition
   * Dummy scans: first 14 time points (7 s)
   * Total scan time = 10.3 min
 
* [SIM_BOLD 22_EVENT [700ms_TR] [100ms_DUR] [6s_ISI] [70_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_22_EVENT_[700ms_TR]_[100ms_DUR]_[6s_ISI]_[70_TRIALS].mat)
   * Event-related design
   * TR = 700 ms
   * Event duration = 100 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 70 events per condition
   * Dummy scans: first 10 time points (7 s)
   * Total scan time = 14.4 min
 
* [SIM_BOLD 23_EVENT [1s_TR] [100ms_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series/files/SIM_BOLD_23_EVENT_[1s_TR]_[100ms_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Event-related design
   * TR = 1 s
   * Event duration = 100 ms
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 6 time points (6 s)
   * Total scan time = 20.6 min
 
* [SIM_BOLD 24_BLOCK [2s_TR] [20s_DUR] [10_BLOCKS] [**ASYMMETRIC**]](/simulated_BOLD_time_series/files/SIM_BOLD_24_BLOCK_[2s_TR]_[20s_DUR]_[10_BLOCKS]_[ASYMMETRIC].mat)
   * Block design
   * **Asymmetric synaptic weight matrices**
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 10 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 13.3 min
 
* [SIM_BOLD 25_BLOCK [2s_TR] [20s_DUR] [20_BLOCKS] [**ASYMMETRIC**]](/simulated_BOLD_time_series/files/SIM_BOLD_25_BLOCK_[2s_TR]_[20s_DUR]_[20_BLOCKS]_[ASYMMETRIC].mat)
   * Block design
   * **Asymmetric synaptic weight matrices**
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 20 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 26.6 min
 
* [SIM_BOLD 25_BLOCK [2s_TR] [20s_DUR] [20_BLOCKS] [**ASYMMETRIC**] [**Variable_HRF**]](/simulated_BOLD_time_series/files/SIM_BOLD_25_BLOCK_[2s_TR]_[20s_DUR]_[20_BLOCKS]_[ASYMMETRIC]_[Variable_HRF].rar)
   * Block design
   * **Asymmetric synaptic weight matrices**
   * **Variable parameters of the Balloon-Windkessel hemodynamic model (time-to-peak = 3-7 s)**
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 20 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 26.6 min
     
* [SIM_BOLD 26_EVENT [2s_TR] [1s_DUR] [6s_ISI] [100_TRIALS] [**ASYMMETRIC**]](/simulated_BOLD_time_series/files/SIM_BOLD_26_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS]_[ASYMMETRIC].mat)
   * Event-related design
   * **Asymmetric synaptic weight matrices**
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 23.6 min
 
* [SIM_BOLD 27_EVENT [2s_TR] [1s_DUR] [6s_ISI] [200_TRIALS] [**ASYMMETRIC**]](/simulated_BOLD_time_series/files/SIM_BOLD_27_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[200_TRIALS]_[ASYMMETRIC].mat)
   * Event-related design
   * **Asymmetric synaptic weight matrices**
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 200 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 40.9 min
 
* [SIM_BOLD 27_EVENT [2s_TR] [1s_DUR] [6s_ISI] [200_TRIALS] [**ASYMMETRIC**] [**Variable_HRF**]](/simulated_BOLD_time_series/files/SIM_BOLD_27_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[200_TRIALS]_[ASYMMETRIC]_[Variable_HRF].part1.rar)
   * Event-related design
   * **Asymmetric synaptic weight matrices**
   * **Variable parameters of the Balloon-Windkessel hemodynamic model (time-to-peak = 3-7 s)**
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 200 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 40.9 min

 * [SIM_BOLD 28_BLOCK [720ms_TR] [HCP_WM]](/simulated_BOLD_time_series/files/SIM_BOLD_28_BLOCK_[720ms_TR]_[HCP_WM].mat)
   * Block design (Human Connectome Project, Working Memory Task)
   * Repetition time (TR) = 720 ms
   * Task block duration = 27.5 s (2-back, 1-back)
   * Rest (R) = 15 s
   * 8 blocks per condition
   * Block sequence: [2-b, 1-b, R, 2-b, 1-b, R, 2-b, 2-b, R, 1-b, 1-b, R, 2-b, 1-b, R, 2-b, 1-b, R, 1-b, 2-b, R, 1-b, 2-b]
   * Dummy scans: first 25 time points (18 s)
   * Total scan time = 10 min

* [SIM_BOLD 29_BLOCK [720ms_TR] [HCP_SOCIAL]](/simulated_BOLD_time_series/files/SIM_BOLD_29_BLOCK_[720ms_TR]_[HCP_SOCIAL].mat)
   * Block design (Human Connectome Project, Social Cognition Task)
   * Repetition time (TR) = 720 ms
   * Task block duration = 23 s (Social, Random)
   * Rest (R) = 15 s
   * 5 blocks per condition
   * Block sequence: [Soc, R, Soc, R, Rnd, R, Soc, R, Rnd, Soc, R, Rnd, R, Rnd, R, Soc, R, Rnd, R]
   * Dummy scans: first 25 time points (18 s)
   * Total scan time = 7 min
 
* [SIM_BOLD 30_EVENT [2s_TR] [CNP_SST] [All_stop_trials]](/simulated_BOLD_time_series/files/SIM_BOLD_30_EVENT_[2s_TR]_[CNP_SST]_[All_stop_trials].mat)
   * Event-related design (UCLA Consortium for Neuropsychiatric Phenomics, Stop Signal Task)
   * TR = 2 s
   * Event duration = 1.5 s
   * Random interstimulus interval (ISI) = 0.5-4 s (mean ISI = 1 s)
   * 96 GO trials
   * 32 STOP-signal trials
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 7 min

* [SIM_BOLD 31_EVENT [2s_TR] [CNP_SST] [Correct_stop_trials]](/simulated_BOLD_time_series/files/SIM_BOLD_31_EVENT_[2s_TR]_[CNP_SST]_[Correct_stop_trials].mat)
   * Unbalanced event-related design (UCLA Consortium for Neuropsychiatric Phenomics, Stop Signal Task)
   * TR = 2 s
   * Event duration = 1.5 s
   * Random interstimulus interval (ISI) = 0.5-4 s (mean ISI = 1 s)
   * 94 GO trials
   * 17 Correct STOP-signal trials
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 7 min

 * [SIM_BOLD 32_EVENT [2s_TR] [CNP_TST] [Switch_vs_NoSwitch]](/simulated_BOLD_time_series/files/SIM_BOLD_32_EVENT_[2s_TR]_[CNP_TST]_[Switch_vs_NoSwitch].mat)
   * Unbalanced event-related design (UCLA Consortium for Neuropsychiatric Phenomics, Task Switching Task)
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI), mean ISI = 3 s
   * 24 Switch trials
   * 72 No Switch trials
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 6 min
  
* [SIM_BOLD 33_BLOCK [2s_TR] [20s_DUR] [20_BLOCKS] [**Variable_HRF**]](/simulated_BOLD_time_series/files/SIM_BOLD_33_BLOCK_[2s_TR]_[20s_DUR]_[20_BLOCKS]_[Variable_HRF].rar)
   * Block design
   * **Variable parameters of the Balloon-Windkessel hemodynamic model (time-to-peak = 3-7 s)**
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 20 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 26.6 min

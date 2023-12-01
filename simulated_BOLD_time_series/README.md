# Simulated BOLD time series 

## Experiments:

1. [SIM_BOLD 01_BLOCK [2s_TR] [20s_DUR] [10_BLOCKS]](/simulated_BOLD_time_series/SIM_BOLD_01_BLOCK_[2s_TR]_[20s_DUR]_[10_BLOCKS].mat)
   * Block design
   * Repetition time (TR) = 2 s
   * Block duration = 20 s
   * 10 blocks per condition
   * Block sequence: [Cond_A, Rest, Cond_B, Rest, Cond_B, ... ]
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 13.3 min

2. [SIM_BOLD 02_EVENT [2s_TR] [1s_DUR] [6s_ISI] [100_TRIALS]](/simulated_BOLD_time_series//SIM_BOLD_02_EVENT_[2s_TR]_[1s_DUR]_[6s_ISI]_[100_TRIALS].mat)
   * Default event-related design
   * TR = 2 s
   * Event duration = 1 s
   * Random interstimulus interval (ISI) = 4-8 s (mean ISI = 6 s)
   * 100 events per condition
   * Dummy scans: first 3 time points (6 s)
   * Total scan time = 23.6 min

# Comparison of Task-Modulated Functional Connectivity (TMFC) Methods

Repository for code and simulations from  <br/> 
**"Comparison of whole-brain task-modulated functional connectivity methods for fMRI task connectomics"** <br/>
by Masharipov, R., Knyazeva, I., Korotkov, A., Cherednichenko, D. &amp; Kireev, M.

Use the repository [Discussions](https://github.com/Masharipov/TMFC_simulations/discussions) for questions or email masharipov@ihb.spb.ru

## Overview

Here, we provide:

1. [Task design files](task_designs) (.*mat format) containing stimulus onsets, durations, condition names, and weighting factors for synaptic matrices. <br/>
   **Onsets**, **durations** and **condition names** are defined in the same way as for **multiple conditions** *.mat file for SPM12.
  
3. [Python code](python_code) for TMFC simulations based on **large-scale Wilson-Cowan neural mass model** and **Ballon-Windkessel haemodynamic model**.

4. User-friendly [Jupyter notebooks](jupyter_notebooks) for reproducing our simulations.

5. [Simulated BOLD time series files](simulated_BOLD_time_series) (*.mat format) for all experiments presented in the paper.

6. SPM12-based [MATLAB code](matlab_code) for TMFC analysis using:
    * correlation difference approach (**CorrDiff**),
    * standard psychophysiological interactions (**sPPI**),
    * generalised psychophysiological interactions (**gPPI**),
    * correlational psychophysiological interactions (**cPPI**),
    * beta-series correlations based on least-squares-all approach (**BSC-LSA**),
    * beta-series correlations based on least-squares-separate approach (**BSC-LSS**).

## TMFC Simulations

## TMFC Analysis




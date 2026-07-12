# Reproduce — Trial-design simulation (Parikh et al., Ann Clin Transl Neurol 2024)

MATLAB (R2018+, verified on R2026a); Python (for `Curve_fit_Python.py`, called by Step 1).

```matlab
Step_2_Calculate_Sample_Effect_Size   % sample sizes + effect sizes from the 48-patient parameters
```
`Step_2` reads the committed `Parameters 48 patients final.xlsx` and, via `Drug_Outcomes.m`,
simulates the trial designs to produce the required sample size / effect size for each outcome
definition (complete seizure cessation, etc.) — the paper's headline power/sample-size results.

To rebuild the parameters from the per-patient time series (fuller pipeline):
```matlab
Step_1_Calculate_parameters           % fits log-normal + drug weights from Spike/CNN/Drug data -> parameters xlsx
```
Figure code is under `Code_To_Generate_Figures_Paper/`; `Sim_Results_Store/` holds the paper's
committed simulation outputs. See `DATA_SOURCE.md`.

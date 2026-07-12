# Data source & provenance — Trial-design simulation for anti-seizure treatment (Parikh et al. 2024)

## Paper
Parikh H, Sun H, Amerineni R, Rosenthal ES, Volfovsky A, Rudin C, **Westover MB**, Zafar SF.
*How many patients do you need? Investigating trial designs for anti-seizure treatment in
acute brain injury patients.* **Ann Clin Transl Neurol 2024;11(7):1681-1690.**
doi:10.1002/acn3.52059 · PMID 38867375.

## Committed data (de-identified)
Analysis inputs are keyed by de-identified subject id (`sidNNNN`) — **no MRNs, names, or dates**:
- `Parameters 48 patients final.xlsx` — per-patient fitted parameters (log-normal seizure/spike
  dynamics `mu`,`Sigma`,`Peak_Value` + drug-response weights `DC_W1..W4`) for the 48 analyzed patients.
- `Spike Data/sidNNNN_Spike_Artifacts.mat` — per-patient 2-second-window spike-probability +
  artifact signals (`Buff_yp`, `Buff_art`).
- `CNN_Label/`, `Combined_Drug_Normalize/` — per-patient CNN cEEG labels and normalized drug-
  concentration time series (2-s windows; `sidNNNN`).
- `Sim_Results_Store/` — the simulation results reported in the paper.

Cohort: 48-50 patients from the MGH SAH SAGE continuous-EEG database.

## PHI — NOT in this repo
The raw medication-administration-record (MAR) EPIC exports (`sidNNNN_MAR_Pre_Post_Epic.xlsx`,
~100 files) and the `MRN/Name → sid` linkage (`MNR_SID_50_*.xlsx`, `SAGE_MRN.xlsx`,
`Info 50 clients…xlsx`) contain **MRNs, names, and real dates** and are **excluded / git-ignored**.
The committed `Combined_Drug_Normalize` series are the de-identified, normalized outputs derived
from those raw MARs. PHI-scanned: committed data = 0 MRNs/names/dates.

## Lineage
SAH SAGE cEEG + MAR (raw, PHI) -> per-patient spike/label/drug time series (`sidNNNN`) ->
`Step_1_Calculate_parameters.m` (+ `Curve_fit_Python.py`) fits log-normal seizure + drug-response
parameters -> `Parameters 48 patients final.xlsx` -> `Step_2_Calculate_Sample_Effect_Size.m`
(+ `Drug_Outcomes.m`) simulates trial designs and computes required **sample sizes / effect sizes**.

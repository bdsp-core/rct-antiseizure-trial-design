# How many patients do you need? — trial-design simulation for anti-seizure treatment

> **Code repository.** This repo holds the analysis code; the de-identified **data** live in the bdsp.io **credentialed release** (data → AWS, code → GitHub):
> [bdsp.io/content/nu7k0916y63wxk4j2dbi/1.0.0](https://bdsp.io/content/nu7k0916y63wxk4j2dbi/1.0.0/) · DOI [10.60508/chjy-bf16](https://doi.org/10.60508/chjy-bf16) · `s3://bdsp-opendata-credentialed/rct-antiseizure-trial-design/`


MATLAB/Python code and de-identified data for:

> Parikh H, Sun H, Amerineni R, Rosenthal ES, Volfovsky A, Rudin C, **Westover MB**, Zafar SF.
> *How many patients do you need? Investigating trial designs for anti-seizure treatment in
> acute brain injury patients.* **Ann Clin Transl Neurol 2024;11(7):1681-1690.**
> [doi:10.1002/acn3.52059](https://doi.org/10.1002/acn3.52059) · PMID 38867375

A simulation framework that uses fitted patient-level seizure dynamics and drug-response models
(from a 48-patient SAH continuous-EEG cohort) to estimate the **sample sizes and effect sizes**
needed for randomized trials of anti-seizure treatment in critically ill patients.

## Reproduce
MATLAB — see **[REPRODUCE.md](REPRODUCE.md)**:
```matlab
Step_2_Calculate_Sample_Effect_Size
```

## Data
De-identified per-patient inputs keyed by `sidNNNN` (`Parameters 48 patients final.xlsx`, Spike/CNN/Drug
time series); **no MRNs, names, or dates** (raw medication + MRN-linkage files are excluded). See
**[DATA_SOURCE.md](DATA_SOURCE.md)**.

## License
BDSP credentialed data terms.

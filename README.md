# nes-lter-rate-calculation-chl
You can find here all the clean and quality checkd Chl-a data from all nes-lter cruises and the matlab scripts to perform the automated calculation of the apparent growth rates k for each treatment and filter types, and the rates calculation. Each cruise has separated csv files and data for all cruises are processed simultaneously.
The automated Quality Check (QC) of the data is detailled in https://github.com/pmarrec/nes-lter-chl-cleaning

The **chl_grazing_experiment_k_values.m** Matlab script is used to compute the apparent growth rates k obtained during the dilution (grazing) experiments.
The duration of each incuabtion is calculated from the date_time_utc_end and the date_time_utc_start for each experiment.
All the triplicate T0 Chl-a conc are calculated (only the one with iode_quality_flag = 1). k are caluclated from each TF values. Up to 6 k values are then obatined for each treatment (dilution/nutrient/light) and filter type (>0&<200, >10&<200, >0 and >0&<10 ).
Based on >0&200 (GFF) and >10&<200 (10um filters, u10 = up 10), >0&<10 (d10 = down 10) Chl-a values are calculated and then corresponding k values. For each triplicate, Chl-a d10 triplcate values are calculated as the diffrenece between the mean Chl-a value on >0&<200 and individual triplicate Chl-a values of >10&<200. Only dat with QC = 1 are considered. If Chl-a d10 <0, then Chl-a conc and k = NaN.

The **chl-grazing_experiment_rates.m** Matlab script is used to compute phytoplnkton growth rates (mu0), microzooplankton grazing rates (g), apparent growth rates in nonamended nutrient treatment (wsw NoN, kNoN) and phytoplankton growth rates in nutrient amended treatments (muN = g + kN). Associted errors (std) were estimated for all these rates. Note that kNoN and muN were estimated only when there was apparent nutrient limitation.
More details to come.

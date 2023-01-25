# nes-lter-rate-calculation-chl
You can find here all the clean and quality checkd Chl-a data from all nes-lter cruises and the matlab scripts to perform the automated calculation of the apparent growth rates k for each treatment and filter types, and the rates calculation. Each cruise has separated csv files and data for all cruises are processed simultaneously.\
The automated Quality Check (QC) of the data is detailled in https://github.com/pmarrec/nes-lter-chl-cleaning

The **chl_grazing_experiment_k_values.m** Matlab script is used to compute the apparent growth rates k obtained during the dilution (grazing) experiments.\
Create a new table for each cruise gathering all the caluclated k-values (apparent groth rates) for each tretment and flter type.\
The duration of each incuabtion is calculated from the date_time_utc_end and the date_time_utc_start for each experiment.\
All the triplicate T0 Chl-a conc are calculated (only the one with iode_quality_flag = 1).\
T0 WSW mean Chla (Total = Chla, >10um = Chlau10, <10um = Chlad10) are reported in the new table, in addition to the % of Chl-a </>10um (Chlad10per and Chlau10per). Chlad10 = Chla - Chlau10. If Chlad10 < 0, Chlad10 = 0, Chlad10per = 0% and Chlau10per = 100%. \
k are caluclated from each TF values. Up to 6 k (9 k values from ar66b, new experimental design) values are then obatined for each treatment (dilution/nutrient/light) and filter type (>0&<200, >10&<200, >0 and >0&<10 ).\
Based on >0&200 (GFF) and >10&<200 (10um filters), >0&<10 (from difference between Chl-a conc. in >0&<200 and >10&<200 size fractions) Chl-a values are calculated and then corresponding k values. For each triplicate, Chl-a d10 triplcate values are calculated as the diffrenece between the mean Chl-a value on >0&<200 and individual triplicate Chl-a values of >10&<200.\ Only dat with QC = 1 are considered. If Chl-a d10 <0, then Chl-a conc and k = NaN.\
*Input: CRUSIE-chl-grazing-experiments-clean.csv files*\
*Outputs: CRUISE-chla-grazing-experiments-k-values.csv files.*\

The **chl-grazing_experiment_rates.m** Matlab script is used to compute phytoplnkton growth rates (mu0), microzooplankton grazing rates (g), apparent growth rates in nonamended nutrient treatment (wsw NoN, kNoN) and phytoplankton growth rates in nutrient amended treatments (muN = g + kN). Associted errors (std) were estimated for all these rates. Note that kNoN and muN were estimated only when there was apparent nutrient limitation.\
Create a new table for each cruise gathering all the caluclated rates for each tretment and flter type from the k (apparent growth rates) values and the associated data and metadata.\
Phytoplankton growth rates and protist grazing rates were estimated from 24 h changes in Chl-a concentration. For each incubation bottle, the apparent growth rates (k, d^-1) were calculated as:\ 
k=1⁄t×ln(C_t⁄C_0) \
where t is the incubation time (d) and C_t and C_0 the final and initial Chl-a concentration (µg L^-1), respectively.\
Protist grazing rates (g, d^-1) were estimated with the equation:\
g=((k_d-k_N))⁄((1-x))\
where k_d and k_N are the apparent growth rates in 20WSW and WSW nutrient amended treatments, respectively, and x is the achieved fraction of WSW in the diluted treatment calculated from T0 Chl-a in 20WSW and WSW. \
Accordingly, the instantaneous, or in situ, growth rate (mu_0, d^-1) was estimated as:\
mu_0=g+k_NoN\
where k_NoN is apparent phytoplankton growth rate k without nutrient addition.\
The potential for nutrient limitation was assessed by comparing apparent phytoplankton growth rates k in nutrient amended (k_N) and nonamended (k_NoN) replicates using a paired t-test. If a significant difference was found (p below 0.05) between k_N and k_NoN, nutrient-amended growth rates (mu_N, d^-1) were also calculated as:\
mu_N = g + k_N. \
Otherwise, all k_N and k_NoN triplicate of replicate values were used to calculate both g and mu_0.\
When size fractionation at 10 µm was performed only on nutrient amended samples, growth rates reported on greater than 10 µm and less than 10 µm fractions in this study were nutrient-amended growth rates (mu_N) when nutrient limitation was observed. If no nutrient limitation was observed, mu_N obtained is equivalent to mu_0.\
The uncertainty of g estimates was quantified using the standard error of the slope fit from a linear regression between replicate k values and dilution levels. When the slope obtained was not significantly different from zero (p higher than 0.05), g was set to 0. \
Thus, the average k_N represented mu_N and the average k_NoN represented mu_0. \
A significant positive slope (i.e. higher growth in the WSW treatment than in the diluted) represents a violation of the method’s assumption. \
In such cases, g was reported as undetermined, and k in the undiluted bottles represented mu_N and mu_0. Uncertainties relative to mu_N and mu_0 were estimated from the standard deviations observed on k_N and k_NoN triplicate values.\
*Input: CRUSIE-chl-grazing-experiments-k-values.csv files*\
*Outputs: CRUISE-chla-grazing-experiments-rates.csv files.*\

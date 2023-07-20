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

The **chl-grazing_experiment_rates.m** Matlab script is used to compute phytoplankton growth rates (mu0), microzooplankton grazing rates (g), apparent growth rates in nonamended nutrient treatment (wsw NoN, kNoN) and phytoplankton growth rates in nutrient amended treatments (muN = g + kN). Associted errors (std) were estimated for all these rates. Note that kNoN and muN were estimated only when there was apparent nutrient limitation.\
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
The uncertainty of g estimates was quantified from the mean standard deviation of k_dil, k_NoN and k_N. \
g was set to 0 when negative g values not significantly different from 0 (t-test between k_dil and k_N) were obtained. \
Thus, the average k_N represented mu_N and the average k_NoN represented mu_0 (Murrell et al., 2002; Chen et al. 2009). \
Negative g values, significantly different from 0 (i.e. higher growth in the WSW treatment than in the diluted) represents a violation of the method's assumption. \ 
In such cases, g was reported as undetermined, and k in the undiluted bottles represented mu_N and mu_0. Uncertainties relative to mu_N and mu_0 were estimated from the standard deviations of the k_N and k_NoN triplicate values. \
Thus, the average k_N represented mu_N and the average k_NoN represented mu_0. \
A significant positive slope (i.e. higher growth in the WSW treatment than in the diluted) represents a violation of the method’s assumption. \
In such cases, g was reported as undetermined, and k in the undiluted bottles represented mu_N and mu_0. Uncertainties relative to mu_N and mu_0 were estimated from the standard deviations observed on k_N and k_NoN triplicate values.\
*Input: CRUSIE-chl-grazing-experiments-k-values.csv files*\
*Outputs: CRUISE-chla-grazing-experiments-rates.csv files.*\

The **chl_grazing_experiment_concatenate_order_rates.m** Matlab script is used to concatenation of the CRUISE rate table to get one final table with all the cruises together, and reorganization of the table by cruise and by inverse latitude (station nb).\
Concatenation of the 11 tables with rates.\
Rearragement of the table by cruise, from the oldest (EN608) to the most recent (EN687), by inverse latitude (station nb) and by size fraction (>0&<200, >10&<200, >0&<10 (from difference between >0&<200 and >10&<200 Chl-a conc), >0 (no mesh, en627) and >0&<10sf (size fractionated dilution experiments during en668).\
*Input: CRUSIE-chl-grazing-experiments-rates.csv files*\
*Outputs: NES-LTER-chla-grazing-experiments-rates.csv file*\

The **chl_grazing_experiment_rates_QC.m** Matlab file is used to Quality Check (QC) and rename of some values of the rate data based on the following criteria:
1) Change grazing rates < 0 (and g_std) to n/d and change grazing rates = NaN to n/d
2) Change muN = NaN (and mu_N_std) to n/n
3) Change mu0 = NaN (and mu_N_std) to n/d. Most of the occurence are from en608,en617 and en627  </> 10 um size fractions. </> 10 um size fractions only from nutrient amended samples, so only mu_N and g, no mu_0. Note that when no nutrient limitation for >0&<200 size fraction, mu_0 can be considered equal to mu_N. For these cruise and these size fraction mu_0 = n/d. Few other occurence (5) for other cruises for the >0&<10 size fraction, because no Chl-a data with QC = 1 for these samples
4) if temp_diff (temperature difference between sampling temperature and temperature in incubator <-4oC or > 4oC, iode_quality_flag (QC flag) = 3 (questionable)
5) For light level = 100% (en644), iode_quality_flag (QC flag) = 3 (questionable). It was too much light intensity for the experiment and the phytoplankton got "burned/fried"
6) if dilution (dilution level for the dilution experiment) > 0.4 (40%), iode_quality_flag (QC flag) = 3 (questionable). The optimal dilution level for the 2-points method is <40% (Morison and Menden-Deuer, 2017)
7) if Chlad10 (<10um) or Chlau10 (>10um) concentrations are < 0.02 mg m-3, the rates for these size fractions are considered questionable (iode_quality_flag (QC flag) = 3)
8) if Chlad10per (<10um) or Chlau10per (>10um) relative contribution to total Chl-a are < 0.02 (2%), the rates for these size fractions are considered questionable (iode_quality_flag (QC flag) = 3)\
*Input: NES-LTER-chla-grazing-experiments-rates.csv file*\
*Outputs: NES-LTER-chla-grazing-experiments-rates-qc.csv file*

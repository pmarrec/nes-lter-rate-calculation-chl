%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for the Quality Check (QC) and renaming  of some values 
% of the rate data based on the following criteria:
%
% 1) Change grazing rates < 0 (and g_std) to n/d
%
% 2) Change muN = NaN (and mu_N_std) to n/n
%
% 3) Change mu0 = NaN (and mu_N_std) to n/d. Most of the occurence are from 
% en608,en617 and en627  </> 10 um size fractions. 
% </> 10 um size fractions only from nutrient amended samples, 
% so only mu_N and g, no mu_0. Note that when no nutrient
% limitation for >0&<200 size fraction, mu_0 can be considered equal to
% mu_N. For these cruise and these size fraction mu_0 = n/d.
% Few other occurence (5) for other cruises for the >0&<10 size fraction, 
% because no Chl-a data with QC = 1 for these samples  
%
% 4) if temp_diff (temperature difference between sampling temperature and 
% temperature in incubator <-4oC or > 4oC, iode_quality_flag (QC flag) = 3
% (questionable)
%
% 5) For light level = 100% (en644), iode_quality_flag (QC flag) = 3
% (questionable). It was too much light intensity for the experiment and
% the phytoplankton got "burned/fried"
%
% 6) if dilution (dilution level for the dilution experiment) > 0.4 (40%), 
% iode_quality_flag (QC flag) = 3 (questionable). The optimal dilution
% level for the 2-points method is <40% (Morison and Menden-Deuer, 2017)
%
% 7) if Chlad10 (<10um) or Chlau10 (>10um) concentrations are < 0.02 mg m-3, the rates for
% these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
%
% 8) if Chlad10per (<10um) or Chlau10per (>10um) relative contribution to total Chl-a 
% are < 0.02 (2%), the rates for these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
%
%
% Input: NES-LTER-chla-grazing-experiments-rates.csv file
%
% Outputs: NES-LTER-chla-grazing-experiments-rates-qc.csv file.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 2/28/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\NES-LTER_Chla_Cleaning_Rates_Computation\';

% Find all the *.csv files
tablename=strcat(rep,'NES-LTER-chla-grazing-experiments-rates.csv');
T1=readtable(tablename);

%Create a new column "iode_quality_flag" with flag = 1
T1.iode_quality_flag=ones(height(T1),1);


% 1) Change grazing rates < 0 (and g_std) to n/d
a1=find(T1.g<0);
T1.g=num2cell(T1.g);
T1.g_std=num2cell(T1.g_std);
for n1=1:length(a1)
T1.g{a1(n1)}={'n/d'};
T1.g_std{a1(n1)}={'n/d'};
end

% 2) Change mu_N = NaN (an mu_N_std) to n/n
a2=find(isnan(T1.muN));
T1.muN=num2cell(T1.muN);
T1.muN_std=num2cell(T1.muN_std);
for n2=1:length(a2)
T1.muN{a2(n2)}={'n/n'};
T1.muN_std{a2(n2)}={'n/n'};
end

% 3) Change mu0 = NaN (and mu_N_std) to n/d. Most of the occurence are from 
% en608,en617 and en627  </> 10 um size fractions. 
% </> 10 um size fractions only from nutrient amended samples, 
% so only mu_N and g, no mu_0. Note that when no nutrient
% limitation for >0&<200 size fraction, mu_0 can be considered equal to
% mu_N. For these cruise and these size fraction mu_0 = n/d.
% Few other occurence (5) for other cruises for the >0&<10 size fraction, 
% because no Chl-a data with QC = 1 for these samples
a3=find(isnan(T1.mu0));
T1.mu0=num2cell(T1.mu0);
T1.mu0_std=num2cell(T1.mu0_std);
for n3=1:length(a3)
T1.mu0{a3(n3)}={'n/d'};
T1.mu0_std{a3(n3)}={'n/d'};
end

% 4) if temp_diff (temperature difference between sampling temperature and 
% temperature in incubator <-4oC or > 4oC, iode_quality_flag (QC flag) = 3
% (questionable)
dT=T1.temperature_sampling - T1.temperature_incubation_avg;
a4=(dT<-4)|(dT>4);
T1.iode_quality_flag(a4)=3;

% 5) For light level = 100% (en644), iode_quality_flag (QC flag) = 3
% (questionable). It was too much light intensity for the experiment and
% the phytoplankton got "burned/fried"
a5=strcmp(T1.light_level,'100%');
T1.iode_quality_flag(a5)=3;

% 6) if dilution (dilution level for the dilution experiment) > 0.4 (40%), 
% iode_quality_flag (QC flag) = 3 (questionable). The optimal dilution
% level for the 2-points method is <40% (Morison and Menden-Deuer, 2017)
a6=(T1.dilution>0.4);
T1.iode_quality_flag(a6)=3;

% 7) if Chlad10 (<10um) or Chlau10 (>10um) concentrations are < 0.02 mg m-3, the rates for
% these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
a71=(strcmp(T1.size_fraction,'>0&<10'))&(T1.Chlad10<0.02);
T1.iode_quality_flag(a71)=3;
a72=(strcmp(T1.size_fraction,'>10&<200'))&(T1.Chlau10<0.02);
T1.iode_quality_flag(a72)=3;

% 8) if Chlad10per (<10um) or Chlau10per (>10um) relative contribution to total Chl-a 
% are < 0.02 (2%), the rates for these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
a81=(strcmp(T1.size_fraction,'>0&<10'))&(T1.Chlad10per<0.02);
T1.iode_quality_flag(a81)=3;
a82=(strcmp(T1.size_fraction,'>10&<200'))&(T1.Chlau10per<0.02);
T1.iode_quality_flag(a82)=3;

%Save the new table csv file
newtablename=strcat(rep,'NES-LTER-chla-grazing-experiments-rates-qc.csv');
writetable(T1,newtablename);
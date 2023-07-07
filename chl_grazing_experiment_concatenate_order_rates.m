%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for concatenation of the CRUISE rate table to get one final
% table with all the cruises together, and reorganization of the table by
% cruise and by inverse latitude (station nb).
%
% 1) concatenation of the 11 tables with rates
% 2) rearragement of the table by cruise, from the oldest (EN608) to the most
% recent (EN687), by inverse latitude (station nb) and by size fraction (>0&<200, >10&<200,
% >0&<10 (from difference between >0&<200 and >10&<200 Chl-a conc),
% 
%
% Input: CRUSIE-chl-grazing-experiments-rates.csv files
%
% Outputs: NES-LTER-chla-grazing-experiments-rates.csv file.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 3/7/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\EDI_Growth_Grazing\DataPackage_GFF_and_10um\';
% Set the directory where the input rates data are
rep1 = strcat(rep,'chl-grazing-experiment-rates\');

% Find all the *.csv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);% List all files of interest in the directory

all_data = cell(1,numel(list));
for n1 = 1:numel(list)
    tablename=strcat(rep1,list(n1).name);
    all_data{n1} = readtable(tablename);
end



% concatenate all the tables into one big table, 
T1=cat(1,all_data{:});

% sort with inverse latitude order (from L1 to L11)
T2=sortrows(T1,'latitude','descend');

% sort according to the size fraction order
%set up the desired order in terms of cruise
D2={'en608';'en617';'en627';'en644';'en649';'en655';'en657';'en661';...
    'en668';'at46';'ar66b';'en687'};
[X2,Y2] = ismember(T2.cruise,D2);
[~,Z2] = sort(Y2);
T3 = T2(Z2,:);

%get only 4 decimal digits for latitude and longitude in table
T3.latitude=round(T3.latitude,4);
T3.longitude=round(T3.longitude,4);
%get only 2 decimal digits for the other numbers in table (except
%incubation_tank
T3.depth=round(T3.depth,2);
T3.temperature_sampling=round(T3.temperature_sampling,2);
T3.temperature_incubation_avg=round(T3.temperature_incubation_avg,2);
T3.temperature_incubation_std=round(T3.temperature_incubation_std,2);
T3.Chla=round(T3.Chla,2);
T3.Chlad10=round(T3.Chlad10,2);
T3.Chlau10=round(T3.Chlau10,2);
T3.Chlad10per=round(T3.Chlad10per,2);
T3.Chlau10per=round(T3.Chlau10per,2);
T3.duration_incubation=round(T3.duration_incubation,2);
T3.dilution=round(T3.dilution,2);
T3.mu_0=round(T3.mu_0,2);
T3.mu_0_std=round(T3.mu_0_std,2);
T3.grazing=round(T3.grazing,2);
T3.grazing_std=round(T3.grazing_std,2);
T3.mu_N=round(T3.mu_N,2);
T3.mu_N_std=round(T3.mu_N_std,2);

% and write it to
% output_file:
newtablename=strcat(rep,'NES-LTER-chla-grazing-experiments-rates.csv');
writetable(T3,newtablename);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for concatenation of the CRUISE rate table to get one final
% table with all the cruises together, and reorganization of the table by
% cruise and by inverse latitude (station nb)
%
% 1) concatenation of the 11 tables with rates
% 2) rearragement of the table by cruise, from the oldest (EN608) to the most
% recent (EN687), by inverse latitude (station nb) and by size fraction (>0&<200, >10&<200,
% >0&<10 (from difference between >0&<200 and >10&<200 Chl-a conc),
% >0 (no mesh, en627) and >0&<10sf (size fractionated dilution experiments
% during en668)
%
% Input: CRUSIE-chl-grazing-experiments-rates.csv files
%
% Outputs: NES-LTER-chla-grazing-experiments-rates.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 1/25/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\NES-LTER_Chla_Cleaning_Rates_Computation\';
% Set the directory where the input rates data are
rep1 = strcat(rep,'chl-grazing-experiment-rates\');

% Find all the *.cnv files
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

% sort according to the cruise desired order
%set up the desired order in terms of cruise
D3={'en608';'en617';'en627';'en644';'en649';'en655';'en657';'en661';...
    'en668';'at46';'ar66b';'en687'};
[X3,Y3] = ismember(T3.cruise,D3);
[~,Z3] = sort(Y3);
T4 = T3(Z3,:);

%get only 4 decimal digits for latitude and longitude in table
T4.latitude=round(T4.latitude,4);
T4.longitude=round(T4.longitude,4);
%get only 2 decimal digits for the other numbers in table (except
%incubation_tank
T4.depth=round(T4.depth,2);
T4.temperature_sampling=round(T4.temperature_sampling,2);
T4.temperature_incubation_avg=round(T4.temperature_incubation_avg,2);
T4.temperature_incubation_std=round(T4.temperature_incubation_std,2);
T4.Chla=round(T4.Chla,2);
T4.Chlad10=round(T4.Chlad10,2);
T4.Chlau10=round(T4.Chlau10,2);
T4.Chlad10per=round(T4.Chlad10per,2);
T4.Chlau10per=round(T4.Chlau10per,2);
T4.duration_incubation=round(T4.duration_incubation,2);
T4.dilution=round(T4.dilution,2);
T4.mu0=round(T4.mu0,2);
T4.mu0_std=round(T4.mu0_std,2);
T4.g=round(T4.g,2);
T4.g_std=round(T4.g_std,2);
T4.muN=round(T4.muN,2);
T4.muN_std=round(T4.muN_std,2);

% and write it to
% output_file:
newtablename=strcat(rep,'NES-LTER-chla-grazing-experiments-rates.csv');
writetable(T4,newtablename);



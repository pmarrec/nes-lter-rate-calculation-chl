%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for removing the single quotes from the cast, niskin and second niskin columns.
%
% Input: NES-LTER-chla-grazing-experiments-rates-qc.csv files
%
% Outputs: NES-LTER-chla-grazing-experiments-rates-qc-noquote.csv file.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 3/7/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\EDI_Growth_Grazing\';

% Find all the *.csv files
tablename=strcat(rep,'NES-LTER-chla-grazing-experiments-rates-qc.csv');
T1=readtable(tablename);

% Erase the extra " ' " in T1.cast and T1.niskin
T1.cast=erase(T1.cast,"'");
T1.niskin=erase(T1.niskin,"'");
T1.niskin_second_cast=erase(T1.niskin_second_cast,"'");

%Save the new table csv file
newtablename=strcat(rep,'NES-LTER-chla-grazing-experiments-rates-qc-noquote.csv');
writetable(T1,newtablename);

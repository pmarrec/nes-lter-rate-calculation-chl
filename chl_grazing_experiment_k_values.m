%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for calculation of the k values (apparent growth rate)
% obtained during dilution (grazing) experiments.
%
% Create a new table for each cruise gathering all the caluclated k-values
% (apparent groth rates) for each tretment and flter type.
% The duration of each incuabtion is calculated from the date_time_utc_end
% and the date_time_utc_start for each experiment.
% All the triplicate T0 Chl-a conc are calculated (only the one with
% iode_quality_flag = 1).
% T0 WSW mean Chla (Total = Chla, >10um = Chlau10, <10um = Chlad10) are 
% reported in the new table, in addition to the % of Chl-a </>10um 
% (Chlad10per and Chlau10per). 
% Chlad10 = Chla - Chlau10. If Chlad10 < 0, Chlad10 = 0, Chlad10per = 0% 
% and Chlau10per = 100%. 
% k are caluclated from each TF values. Up to 6 k values are then obatined
% for each treatment (dilution/nutrient/light) and filter type (>0&<200,
% >10&<200, >0 and >0&<10 ).
% Based on >0&200 (GFF) and >10&<200 (10um filters, u10 = up 10),
% >0&<10 (d10 = down 10) Chl-a values are calculated and then corresponding
% k values. For each triplicate, Chl-a d10 triplcate values are calculated
% as the diffrenece between the mean Chl-a value on >0&<200 and individual
% triplicate Chl-a values of >10&<200. Only dat with QC = 1 are considered.
% If Chl-a d10 <0, then Chl-a conc and k = NaN.
%
%
% Input: CRUSIE-chl-grazing-experiments-clean.csv files
%
% Outputs: CRUISE-chla-grazing-experiments-k-values.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 5/26/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\My Drive\NES-LTER_Chla_Cleaning_Rates_Computation\';
% Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-clean\');
% Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-k-values\');

% Find all the *.cnv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);% List all files of interest in the directory

for n1=1:numel(list)
    % load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);

    % Create a table for each cruise to store the individual k-values
    % For each cast/niskin, 6 rows if replicate bottles a and b, 9 rows if
    % replicate bottles a, b anc c (from AR66).
    % The number of colums depends on the number of treatments (nutrient, light)
    % and on the filter types (>0&<200, >0&<200, but also >0 for EN627-L11-B and
    % >0&<10 for EN668.

    % Only consider the TF values
    b0=strcmp(T1.T0_TF,'TF');
    % Get the number of station/depth (cast/niskin) sampled during the
    % cruise
    [c1,C1,C2]=findgroups(T1.cast(b0),T1.niskin(b0));
    % Get the number of replicate bottles (if a0=2, 6 rows per station/depth, if
    % a0=3, 9 columns per station/depth
    a0=unique(T1.replicate_bottle(b0));
    if length(a0)==2
        nrow=6;
    elseif length(a0)==3
        nrow=9;
    end

    % Start with a table T2 with the good nb of rows and with only 4 variables:
    % crusie, cast. niskin, and dilution = observed dilution based on GFF
    % samples. Extra columns for the k values for each treatment will be
    % added later
    T2=table('Size',[length(C1)*nrow 25],'VariableTypes',...
        {'string','string','string','string',...
        'string','double','string','string','string',...
        'double','double','double','double',...
        'double','double','double','string','string',...
        'double','double','double','double','double',...
        'double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_other_method',...
        'nearest_station','distance','date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','temperature_sampling','incubation_tank',...
        'temperature_incubation_avg','temperature_incubation_std', ...
        'light_level_HL','light_level_LL','Chla','Chlad10','Chlau10','Chlad10per','Chlau10per',...
        'duration_incubation','dilution'});


    % Erase the extra " ' " in T1.cast and T1.niskin
    T1.cast=erase(T1.cast,"'");
    T1.niskin=erase(T1.niskin,"'");

    % Identify each unique cast
    a1=unique(T1.cast);

    % Set up a loop counter for indexing the unique cast/depth
    cnt1 = 0;

    % Find the rows corresponding to the corresponding cast
    for n2=1:length(a1)
        b1=strcmp(T1.cast,a1(n2));
        % for each cast, identify the unique sampling depth
        a2=unique(T1.niskin(b1));
        for n3=1:length(a2)

            cnt1 = cnt1 +1;
            %Get the index of all the values from a given depth
            b2=b1 & strcmp(T1.niskin,a2(n3));
            % Duration of the incubation
            B2 = find(b2, 1, 'first');%find the first occurence of b2=1
            T2_start=(cnt1*nrow-(nrow-1));%Define where to store the first value for this given cast/depth
            T2_end=(cnt1*nrow);%Define where to store the last value for this given cast/depth

            Tinc=datenum(T1.date_time_utc_end(B2),'yyyy-mm-dd hh:MM:ss')-datenum(T1.date_time_utc_start(B2),'yyyy-mm-dd hh:MM:ss');
            % Correpsonding cruise/cast/niskin/date_time...
            T2.cruise(T2_start:T2_end)=T1.cruise(B2);
            T2.cast(T2_start:T2_end)=T1.cast(B2);
            T2.niskin(T2_start:T2_end)=T1.niskin(B2);
            T2.niskin_other_method(T2_start:T2_end)=T1.niskin_other_method(B2);
            T2.nearest_station(T2_start:T2_end)=T1.nearest_station(B2);
            T2.distance(T2_start:T2_end)=T1.distance(B2);
            T2.date_time_utc_sampling(T2_start:T2_end)=T1.date_time_utc_sampling(B2);
            T2.date_time_utc_start(T2_start:T2_end)=T1.date_time_utc_start(B2);
            T2.date_time_utc_end(T2_start:T2_end)=T1.date_time_utc_end(B2);
            T2.latitude(T2_start:T2_end)=T1.latitude(B2);
            T2.longitude(T2_start:T2_end)=T1.longitude(B2);
            T2.depth(T2_start:T2_end)=T1.depth(B2);
            T2.temperature_sampling(T2_start:T2_end)=T1.temperature_sampling(B2);
            T2.incubation_tank(T2_start:T2_end)=T1.incubation_tank(B2);
            T2.temperature_incubation_avg(T2_start:T2_end)=T1.temperature_incubation_avg(B2);
            T2.temperature_incubation_std(T2_start:T2_end)=T1.temperature_incubation_std(B2);
            T2.duration_incubation(T2_start:T2_end)=Tinc;

            %Get the light level for HL and LL
            d1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'));
            D1=find(d1==1);
            if ~isempty(D1)
                T2.light_level_HL(T2_start:T2_end)=T1.light_level(D1);%Get the HL light level
            else
                T2.light_level_HL(T2_start:T2_end)='NA';%Nan if no HL
            end
            d2=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%'));
            D2=find(d2==1);
            if ~isempty(D2)
                T2.light_level_LL(T2_start:T2_end)=T1.light_level(D2);%Get the HL light level
            else
                T2.light_level_LL(T2_start:T2_end)='NA';%Nan if no HL
            end

            %Get T0 WSW Chl-a Total, Chl-a<10um and Chl-a>10um, %<10um and
            % %>10um
            e1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            CHLA=mean(T1.chl(e1));
            T2.Chla(T2_start:T2_end)=CHLA;
            e2=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            CHLAu10=mean(T1.chl(e2));
            T2.Chlau10(T2_start:T2_end)=CHLAu10;
            CHLAd10=CHLA-CHLAu10;
            if CHLAd10<0
                T2.Chlad10(T2_start:T2_end)=0;
                T2.Chlau10per(T2_start:T2_end)=100;
                T2.Chlad10per(T2_start:T2_end)=0;
            else
                T2.Chlad10(T2_start:T2_end)=CHLAd10;
                T2.Chlau10per(T2_start:T2_end)=CHLAu10/CHLA;
                T2.Chlad10per(T2_start:T2_end)=CHLAd10/CHLA;
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0&<200 filters (GFF)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mean Chl-a T0 dil >0&<200
            % Identify all values obtained with >0&<200 filters at T0 dil
            % with a iode_quality_flag = 1
            c1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
            chl_T0_dil=mean(T1.chl(c1));
            % Mean Chl-a T0 wsw >0&<200
            c2=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            chl_T0_wsw=mean(T1.chl(c2));
            % Dilution level
            T2.dilution(T2_start:T2_end)=chl_T0_dil/chl_T0_wsw;
            % k-values dil High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_HL(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_HL(T2_start+C4-1)=nan;
            else
                T2.k_dil_HL(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                if length(a0)==2%When only 2 replicate bottles
                    T2.k_wsw_NoN_HL(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_HL(T2_start+C4-1)=nan;
                elseif length(a0)==3%When only 3 replicate bottles
                    T2.k_wsw_NoN_HL(T2_start:T2_end-3)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_HL(T2_start+C4-1)=nan;
                    T2.k_wsw_NoN_HL(T2_start+6:T2_end)=nan;
                end
            else
                T2.k_wsw_NoN_HL(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw + N High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_HL(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_HL(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_HL(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values dil Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_LL(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_LL(T2_start+C4-1)=nan;
            else
                T2.k_dil_LL(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                if length(a0)==2%When only 2 replicate bottles
                    T2.k_wsw_NoN_LL(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_LL(T2_start+C4-1)=nan;
                elseif length(a0)==3%When only 3 replicate bottles
                    T2.k_wsw_NoN_LL(T2_start:T2_end-3)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_LL(T2_start+C4-1)=nan;
                    T2.k_wsw_NoN_LL(T2_start+6:T2_end)=nan;
                end
            else
                T2.k_wsw_NoN_LL(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw N Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_LL(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_LL(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_LL(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >10&<200 filters (10um) k_u10
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mean Chl-a T0 dil >10&<200
            % Identify all values obtained with >10&<200 filters at T0 dil
            % with a iode_quality_flag = 1
            c1=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
            chl_T0_dil_u10=mean(T1.chl(c1));
            % Mean Chl-a T0 wsw >10&<200
            c2=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            chl_T0_wsw_u10=mean(T1.chl(c2));
            % k-values dil High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_HL_u10(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil_u10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_HL_u10(T2_start+C4-1)=nan;
            else
                T2.k_dil_HL_u10(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                if length(a0)==2%When only 2 replicate bottles
                    T2.k_wsw_NoN_HL_u10(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_u10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_HL_u10(T2_start+C4-1)=nan;
                elseif length(a0)==3%When only 3 replicate bottles
                    T2.k_wsw_NoN_HL_u10(T2_start:T2_end-3)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_u10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_HL_u10(T2_start+C4-1)=nan;
                    T2.k_wsw_NoN_HL_u10(T2_start+6:T2_end)=nan;
                end
            else
                T2.k_wsw_NoN_HL_u10(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw + N High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_HL_u10(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_u10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_HL_u10(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_HL_u10(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values dil Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_LL_u10(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil_u10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_LL_u10(T2_start+C4-1)=nan;
            else
                T2.k_dil_LL_u10(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                if length(a0)==2%When only 2 replicate bottles
                    T2.k_wsw_NoN_LL_u10(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_u10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_LL_u10(T2_start+C4-1)=nan;
                elseif length(a0)==3%When only 3 replicate bottles
                    T2.k_wsw_NoN_LL_u10(T2_start:T2_end-3)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_u10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_LL_u10(T2_start+C4-1)=nan;
                    T2.k_wsw_NoN_LL_u10(T2_start+6:T2_end)=nan;
                end
            else
                T2.k_wsw_NoN_LL_u10(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw N Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_LL_u10(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_u10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_LL_u10(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_LL_u10(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % k_d10 - Apparent growth rates of <10um size fraction from
            % the difference between >0&<200 and >10&<200 filters
            % Different from >0&<10 Chl-a conc obatined during en668
            % k_d10_sf (for size fraction)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mean Chl-a T0 dil <10um (d10)
            chl_T0_dil_d10=chl_T0_dil-chl_T0_dil_u10;
            if chl_T0_dil_d10<0
                chl_T0_dil_d10=nan;
            end
            % Mean Chl-a T0 dil <10um (d10)
            chl_T0_wsw_d10=chl_T0_wsw-chl_T0_wsw_u10;
            if chl_T0_wsw_d10<0
                chl_T0_wsw_d10=nan;
            end
            % k-values dil High Light (65% or 100% for EN644)
            % Rewrite it for replicate bottle a, b and c if c
            c31= b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'));
            c32= b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) & T1.iode_quality_flag==1;
            chl_TF=mean(T1.chl(c32));
            C3=find(c31==1);
            if ~isempty(C3)
                chl_TF_d10=chl_TF-T1.chl(C3);
                b3=chl_TF_d10<0;
                chl_TF_d10(b3)=nan;
                T2.k_dil_HL_d10(T2_start:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_HL_d10(T2_start+C4-1)=nan;
            else
                T2.k_dil_HL_d10(T2_start:T2_end)=nan;
            end
            clear c31 c32 chl_TF C3 C4 chl_TF_d10 b3
            % k-values wsw NoN High Light (65% or 100% for EN644)
            c31= b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'))&...
                strcmp(T1.nutrient_treatment,'NoN');
            c32= b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) & T1.iode_quality_flag==1 &...
                strcmp(T1.nutrient_treatment,'NoN');
            chl_TF=mean(T1.chl(c32));
            C3=find(c31==1);
            if ~isempty(C3)
                chl_TF_d10=chl_TF-T1.chl(C3);
                b3=chl_TF_d10<0;
                chl_TF_d10(b3)=nan;
                if length(a0)==2%When only 2 replicate bottles
                    T2.k_wsw_NoN_HL_d10(T2_start:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_HL_d10(T2_start+C4-1)=nan;
                elseif length(a0)==3%When only 3 replicate bottles
                    T2.k_wsw_NoN_HL_d10(T2_start:T2_end-3)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_HL_d10(T2_start+C4-1)=nan;
                    T2.k_wsw_NoN_HL_d10(T2_start+6:T2_end)=nan;
                end
            else
                T2.k_wsw_NoN_HL_d10(T2_start:T2_end)=nan;
            end
            clear c31 c32 chl_TF C3 C4 chl_TF_d10 b3
            % k-values wsw + N High Light (65% or 100% for EN644)
            c31= b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'N');
            c32= b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) & T1.iode_quality_flag==1 &...
                strcmp(T1.nutrient_treatment,'N');
            chl_TF=mean(T1.chl(c32));
            C3=find(c31==1);
            if ~isempty(C3)
                chl_TF_d10=chl_TF-T1.chl(C3);
                b3=chl_TF_d10<0;
                chl_TF_d10(b3)=nan;
                T2.k_wsw_N_HL_d10(T2_start:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_HL_d10(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_HL_d10(T2_start:T2_end)=nan;
            end
            clear c31 c32 chl_TF C3 C4 chl_TF_d10 b3
            % k-values dil Low Light (30% or 15% or 5% or 3% or 1%)
            c31= b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%'));
            c32= b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) & ...
                T1.iode_quality_flag==1;
            chl_TF=mean(T1.chl(c32));
            C3=find(c31==1);
            if ~isempty(C3)
                chl_TF_d10=chl_TF-T1.chl(C3);
                b3=chl_TF_d10<0;
                chl_TF_d10(b3)=nan;
                T2.k_dil_LL_d10(T2_start:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_LL_d10(T2_start+C4-1)=nan;
            else
                T2.k_dil_LL_d10(T2_start:T2_end)=nan;
            end
            clear c31 c32 chl_TF C3 C4 chl_TF_d10 b3
            % k-values wsw NoN Low Light (30% or 15% or 5% or 3% or 1%)
            c31= b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            c32= b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                T1.iode_quality_flag==1 & strcmp(T1.nutrient_treatment,'NoN');
            chl_TF=mean(T1.chl(c32));
            C3=find(c31==1);
            if ~isempty(C3)
                chl_TF_d10=chl_TF-T1.chl(C3);
                b3=chl_TF_d10<0;
                chl_TF_d10(b3)=nan;
                if length(a0)==2%When only 2 replicate bottles
                    T2.k_wsw_NoN_LL_d10(T2_start:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_LL_d10(T2_start+C4-1)=nan;
                elseif length(a0)==3%When only 3 replicate bottles
                    T2.k_wsw_NoN_LL_d10(T2_start:T2_end-3)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                    C4=find(T1.iode_quality_flag(C3)==3);
                    T2.k_wsw_NoN_LL_d10(T2_start+C4-1)=nan;
                    T2.k_wsw_NoN_LL_d10(T2_start+6:T2_end)=nan;
                end
            else
                T2.k_wsw_NoN_LL_d10(T2_start:T2_end)=nan;
            end
            clear c31 c32 chl_TF C3 C4 chl_TF_d10 b3
            % k-values wsw N Low Light (30% or 15% or 5% or 3% or 1%)
            c31= b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'N');
            c32= b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                T1.iode_quality_flag==1 & strcmp(T1.nutrient_treatment,'N');
            chl_TF=mean(T1.chl(c32));
            C3=find(c31==1);
            if ~isempty(C3)
                chl_TF_d10=chl_TF-T1.chl(C3);
                b3=chl_TF_d10<0;
                chl_TF_d10(b3)=nan;
                T2.k_wsw_N_LL_d10(T2_start:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_LL_d10(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_LL_d10(T2_start:T2_end)=nan;
            end
            clear c31 c32 chl_TF C3 C4 chl_TF_d10 b3

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0&<10 filters (10um size fractionation EN668)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mean Chl-a T0 dil >0&<10
            % Identify all values obtained with >0&<10 filters at T0 dil
            % with a iode_quality_flag = 1
            c1=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
            chl_T0_dil_d10_sf=mean(T1.chl(c1));
            % Mean Chl-a T0 wsw >0&<10
            c2=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            chl_T0_wsw_d10_sf=mean(T1.chl(c2));
            % k-values dil High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_HL_d10_sf(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil_d10_sf);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_HL_d10_sf(T2_start+C4-1)=nan;
            else
                T2.k_dil_HL_d10_sf(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_NoN_HL_d10_sf(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_d10_sf);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_NoN_HL_d10_sf(T2_start+C4-1)=nan;
            else
                T2.k_wsw_NoN_HL_d10_sf(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw + N High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_HL_d10_sf(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_d10_sf);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_HL_d10_sf(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_HL_d10_sf(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values dil Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_LL_d10_sf(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil_d10_sf);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_LL_d10_sf(T2_start+C4-1)=nan;
            else
                T2.k_dil_LL_d10_sf(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_NoN_LL_d10_sf(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_d10_sf);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_NoN_LL_d10_sf(T2_start+C4-1)=nan;
            else
                T2.k_wsw_NoN_LL_d10_sf(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw N Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_LL_d10_sf(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_d10_sf);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_LL_d10_sf(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_LL_d10_sf(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0 filters (no 200um screening EN627 L11-B)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mean Chl-a T0 dil >0
            % Identify all values obtained with >0&<10 filters at T0 dil
            % with a iode_quality_flag = 1
            c1=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
            chl_T0_dil_no_mesh=mean(T1.chl(c1));
            % Mean Chl-a T0 wsw >0
            c2=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            chl_T0_wsw_no_mesh=mean(T1.chl(c2));
            % k-values dil High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_HL_no_mesh(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil_no_mesh);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_HL_no_mesh(T2_start+C4-1)=nan;
            else
                T2.k_dil_HL_no_mesh(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_NoN_HL_no_mesh(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_no_mesh);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_NoN_HL_no_mesh(T2_start+C4-1)=nan;
            else
                T2.k_wsw_NoN_HL_no_mesh(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw + N High Light (65% or 100% for EN644)
            c3=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'65%')|...
                strcmp(T1.light_level,'100%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_HL_no_mesh(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_no_mesh);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_HL_no_mesh(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_HL_no_mesh(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values dil Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'dil') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%'));
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_dil_LL_no_mesh(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil_no_mesh);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_dil_LL_no_mesh(T2_start+C4-1)=nan;
            else
                T2.k_dil_LL_no_mesh(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw NoN Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'NoN');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_NoN_LL_no_mesh(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_no_mesh);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_NoN_LL_no_mesh(T2_start+C4-1)=nan;
            else
                T2.k_wsw_NoN_LL_no_mesh(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4
            % k-values wsw N Low Light (30% or 15% or 5% or 3% or 1%)
            c3=b2 & strcmp(T1.filter_size,'>0') & strcmp(T1.T0_TF,'TF') ...
                & strcmp(T1.dilution,'wsw') & (strcmp(T1.light_level,'30%')|...
                strcmp(T1.light_level,'15%')|strcmp(T1.light_level,'5%')|...
                strcmp(T1.light_level,'3%')|strcmp(T1.light_level,'1%')) &...
                strcmp(T1.nutrient_treatment,'N');
            C3=find(c3==1);
            if ~isempty(C3)
                T2.k_wsw_N_LL_no_mesh(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw_no_mesh);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_N_LL_no_mesh(T2_start+C4-1)=nan;
            else
                T2.k_wsw_N_LL_no_mesh(T2_start:T2_end)=nan;
            end
            clear c3 C3 C4

        end
    end

    clear cnt1 Tinc



    % Make sure cast and niskin are in text format in the table
    T2.cast=strcat(T2.cast,"'");
    T2.niskin=strcat(T2.niskin,"'");

    % Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    % cruise
    newname=strrep(list(n1).name,'clean','k-values');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path
    writetable(T2,newtablename)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for calculation of the k values (apparent growth rate)
% obtained during dilution (grazing) experiments.
%
% Create a new table for each cruise gathering all the caluclated k-values
% (apparent groth rates) for each tretment and flter type.
% The duration of each incuabtion is calculated from the date_time_utc_end
% and the date_time_utc_start for each experiment.
% Only Triplicate (or x6 for some cruises) with iode_quality_flag = 1) T0
% Chl-a conc.
% T0 WSW mean Chla (Total = Chla, >10um = Chlau10, <10um = Chlad10) are
% reported in the new table, in addition to the % of Chl-a </>10um
% (Chlad10per and Chlau10per).
% Chlad10 = Chla - Chlau10. If Chlad10 < 0, Chlad10 = 0, Chlad10per = 0%
% and Chlau10per = 100%.
% k are caluclated from each TF values.
% Up to 6 k (9 k values from AR66b) values are then obatined
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
% 3/7/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\EDI_Growth_Grazing\';
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
    % cruise and the number of different treatment (light/filter size)
    [c1,C1,C2,C3,C4]=findgroups(T1.cast(b0),T1.niskin(b0),T1.filter_size(b0),T1.light_level(b0));

    % Get the number of replicate bottles (if a0=2, 6 rows per station/depth, if
    % a0=3, 9 columns per station/depth
    a0=unique(T1.replicate_bottle(b0));
    if length(a0)==2
        nrow=6;
    elseif length(a0)==3
        nrow=9;
    end

    %find the number of >10&<200 "groups" in C3 to then add the requested
    %amount of rows for >0&<10 k values
    f0=strcmp(C3,'>10&<200');
    F0=sum(f0);

    % Start with a table T2 with the good nb of rows
    % Nb of rows = number of of different station/depth/filter size/light
    % treatment + F0 nb of >0&<10  multiply by nrows (6 or 9 depending on
    % cruises)
    T2=table('Size',[(length(C1)+F0)*nrow 29],'VariableTypes',...
        {'string','string','string','string',...
        'string','string','string','string',...
        'double','double','double','double',...
        'double','double','double','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_second_cast','station',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','temperature_sampling','incubation_tank',...
        'temperature_incubation_avg','temperature_incubation_std', ...
        'light_level','size_fraction','replicate_bottle','replicate_chl',...
        'Chla','Chlad10','Chlau10','Chlad10per','Chlau10per',...
        'duration_incubation','dilution','k_dil','k_wsw_NoN','k_wsw_N'});


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
            b2=b1 & strcmp(T1.niskin,a2(n3));

            %Get T0 WSW Chl-a Total, Chl-a<10um and Chl-a>10um, %<10um and
            % %>10um
            e1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            CHLA=mean(T1.chl(e1));

            e2=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
                & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
            CHLAu10=mean(T1.chl(e2));

            %for each sampling depth, indetify the unique light levels
            a3=unique(T1.light_level(b2));

            for n4=1:(length(a3)-1) %(-1) to don't consider Light leveles = NA (e.g. at T0)
                %for eahc type of filter (+1 if >10&<200 roo count for
                %>0&<200)
                b3=b2 & strcmp(T1.light_level,a3(n4));
                %for each sampling depth, indetify the unique filter sizes
                a4=unique(T1.filter_size(b3));

                for n5=1:length(a4)

                    cnt1 = cnt1 +1;
                    %Get the index of all the values from a given size
                    %fraction
                    b4=b3 & strcmp(T1.filter_size,a4(n5));
                    % Duration of the incubation
                    B4 = find(b4, 1, 'first');%find the first occurence of b4=1
                    Tinc=datenum(T1.date_time_utc_end(B4),'yyyy-mm-dd hh:MM:ss')-datenum(T1.date_time_utc_start(B4),'yyyy-mm-dd hh:MM:ss');
                    %Define where to store the first value for this given cast/depth
                    T2_start=(cnt1*nrow-(nrow-1));
                    %Define where to store the last value for this given cast/depth
                    T2_end=(cnt1*nrow);
                    % Correpsonding cruise/cast/niskin/date_time...
                    T2.cruise(T2_start:T2_end)=T1.cruise(B4);
                    T2.cast(T2_start:T2_end)=T1.cast(B4);
                    T2.niskin(T2_start:T2_end)=T1.niskin(B4);
                    T2.niskin_second_cast(T2_start:T2_end)=T1.niskin_second_cast(B4);
                    T2.station(T2_start:T2_end)=T1.station(B4);
                    T2.date_time_utc_sampling(T2_start:T2_end)=T1.date_time_utc_sampling(B4);
                    T2.date_time_utc_start(T2_start:T2_end)=T1.date_time_utc_start(B4);
                    T2.date_time_utc_end(T2_start:T2_end)=T1.date_time_utc_end(B4);
                    T2.latitude(T2_start:T2_end)=T1.latitude(B4);
                    T2.longitude(T2_start:T2_end)=T1.longitude(B4);
                    T2.depth(T2_start:T2_end)=T1.depth(B4);
                    T2.temperature_sampling(T2_start:T2_end)=T1.temperature_sampling(B4);
                    T2.incubation_tank(T2_start:T2_end)=T1.incubation_tank(B4);
                    T2.temperature_incubation_avg(T2_start:T2_end)=T1.temperature_incubation_avg(B4);
                    T2.temperature_incubation_std(T2_start:T2_end)=T1.temperature_incubation_std(B4);
                    T2.duration_incubation(T2_start:T2_end)=Tinc;
                    if nrow==6
                        T2.replicate_bottle(T2_start:T2_end)={'a';'a';'a';'b';'b';'b'};
                        T2.replicate_chl(T2_start:T2_end)={'a';'b';'c';'a';'b';'c'};
                    elseif nrow==9
                        T2.replicate_bottle(T2_start:T2_end)={'a';'a';'a';'b';'b';'b';'c';'c';'c'};
                        T2.replicate_chl(T2_start:T2_end)={'a';'b';'c';'a';'b';'c';'a';'b';'c'};
                    end
                    T2.size_fraction(T2_start:T2_end)=T1.filter_size(B4);
                    T2.light_level(T2_start:T2_end)=T1.light_level(B4);

                    T2.Chla(T2_start:T2_end)=CHLA;
                    T2.Chlau10(T2_start:T2_end)=CHLAu10;
                    CHLAd10=CHLA-CHLAu10;
                    %correct if Chla '>0&<200' by difference < 0 and then
                    %compute the % of Chl-a in each size fraction
                    if CHLAd10<0
                        T2.Chlad10(T2_start:T2_end)=0;
                        T2.Chlau10per(T2_start:T2_end)=1;
                        T2.Chlad10per(T2_start:T2_end)=0;
                    else
                        T2.Chlad10(T2_start:T2_end)=CHLAd10;
                        T2.Chlau10per(T2_start:T2_end)=CHLAu10/CHLA;
                        T2.Chlad10per(T2_start:T2_end)=CHLAd10/CHLA;
                    end
                    %if Chla '>10&<200' > Chla '>0&<200'
                    if CHLAu10>CHLA
                        CHLAu10=CHLA;
                        T2.Chlau10(T2_start:T2_end)=CHLAu10;
                    end

                    %Get the dilution level from >0&<200 filters at T0
                    %This dilution level also used for size fractions
                    % Identify all values obtained with >0&<200 filters at T0 dil
                    % with a iode_quality_flag = 1
                    c1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
                    CHL_T0_dil=mean(T1.chl(c1));
                    % Mean Chl-a T0 wsw >0&<200
                    c2=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
                    CHL_T0_wsw=mean(T1.chl(c2));
                    % Dilution level
                    T2.dilution(T2_start:T2_end)=CHL_T0_dil/CHL_T0_wsw;

                    clear c1 c2


                    % Define the T0 Chl-a concentrations to use depending on
                    % the filter type
                    % Identify all values obtained with these filters at T0 dil
                    % with a iode_quality_flag = 1
                    c1=b2 & strcmp(T1.filter_size,a4(n5)) & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
                    chl_T0_dil=mean(T1.chl(c1));
                    % Mean Chl-a T0 wsw >0&<200
                    c2=b2 & strcmp(T1.filter_size,a4(n5)) & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
                    chl_T0_wsw=mean(T1.chl(c2));

                    % k-values dil
                    c3= b0 & b4 & strcmp(T1.dilution,'dil');
                    C3=find(c3==1);
                    if ~isempty(C3)
                        T2.k_dil(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_dil(T2_start+C4-1)=nan;
                    else
                        T2.k_dil(T2_start:T2_end)=nan;
                    end
                    clear c3 C3 C4
                    % k-values wsw
                    c3=b0 & b4 & strcmp(T1.dilution,'wsw') & strcmp(T1.nutrient_treatment,'NoN');
                    C3=find(c3==1);
                    if ~isempty(C3)
                        if length(a0)==2%When only 2 replicate bottles
                            T2.k_wsw_NoN(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                            C4=find(T1.iode_quality_flag(C3)==3);
                            T2.k_wsw_NoN(T2_start+C4-1)=nan;
                        elseif length(a0)==3%When only 3 replicate bottles
                            T2.k_wsw_NoN(T2_start:T2_end-3)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                            C4=find(T1.iode_quality_flag(C3)==3);
                            T2.k_wsw_NoN(T2_start+C4-1)=nan;
                            T2.k_wsw_NoN(T2_start+6:T2_end)=nan;
                        end
                    else
                        T2.k_wsw_NoN(T2_start:T2_end)=nan;
                    end
                    clear c3 C3 C4
                    % k-values wsw + N High Light (65% or 100% for EN644)
                    c3=b0 & b4 & strcmp(T1.dilution,'wsw') & strcmp(T1.nutrient_treatment,'N');
                    C3=find(c3==1);
                    if ~isempty(C3)
                        T2.k_wsw_N(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_wsw_N(T2_start+C4-1)=nan;
                    else
                        T2.k_wsw_N(T2_start:T2_end)=nan;
                    end
                    clear c3 C3 C4

                end

                %if >10&<200 data, consider also >0&<10 size fraction,
                %which is obtained by the difference between >0&<200 and
                %>10&<200.
                %DIFFERENT FROM THE '>0&<10' FILTER_SIZE USED DURING EN668
                %the difference between >0&<200 and
                %>10&<200 size fraction will be named >0&<10 size_fraction
                %instead of the previous filter type (cf end of the script)
                if sum(ismember(a4,'>10&<200'))==1
                    cnt1 = cnt1 +1;

                    %Get the index of all the values from the >0&<200
                    %filter_size
                    b41=b3 & strcmp(T1.filter_size,'>0&<200');
                    % Duration of the incubation
                    B41 = find(b41, 1, 'first');%find the first occurence of b4=1
                    Tinc=datenum(T1.date_time_utc_end(B41),'yyyy-mm-dd hh:MM:ss')-datenum(T1.date_time_utc_start(B41),'yyyy-mm-dd hh:MM:ss');
                    %Define where to store the first value for this given cast/depth
                    T2_start=(cnt1*nrow-(nrow-1));
                    %Define where to store the last value for this given cast/depth
                    T2_end=(cnt1*nrow);
                    % Correpsonding cruise/cast/niskin/date_time...
                    T2.cruise(T2_start:T2_end)=T1.cruise(B41);
                    T2.cast(T2_start:T2_end)=T1.cast(B41);
                    T2.niskin(T2_start:T2_end)=T1.niskin(B41);
                    T2.niskin_second_cast(T2_start:T2_end)=T1.niskin_second_cast(B41);
                    T2.station(T2_start:T2_end)=T1.station(B41);
                    T2.date_time_utc_sampling(T2_start:T2_end)=T1.date_time_utc_sampling(B41);
                    T2.date_time_utc_start(T2_start:T2_end)=T1.date_time_utc_start(B41);
                    T2.date_time_utc_end(T2_start:T2_end)=T1.date_time_utc_end(B41);
                    T2.latitude(T2_start:T2_end)=T1.latitude(B41);
                    T2.longitude(T2_start:T2_end)=T1.longitude(B41);
                    T2.depth(T2_start:T2_end)=T1.depth(B41);
                    T2.temperature_sampling(T2_start:T2_end)=T1.temperature_sampling(B41);
                    T2.incubation_tank(T2_start:T2_end)=T1.incubation_tank(B41);
                    T2.temperature_incubation_avg(T2_start:T2_end)=T1.temperature_incubation_avg(B41);
                    T2.temperature_incubation_std(T2_start:T2_end)=T1.temperature_incubation_std(B41);
                    T2.duration_incubation(T2_start:T2_end)=Tinc;
                    T2.size_fraction(T2_start:T2_end)='<10';%Name '<10' for the moment, will be rename at the end of the script
                    T2.light_level(T2_start:T2_end)=T1.light_level(B41);
                    if nrow==6
                        T2.replicate_bottle(T2_start:T2_end)={'a';'a';'a';'b';'b';'b'};
                        T2.replicate_chl(T2_start:T2_end)={'a';'b';'c';'a';'b';'c'};
                    elseif nrow==9
                        T2.replicate_bottle(T2_start:T2_end)={'a';'a';'a';'b';'b';'b';'c';'c';'c'};
                        T2.replicate_chl(T2_start:T2_end)={'a';'b';'c';'a';'b';'c';'a';'b';'c'};
                    end

                    T2.Chla(T2_start:T2_end)=CHLA;
                    T2.Chlau10(T2_start:T2_end)=CHLAu10;
                    CHLAd10=CHLA-CHLAu10;
                    %correct if Chla '>0&<200' by difference < 0 and then
                    %compute the % of Chl-a in each size fraction
                    if CHLAd10<0
                        T2.Chlad10(T2_start:T2_end)=0;
                        T2.Chlau10per(T2_start:T2_end)=1;
                        T2.Chlad10per(T2_start:T2_end)=0;
                    else
                        T2.Chlad10(T2_start:T2_end)=CHLAd10;
                        T2.Chlau10per(T2_start:T2_end)=CHLAu10/CHLA;
                        T2.Chlad10per(T2_start:T2_end)=CHLAd10/CHLA;
                    end
                    %if Chla '>10&<200' > Chla '>0&<200'
                    if CHLAu10>CHLA
                        CHLAu10=CHLA;
                        T2.Chlau10(T2_start:T2_end)=CHLAu10;
                    end

                    %Get the dilution level from >0&<200 filters at T0
                    %This dilution level also used for size fractions
                    % Identify all values obtained with >0&<200 filters at T0 dil
                    % with a iode_quality_flag = 1
                    c1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
                    CHL_T0_dil=mean(T1.chl(c1));
                    % Mean Chl-a T0 wsw >0&<200
                    c2=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
                    CHL_T0_wsw=mean(T1.chl(c2));
                    % Dilution level
                    T2.dilution(T2_start:T2_end)=CHL_T0_dil/CHL_T0_wsw;

                    clear c1 c2


                    % Identify all values obtained with >0&<200 filters at T0 dil
                    % with a iode_quality_flag = 1
                    c1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
                    chl_T0_dil=mean(T1.chl(c1));
                    % Mean Chl-a T0 wsw >0&<200
                    c2=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
                    chl_T0_wsw=mean(T1.chl(c2));

                    % Identify all values obtained with >0&<200 filters at T0 dil
                    % with a iode_quality_flag = 1
                    c1=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
                    chl_T0_dil_u10=mean(T1.chl(c1));
                    % Mean Chl-a T0 wsw >0&<200
                    c2=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
                        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
                    chl_T0_wsw_u10=mean(T1.chl(c2));

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

                    % k-values dil for bottle a
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'a');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'a');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        T2.k_dil(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_dil(T2_start+C4-1)=nan;
                    else
                        T2.k_dil(T2_start:T2_start+2)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

                    % k-values dil for bottle b
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'b');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'b');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        T2.k_dil(T2_start+3:T2_start+5)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_dil(T2_start+C4-1)=nan;
                    else
                        T2.k_dil(T2_start+3:T2_start+5)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

                    % k-values dil for bottle c
                    if nrow==9
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'c');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'c');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        T2.k_dil(T2_start+6:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_dil(T2_start+C4-1)=nan;
                    else
                        T2.k_dil(T2_start+6:T2_end)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10
                    end


                    % k-values wsw NoN for bottle a
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
                        & strcmp(T1.nutrient_treatment,'NoN');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
                        & strcmp(T1.nutrient_treatment,'NoN');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        if length(a0)==2%When only 2 replicate bottles
                            T2.k_wsw_NoN(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                            C4=find(T1.iode_quality_flag(C3)==3);
                            T2.k_wsw_NoN(T2_start+C4-1)=nan;
                        elseif length(a0)==3%When only 3 replicate bottles
                            T2.k_wsw_NoN(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                            C4=find(T1.iode_quality_flag(C3)==3);
                            T2.k_wsw_NoN(T2_start+C4-1)=nan;
                            T2.k_wsw_NoN(T2_start+6:T2_end)=nan;
                        end
                    else
                        T2.k_wsw_NoN(T2_start:T2_end)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

                    % k-values wsw NoN for bottle b
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
                        & strcmp(T1.nutrient_treatment,'NoN');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
                        & strcmp(T1.nutrient_treatment,'NoN');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        if length(a0)==2%When only 2 replicate bottles
                            T2.k_wsw_NoN(T2_start+3:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                            C4=find(T1.iode_quality_flag(C3)==3);
                            T2.k_wsw_NoN(T2_start+C4-1)=nan;
                        elseif length(a0)==3%When only 3 replicate bottles
                            T2.k_wsw_NoN(T2_start+3:T2_end-3)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                            C4=find(T1.iode_quality_flag(C3)==3);
                            T2.k_wsw_NoN(T2_start+C4-1)=nan;
                            T2.k_wsw_NoN(T2_start+3:T2_start+5)=nan;
                        end
                    else
                        T2.k_wsw_NoN(T2_start+3:T2_end)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

                    % k-values wsw + N for bottle a
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
                        & strcmp(T1.nutrient_treatment,'N');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
                        & strcmp(T1.nutrient_treatment,'N');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        T2.k_wsw_N(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_wsw_N(T2_start+C4-1)=nan;
                    else
                        T2.k_wsw_N(T2_start:T2_start+2)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

                    % k-values wsw + N for bottle b
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
                        & strcmp(T1.nutrient_treatment,'N');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
                        & strcmp(T1.nutrient_treatment,'N');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        T2.k_wsw_N(T2_start+3:T2_start+5)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_wsw_N(T2_start+C4-1)=nan;
                    else
                        T2.k_wsw_N(T2_start+3:T2_start+5)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

                    % k-values wsw + N for bottle c
                    if nrow==9
                    c31= b0 & b3 & strcmp(T1.filter_size,'>10&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'c')...
                        & strcmp(T1.nutrient_treatment,'N');
                    c32= b0 & b3 & strcmp(T1.filter_size,'>0&<200') ...
                        & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'c')...
                        & strcmp(T1.nutrient_treatment,'N');
                    chl_TF=mean(T1.chl(c32));
                    C3=find(c31==1);
                    if ~isempty(C3)
                        chl_TF_d10=chl_TF-T1.chl(C3);
                        c4=chl_TF_d10<0;
                        chl_TF_d10(c4)=nan;
                        T2.k_wsw_N(T2_start+6:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                        C4=find(T1.iode_quality_flag(C3)==3);
                        T2.k_wsw_N(T2_start+C4-1)=nan;
                    else
                        T2.k_wsw_N(T2_start+6:T2_end)=nan;
                    end
                    clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10
                    end


                end

            end

        end

    end

    % Make sure cast and niskin are in text format in the table
    T2.cast=strcat(T2.cast,"'");
    T2.niskin=strcat(T2.niskin,"'");

    %Change the name of the size fractions
    T2.size_fraction=strrep(T2.size_fraction,'>0&<10','>0&<10sf');
    s1=strcmp(T2.size_fraction,'<10');
    T2.size_fraction(s1)=strrep(T2.size_fraction(s1),'<10','>0&<10');
    
    % Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    % cruise
    newname=strrep(list(n1).name,'clean','k-values');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path
    writetable(T2,newtablename)

end
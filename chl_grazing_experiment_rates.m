%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for calculation of the rates obtained during dilution (grazing)
% experiments.
%
% Create a new table for each cruise gathering all the caluclated rates
% for each tretment and flter type from the k (apparent growth rates) values
% and the associated data and metadata.
% Phytoplankton growth rates and protist grazing rates were estimated from
% 24 h changes in Chl-a concentration. For each incubation bottle,
% the apparent growth rates (k, d^-1) were calculated as: k=1⁄t×ln(C_t⁄C_0)
% where t is the incubation time (d) and C_t and C_0 the final and
% initial Chl-a concentration (µg L^-1), respectively.
% Protist grazing rates (g, d^-1) were estimated with the equation:
% g=((k_d-k_N))⁄((1-x))
% where k_d and k_N are the apparent growth rates in 20WSW and WSW nutrient
% amended treatments, respectively, and x is the achieved fraction of WSW
% in the diluted treatment calculated from T0 Chl-a in 20WSW and WSW. Accordingly,
% the instantaneous, or in situ, growth rate (mu_0, d^-1) was estimated as:
% mu_0=g+k_NoN
% where k_NoN is apparent phytoplankton growth rate k without nutrient addition.
% The potential for nutrient limitation was assessed by comparing apparent
% phytoplankton growth rates k in nutrient amended (k_N) and nonamended
% (k_NoN) replicates using a paired t-test. If a significant difference was
% found (p below 0.05) between k_N and k_NoN, nutrient-amended growth rates
% (mu_N, d^-1) were also calculated as mu_N = g + k_N. Otherwise, all k_N
% and k_NoN triplicate of replicate values were used to calculate both g and mu_0.
% When size fractionation at 10 µm was performed only on nutrient amended samples,
% growth rates reported on greater than 10 µm and less than 10 µm fractions
% in this study were nutrient-amended growth rates (mu_N) when nutrient
% limitation was observed. If no nutrient limitation was observed,
% mu_N obtained is equivalent to mu_0.
% The uncertainty of g estimates was quantified using the standard error
% of the slope fit from a linear regression between replicate k values
% and dilution levels. When the slope obtained was not significantly
% different from zero (p higher than 0.05), g was set to 0. Thus,
% the average k_N represented mu_N and the average k_NoN represented mu_0.
% A significant positive slope (i.e. higher growth in the WSW treatment
% than in the diluted) represents a violation of the method’s assumption.
% In such cases, g was reported as undetermined, and k in the undiluted
% bottles represented mu_N and mu_0. Uncertainties relative to mu_N and
% mu_0 were estimated from the standard deviations observed on k_N and
% k_NoN triplicate values.
%
%
% Input: CRUSIE-chl-grazing-experiments-k-values.csv files
%
% Outputs: CRUISE-chla-grazing-experiments-rates.csv files.
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
% Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-k-values\');
% Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-rates\');

% Find all the *.csv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);% List all files of interest in the directory


for n1=1:numel(list)
    % load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);

    % Create a table for each cruise to store the rates
    % for each cast/niskin
    % The number of colums depends on the number of treatments (nutrient, light)
    % and on the filter types (>0&<200, >0&<200, but also >0 for EN627-L11-B and
    % >0&<10 for EN668 and will be implemented for each cruise

    % Get the number of station/depth (cast/niskin) sampled during the
    % cruise
    [c1,C1,C2]=findgroups(T1.cast,T1.niskin);

    % Table T2 for with the good nb of rows and columns for rates
    T2=table('Size',[length(C1) 30],'VariableTypes',...
        {'string','string','string','string','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double',...
        'string','string',...
        'double','double','double','double','double',...
        'double','double',...
        'double','double','double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_second_cast','station',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','temperature_sampling','incubation_tank',...
        'temperature_incubation_avg','temperature_incubation_std',...
        'light_level','size_fraction','Chla','Chlad10','Chlau10','Chlad10per','Chlau10per',...
        'duration_incubation','dilution',...
        'mu_0','mu_0_std','grazing','grazing_std','mu_N','mu_N_std'});


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
            %for each sampling depth, indetify the unique light levels
            a3=unique(T1.light_level(b2));

            for n4=1:length(a3) 
                b3=b2 & strcmp(T1.light_level,a3(n4));
                %for each sampling depth, indetify the unique filter sizes
                a4=unique(T1.size_fraction(b3));

                for n5=1:length(a4)

                    cnt1 = cnt1 +1;
                    %Get the index of all the values from a given size
                    %fraction
                    b4=b3 & strcmp(T1.size_fraction,a4(n5));
                    % Duration of the incubation
                    B4 = find(b4, 1, 'first');%find the first occurence of b4=1

                    % Correpsonding cruise/cast/niskin/date_time...
                    T2.cruise(cnt1)=T1.cruise(B4);
                    T2.cast(cnt1)=T1.cast(B4);
                    T2.niskin(cnt1)=T1.niskin(B4);
                    T2.niskin_second_cast(cnt1)=T1.niskin_second_cast(B4);
                    T2.station(cnt1)=T1.station(B4);
                    T2.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B4);
                    T2.date_time_utc_start(cnt1)=T1.date_time_utc_start(B4);
                    T2.date_time_utc_end(cnt1)=T1.date_time_utc_end(B4);
                    T2.latitude(cnt1)=T1.latitude(B4);
                    T2.longitude(cnt1)=T1.longitude(B4);
                    T2.depth(cnt1)=T1.depth(B4);
                    T2.temperature_sampling(cnt1)=T1.temperature_sampling(B4);
                    T2.incubation_tank(cnt1)=T1.incubation_tank(B4);
                    T2.temperature_incubation_avg(cnt1)=T1.temperature_incubation_avg(B4);
                    T2.temperature_incubation_std(cnt1)=T1.temperature_incubation_std(B4);
                    T2.duration_incubation(cnt1)=T1.duration_incubation(B4);
                    T2.dilution(cnt1)=T1.dilution(B4);
                    T2.size_fraction(cnt1)=T1.size_fraction(B4);
                    T2.light_level(cnt1)=T1.light_level(B4);
                    T2.Chla(cnt1)=T1.Chla(B4);
                    T2.Chlau10(cnt1)=T1.Chlau10(B4);
                    T2.Chlad10(cnt1)=T1.Chlad10(B4);
                    T2.Chlau10per(cnt1)=T1.Chlau10per(B4);
                    T2.Chlad10per(cnt1)=T1.Chlad10per(B4);


                    %Create a vector with the dilution values
                    if length(T1.dilution(b4))==6%Replicate botles A and B
                        d1=[T1.dilution(b4);ones(12,1)];%When no nutrient limitation
                        d2=[T1.dilution(b4);ones(6,1)];%When nutrient limitation
                    else%Replicate bottles A, B and C
                        d1=[T1.dilution(b4);ones(18,1)];%When no nutrient limitation
                        d2=[T1.dilution(b4);ones(9,1)];%When nutrient limitation
                    end

                    k_dil=T1.k_dil(b4);
                    k_dil_nan=isnan(k_dil);
                    k_wsw_NoN=T1.k_wsw_NoN(b4);
                    k_wsw_NoN_nan=isnan(k_wsw_NoN);
                    k_wsw_N=T1.k_wsw_N(b4);
                    k_wsw_N_nan=isnan(k_wsw_N);

                    if sum(k_dil_nan)<length(k_dil) && ...
                            (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                        %Test if Nutrient limitation
                        [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','equal');%Need to define the different argument to consider in the paired t-test

                        if h==0%No nutrient Limitation

                             k=[k_dil;k_wsw_NoN;k_wsw_N];
                            g=(mean(k_dil,'omitnan')-mean([k_wsw_NoN;k_wsw_N],'omitnan'))/(1-T1.dilution(B4));

                            if g<0
                                %mu is then equal to k_100
                                mu=mean([k_wsw_NoN;k_wsw_N],'omitnan');
                                
                                %test if g is significantly different from
                                %0 with a ttest between k_dil and k_100
                                [h1,p1]=ttest2(k_dil,[k_wsw_NoN;k_wsw_N],'Tail','both','Vartype','equal');
                                
                                if h1==0% No signiicant diference between k_dil and k_100
                                    T2.grazing(cnt1)=0;%set g = 0
                                else% Significant difference between k_dil and k_100
                                    T2.grazing(cnt1)=g;%keep the <0 values, will be considered as non determined in the rate QC
                                end

                            else%g is positive
                                T2.grazing(cnt1)=g;
                                %mu is then equal to k_100 + g
                                mu=mean([k_wsw_NoN;k_wsw_N],'omitnan')+g;
                            end
                            
                            T2.mu_0(cnt1)=mu;
                            mu_stdev=max(std(k_dil,'omitnan'),std([k_wsw_NoN;k_wsw_N],'omitnan'));
                            %mu_stdev=mdl.Coefficients{1,2};
                            T2.mu_0_std(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                            g_stdev=mean([std(k_dil,'omitnan'),std([k_wsw_NoN;k_wsw_N],'omitnan')],'omitnan');
                            T2.grazing_std(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                            %No-nutrietn limitation, mu_N = NaN
                            T2.mu_N(cnt1)=nan;
                            T2.mu_N_std(cnt1)=nan;

                            clear k mdl mu g mu_stdev g_stdev

                        else %Nutirent Limited


                            %N ammended samples are used to compute muN and g
                            k=[k_dil;k_wsw_N];
                            g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B4));
                            
                            if g<0
                                %mu_N = k_wsw_N
                                muN=mean(k_wsw_N,'omitnan');
                                %test if g is significantly different from
                                %0 with a ttest between k_dil and k_100
                                [h1,p1]=ttest2(k_dil,k_wsw_N,'Tail','both','Vartype','equal');
                                
                                if h1==0% No signiicant diference between k_dil and k_100
                                    T2.grazing(cnt1)=0;%set g = 0
                                else% Significant difference between k_dil and k_100
                                    T2.grazing(cnt1)=g;%keep the <0 values, will be considered as non determined in the rate QC
                                end
                                
                            else%g is positive
                                T2.grazing(cnt1)=g;
                                muN=mean(k_wsw_N,'omitnan')+g;
                            end
                            T2.mu_N(cnt1)=muN;% mu = y-intercept = g + average k(1)
                            muN_stdev=max(std(k_dil,'omitnan'),std(k_wsw_N,'omitnan'));
                            T2.mu_N_std(cnt1)=muN_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                            g_stdev=mean([std(k_dil,'omitnan'),std(k_wsw_N,'omitnan')],'omitnan');
                            T2.grazing_std(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                            clear d mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                            %Computation of muNoN (in-situ growth rate) from k(1)NoN
                            %and g calculated from N amended samples

                            kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                            
                            if g<0
                                muNoN=kNoN;
                            else
                                muNoN=kNoN+g;
                            end

                            T2.mu_0(cnt1)=muNoN;
                            T2.mu_0_std(cnt1)=max(std(k_dil,'omitnan'),std(k_wsw_NoN,'omitnan'));%Find a better way to estimate StdDev on muNoN
                        end

                        clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

                    else%Nan Values if only NaN values for k
                        T2.mu_0(cnt1)=nan;
                        T2.mu_0_std(cnt1)=nan;
                        T2.grazing(cnt1)=nan;
                        T2.grazing_std(cnt1)=nan;
                        T2.mu_N(cnt1)=nan;
                        T2.mu_N_std(cnt1)=nan;

                    end

                end

            end

        end

    end

    % Make sure cast and niskin are in text format in the table
    T2.cast=strcat(T2.cast,"'");
    T2.niskin=strcat(T2.niskin,"'");


    % Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    % cruise
    newname2=strrep(list(n1).name,'k-values','rates');%Replace raw by clean
    newtablename2=strcat(rep2,newname2);%New tablename and path
    writetable(T2,newtablename2)

end
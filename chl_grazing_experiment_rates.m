%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for calculation of the rates obtained during dilution (grazing)
% experiments.
%
% Create a new table for each cruise gathering all the caluclated rates
% for each tretment and flter type from the k (apparent growth rates) values
% and the dilution levels.
% mu0, mu0_std:
% g, g_std:
% kNoN, kNoN_std:
% muN, muN_std:
% h_ttest, p_ttest:
%
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
% 6/9/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

% Set the directory where we work
rep = 'C:\Users\pierr\Desktop\NES-LTER_Chla_Cleaning_Rates_Computation\';
% Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-k-values\');
% Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-rates\');

% Find all the *.cnv files
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

    % Table T2 for with the good nb of rows and columns for rates from
    % >0&<200 size fraction (GFF)
    T2=table('Size',[length(C1) 28],'VariableTypes',...
        {'string','string','string','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_other_method',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','duration_incubation','dilution',...
        'mu0_HL','mu0_std_HL','g_HL','g_std_HL','kNoN_HL','kNoN_std_HL','muN_HL','muN_std_HL',...
        'mu0_LL','mu0_std_LL','g_LL','g_std_LL','kNoN_LL','kNoN_std_LL','muN_LL','muN_std_LL'});

    % Table T3 for with the good nb of rows and columns for rates from
    % >10&<200 size fraction (10um filters, u10, k_u10um)
    T3=table('Size',[length(C1) 28],'VariableTypes',...
        {'string','string','string','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_other_method',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','duration_incubation','dilution',...
        'mu0_HL','mu0_std_HL','g_HL','g_std_HL','kNoN_HL','kNoN_std_HL','muN_HL','muN_std_HL',...
        'mu0_LL','mu0_std_LL','g_LL','g_std_LL','kNoN_LL','kNoN_std_LL','muN_LL','muN_std_LL'});

    % Table T4 for with the good nb of rows and columns for rates from
    % >0&<10 size fraction (GFF - 10um, d10, k_d10um)
    T4=table('Size',[length(C1) 28],'VariableTypes',...
        {'string','string','string','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_other_method',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','duration_incubation','dilution',...
        'mu0_HL','mu0_std_HL','g_HL','g_std_HL','kNoN_HL','kNoN_std_HL','muN_HL','muN_std_HL',...
        'mu0_LL','mu0_std_LL','g_LL','g_std_LL','kNoN_LL','kNoN_std_LL','muN_LL','muN_std_LL'});

    % Table T5 for with the good nb of rows and columns for rates from
    % >0&<10 size fraction (10um size fractionation EN668, k_10um_sf)
    T5=table('Size',[length(C1) 28],'VariableTypes',...
        {'string','string','string','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_other_method',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','duration_incubation','dilution',...
        'mu0_HL','mu0_std_HL','g_HL','g_std_HL','kNoN_HL','kNoN_std_HL','muN_HL','muN_std_HL',...
        'mu0_LL','mu0_std_LL','g_LL','g_std_LL','kNoN_LL','kNoN_std_LL','muN_LL','muN_std_LL'});

    % Table T6 for with the good nb of rows and columns for rates from
    % >0 size fraction (no 200um screening EN627 L11-B)
    T6=table('Size',[length(C1) 28],'VariableTypes',...
        {'string','string','string','string',...
        'string','string','string',...
        'double','double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double',...
        'double','double','double','double'},...
        'VariableNames',{'cruise','cast','niskin','niskin_other_method',...
        'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
        'latitude','longitude','depth','duration_incubation','dilution',...
        'mu0_HL','mu0_std_HL','g_HL','g_std_HL','kNoN_HL','kNoN_std_HL','muN_HL','muN_std_HL',...
        'mu0_LL','mu0_std_LL','g_LL','g_std_LL','kNoN_LL','kNoN_std_LL','muN_LL','muN_std_LL'});

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

            b2=b1 & strcmp(T1.niskin,a2(n3));
            B2 = find(b2, 1, 'first');%find the first occurence of b2=1

            %Create a vector with the dilution values
            if length(T1.dilution(b2))==6%Replicate botles A and B
                d1=[T1.dilution(b2);ones(12,1)];%When no nutrient limitation
                d2=[T1.dilution(b2);ones(6,1)];%When nutrient limitation
            else%Replicate bottles A, B and C
                d1=[T1.dilution(b2);ones(18,1)];%When no nutrient limitation
                d2=[T1.dilution(b2);ones(9,1)];%When nutrient limitation
            end



            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0&<200 filters (GFF)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %High Light (65% or 100% for EN644)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % Correpsonding cruise/cast/niskin
            T2.cruise(cnt1)=T1.cruise(B2);
            T2.cast(cnt1)=T1.cast(B2);
            T2.niskin(cnt1)=T1.niskin(B2);
            T2.niskin_other_method(cnt1)=T1.niskin_other_method(B2);
            T2.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B2);
            T2.date_time_utc_start(cnt1)=T1.date_time_utc_start(B2);
            T2.date_time_utc_end(cnt1)=T1.date_time_utc_end(B2);
            T2.duration_incubation(cnt1)=T1.duration_incubation(B2);
            T2.dilution(cnt1)=T1.dilution(B2);

            k_dil=T1.k_dil_HL(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_HL(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_HL(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation



                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T2.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T2.mu0_std_HL(cnt1)=mu_stdev;
                        g=0;
                        T2.g_HL(cnt1)=g;
                        g_stdev=0;
                        T2.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T2.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T2.mu0_std_HL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T2.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T2.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T2.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T2.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T2.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T2.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end


                    T2.kNoN_HL(cnt1)=nan;
                    T2.kNoN_std_HL(cnt1)=nan;
                    T2.muN_HL(cnt1)=nan;
                    T2.muN_std_HL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T2.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T2.muN_std_HL(cnt1)=muN_stdev;
                        g=0;
                        T2.g_HL(cnt1)=g;
                        g_stdev=0;
                        T2.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T2.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T2.muN_std_HL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T2.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T2.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T2.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T2.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T2.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T2.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear d mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T2.kNoN_HL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T2.kNoN_std_HL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T2.mu0_HL(cnt1)=muNoN;
                    T2.mu0_std_HL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN
                end
                clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

            else%Nan Values if only NaN values for k
                T2.mu0_HL(cnt1)=nan;
                T2.mu0_std_HL(cnt1)=nan;
                T2.g_HL(cnt1)=nan;
                T2.g_std_HL(cnt1)=nan;
                T2.kNoN_HL(cnt1)=nan;
                T2.kNoN_std_HL(cnt1)=nan;
                T2.muN_HL(cnt1)=nan;
                T2.muN_std_HL(cnt1)=nan;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %Low Light (30% or 15% or 5% or 3% or 1%)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            k_dil=T1.k_dil_LL(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_LL(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_LL(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);


                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T2.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T2.mu0_std_LL(cnt1)=mu_stdev;
                        g=0;
                        T2.g_LL(cnt1)=g;
                        g_stdev=0;
                        T2.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T2.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T2.mu0_std_LL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T2.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T2.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T2.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T2.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T2.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T2.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    T2.kNoN_LL(cnt1)=nan;
                    T2.kNoN_std_LL(cnt1)=nan;
                    T2.muN_LL(cnt1)=nan;
                    T2.muN_std_LL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T2.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T2.muN_std_LL(cnt1)=muN_stdev;
                        g=0;
                        T2.g_LL(cnt1)=g;
                        g_stdev=0;
                        T2.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T2.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T2.muN_std_LL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T2.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T2.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T2.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T2.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T2.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T2.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T2.kNoN_LL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T2.kNoN_std_LL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T2.mu0_LL(cnt1)=muNoN;
                    T2.mu0_std_LL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN

                    clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

                end

            else%Nan Values if only NaN values for k
                T2.mu0_LL(cnt1)=nan;
                T2.mu0_std_LL(cnt1)=nan;
                T2.g_LL(cnt1)=nan;
                T2.g_std_LL(cnt1)=nan;
                T2.kNoN_LL(cnt1)=nan;
                T2.kNoN_std_LL(cnt1)=nan;
                T2.muN_LL(cnt1)=nan;
                T2.muN_std_LL(cnt1)=nan;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >10&<200 filters (GFF)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %High Light (65% or 100% for EN644)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % Correpsonding cruise/cast/niskin
            T3.cruise(cnt1)=T1.cruise(B2);
            T3.cast(cnt1)=T1.cast(B2);
            T3.niskin(cnt1)=T1.niskin(B2);
            T3.niskin_other_method(cnt1)=T1.niskin_other_method(B2);
            T3.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B2);
            T3.date_time_utc_start(cnt1)=T1.date_time_utc_start(B2);
            T3.date_time_utc_end(cnt1)=T1.date_time_utc_end(B2);
            T3.duration_incubation(cnt1)=T1.duration_incubation(B2);
            T3.dilution(cnt1)=T1.dilution(B2);

            k_dil=T1.k_dil_HL_u10um(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_HL_u10um(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_HL_u10um(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T3.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T3.mu0_std_HL(cnt1)=mu_stdev;
                        g=0;
                        T3.g_HL(cnt1)=g;
                        g_stdev=0;
                        T3.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T3.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T3.mu0_std_HL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T3.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T3.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T3.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T3.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T3.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T3.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end


                    T3.kNoN_HL(cnt1)=nan;
                    T3.kNoN_std_HL(cnt1)=nan;
                    T3.muN_HL(cnt1)=nan;
                    T3.muN_std_HL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T3.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T3.muN_std_HL(cnt1)=muN_stdev;
                        g=0;
                        T3.g_HL(cnt1)=g;
                        g_stdev=0;
                        T3.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T3.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T3.muN_std_HL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T3.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T3.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T3.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T3.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T3.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T3.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T3.kNoN_HL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T3.kNoN_std_HL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T3.mu0_HL(cnt1)=muNoN;
                    T3.mu0_std_HL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN
                end
                clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

            else%Nan Values if only NaN values for k
                T3.mu0_HL(cnt1)=nan;
                T3.mu0_std_HL(cnt1)=nan;
                T3.g_HL(cnt1)=nan;
                T3.g_std_HL(cnt1)=nan;
                T3.kNoN_HL(cnt1)=nan;
                T3.kNoN_std_HL(cnt1)=nan;
                T3.muN_HL(cnt1)=nan;
                T3.muN_std_HL(cnt1)=nan;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %Low Light (30% or 15% or 5% or 3% or 1%)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            k_dil=T1.k_dil_LL_u10um(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_LL_u10um(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_LL_u10um(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);


                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T3.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T3.mu0_std_LL(cnt1)=mu_stdev;
                        g=0;
                        T3.g_LL(cnt1)=g;
                        g_stdev=0;
                        T3.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T3.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T3.mu0_std_LL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T3.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T3.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T3.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T3.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T3.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T3.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    T3.kNoN_LL(cnt1)=nan;
                    T3.kNoN_std_LL(cnt1)=nan;
                    T3.muN_LL(cnt1)=nan;
                    T3.muN_std_LL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T3.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T3.muN_std_LL(cnt1)=muN_stdev;
                        g=0;
                        T3.g_LL(cnt1)=g;
                        g_stdev=0;
                        T3.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T3.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T3.muN_std_LL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T3.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T3.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T3.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T3.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T3.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T3.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T3.kNoN_LL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T3.kNoN_std_LL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T3.mu0_LL(cnt1)=muNoN;
                    T3.mu0_std_LL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN

                    clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

                end

            else%Nan Values if only NaN values for k
                T3.mu0_LL(cnt1)=nan;
                T3.mu0_std_LL(cnt1)=nan;
                T3.g_LL(cnt1)=nan;
                T3.g_std_LL(cnt1)=nan;
                T3.T3kNoN_LL(cnt1)=nan;
                T3.kNoN_std_LL(cnt1)=nan;
                T3.muN_LL(cnt1)=nan;
                T3.muN_std_LL(cnt1)=nan;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0&<10 size fraction (GFF - 10um, d10, k_d10um)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %High Light (65% or 100% for EN644)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % Correpsonding cruise/cast/niskin
            T4.cruise(cnt1)=T1.cruise(B2);
            T4.cast(cnt1)=T1.cast(B2);
            T4.niskin(cnt1)=T1.niskin(B2);
            T4.niskin_other_method(cnt1)=T1.niskin_other_method(B2);
            T4.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B2);
            T4.date_time_utc_start(cnt1)=T1.date_time_utc_start(B2);
            T4.date_time_utc_end(cnt1)=T1.date_time_utc_end(B2);
            T4.duration_incubation(cnt1)=T1.duration_incubation(B2);
            T4.dilution(cnt1)=T1.dilution(B2);

            k_dil=T1.k_dil_HL_d10um(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_HL_d10um(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_HL_d10um(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T4.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T4.mu0_std_HL(cnt1)=mu_stdev;
                        g=0;
                        T4.g_HL(cnt1)=g;
                        g_stdev=0;
                        T4.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T4.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T4.mu0_std_HL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T4.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T4.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T4.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T4.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T4.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T4.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end


                    T4.kNoN_HL(cnt1)=nan;
                    T4.kNoN_std_HL(cnt1)=nan;
                    T4.muN_HL(cnt1)=nan;
                    T4.muN_std_HL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T4.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T4.muN_std_HL(cnt1)=muN_stdev;
                        g=0;
                        T4.g_HL(cnt1)=g;
                        g_stdev=0;
                        T4.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T4.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T4.muN_std_HL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T4.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T4.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T4.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T4.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T4.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T4.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T4.kNoN_HL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T4.kNoN_std_HL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T4.mu0_HL(cnt1)=muNoN;
                    T4.mu0_std_HL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN
                end
                clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

            else%Nan Values if only NaN values for k
                T4.mu0_HL(cnt1)=nan;
                T4.mu0_std_HL(cnt1)=nan;
                T4.g_HL(cnt1)=nan;
                T4.g_std_HL(cnt1)=nan;
                T4.kNoN_HL(cnt1)=nan;
                T4.kNoN_std_HL(cnt1)=nan;
                T4.muN_HL(cnt1)=nan;
                T4.muN_std_HL(cnt1)=nan;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %Low Light (30% or 15% or 5% or 3% or 1%)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            k_dil=T1.k_dil_LL_d10um(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_LL_d10um(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_LL_d10um(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);


                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T4.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T4.mu0_std_LL(cnt1)=mu_stdev;
                        g=0;
                        T4.g_LL(cnt1)=g;
                        g_stdev=0;
                        T4.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T4.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T4.mu0_std_LL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T4.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T4.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T4.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T4.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T4.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T4.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    T4.kNoN_LL(cnt1)=nan;
                    T4.kNoN_std_LL(cnt1)=nan;
                    T4.muN_LL(cnt1)=nan;
                    T4.muN_std_LL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T4.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T4.muN_std_LL(cnt1)=muN_stdev;
                        g=0;
                        T4.g_LL(cnt1)=g;
                        g_stdev=0;
                        T4.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T4.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T4.muN_std_LL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T4.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T4.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T4.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T4.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T4.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T4.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T4.kNoN_LL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T4.kNoN_std_LL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T4.mu0_LL(cnt1)=muNoN;
                    T4.mu0_std_LL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN

                    clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

                end

            else%Nan Values if only NaN values for k
                T4.mu0_LL(cnt1)=nan;
                T4.mu0_std_LL(cnt1)=nan;
                T4.g_LL(cnt1)=nan;
                T4.g_std_LL(cnt1)=nan;
                T4.kNoN_LL(cnt1)=nan;
                T4.kNoN_std_LL(cnt1)=nan;
                T4.muN_LL(cnt1)=nan;
                T4.muN_std_LL(cnt1)=nan;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0&<10 size fraction (10um size fractionation EN668, k_10um_sf)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %High Light (65% or 100% for EN644)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % Correpsonding cruise/cast/niskin
            T5.cruise(cnt1)=T1.cruise(B2);
            T5.cast(cnt1)=T1.cast(B2);
            T5.niskin(cnt1)=T1.niskin(B2);
            T5.niskin_other_method(cnt1)=T1.niskin_other_method(B2);
            T5.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B2);
            T5.date_time_utc_start(cnt1)=T1.date_time_utc_start(B2);
            T5.date_time_utc_end(cnt1)=T1.date_time_utc_end(B2);
            T5.duration_incubation(cnt1)=T1.duration_incubation(B2);
            T5.dilution(cnt1)=T1.dilution(B2);

            k_dil=T1.k_dil_HL_10um_sf(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_HL_10um_sf(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_HL_10um_sf(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T5.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T5.mu0_std_HL(cnt1)=mu_stdev;
                        g=0;
                        T5.g_HL(cnt1)=g;
                        g_stdev=0;
                        T5.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T5.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T5.mu0_std_HL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T5.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T5.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T5.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T5.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T5.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T5.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end


                    T5.kNoN_HL(cnt1)=nan;
                    T5.kNoN_std_HL(cnt1)=nan;
                    T5.muN_HL(cnt1)=nan;
                    T5.muN_std_HL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T5.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T5.muN_std_HL(cnt1)=muN_stdev;
                        g=0;
                        T5.g_HL(cnt1)=g;
                        g_stdev=0;
                        T5.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T5.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T5.muN_std_HL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T5.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T5.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T5.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T5.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T5.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T5.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T5.kNoN_HL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T5.kNoN_std_HL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T5.mu0_HL(cnt1)=muNoN;
                    T5.mu0_std_HL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN
                end
                clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

            else%Nan Values if only NaN values for k
                T5.mu0_HL(cnt1)=nan;
                T5.mu0_std_HL(cnt1)=nan;
                T5.g_HL(cnt1)=nan;
                T5.g_std_HL(cnt1)=nan;
                T5.kNoN_HL(cnt1)=nan;
                T5.kNoN_std_HL(cnt1)=nan;
                T5.muN_HL(cnt1)=nan;
                T5.muN_std_HL(cnt1)=nan;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %Low Light (30% or 15% or 5% or 3% or 1%)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            k_dil=T1.k_dil_LL_10um_sf(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_LL_10um_sf(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_LL_10um_sf(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);


                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T5.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T5.mu0_std_LL(cnt1)=mu_stdev;
                        g=0;
                        T5.g_LL(cnt1)=g;
                        g_stdev=0;
                        T5.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T5.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T5.mu0_std_LL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T5.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T5.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T5.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T5.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T5.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T5.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    T5.kNoN_LL(cnt1)=nan;
                    T5.kNoN_std_LL(cnt1)=nan;
                    T5.muN_LL(cnt1)=nan;
                    T5.muN_std_LL(cnt1)=nan;

                    clear k mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T5.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T5.muN_std_LL(cnt1)=muN_stdev;
                        g=0;
                        T5.g_LL(cnt1)=g;
                        g_stdev=0;
                        T5.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T5.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T5.muN_std_LL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T5.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T5.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T5.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T5.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T5.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T5.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T5.kNoN_LL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T5.kNoN_std_LL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T5.mu0_LL(cnt1)=muNoN;
                    T5.mu0_std_LL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN

                    clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

                end

            else%Nan Values if only NaN values for k
                T5.mu0_LL(cnt1)=nan;
                T5.mu0_std_LL(cnt1)=nan;
                T5.g_LL(cnt1)=nan;
                T5.g_std_LL(cnt1)=nan;
                T5.kNoN_LL(cnt1)=nan;
                T5.kNoN_std_LL(cnt1)=nan;
                T5.muN_LL(cnt1)=nan;
                T5.muN_std_LL(cnt1)=nan;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % >0 size fraction (no 200um screening EN627 L11-B)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %High Light (65% or 100% for EN644)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % Correpsonding cruise/cast/niskin
            T6.cruise(cnt1)=T1.cruise(B2);
            T6.cast(cnt1)=T1.cast(B2);
            T6.niskin(cnt1)=T1.niskin(B2);
            T6.niskin_other_method(cnt1)=T1.niskin_other_method(B2);
            T6.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B2);
            T6.date_time_utc_start(cnt1)=T1.date_time_utc_start(B2);
            T6.date_time_utc_end(cnt1)=T1.date_time_utc_end(B2);
            T6.duration_incubation(cnt1)=T1.duration_incubation(B2);
            T6.dilution(cnt1)=T1.dilution(B2);

            k_dil=T1.k_dil_HL_no_mesh(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_HL_no_mesh(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_HL_no_mesh(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T6.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T6.mu0_std_HL(cnt1)=mu_stdev;
                        g=0;
                        T6.g_HL(cnt1)=g;
                        g_stdev=0;
                        T6.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T6.mu0_HL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T6.mu0_std_HL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T6.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T6.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T6.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T6.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T6.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T6.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end


                    T6.kNoN_HL(cnt1)=nan;
                    T6.kNoN_std_HL(cnt1)=nan;
                    T6.muN_HL(cnt1)=nan;
                    T6.muN_std_HL(cnt1)=nan;

                    clear d mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T6.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T6.muN_std_HL(cnt1)=muN_stdev;
                        g=0;
                        T6.g_HL(cnt1)=g;
                        g_stdev=0;
                        T6.g_std_HL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T6.muN_HL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T6.muN_std_HL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T6.g_HL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T6.g_std_HL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T6.g_HL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T6.mu0_HL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T6.mu0_std_HL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T6.g_std_HL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T6.kNoN_HL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T6.kNoN_std_HL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T6.mu0_HL(cnt1)=muNoN;
                    T6.mu0_std_HL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN
                end
                clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

            else%Nan Values if only NaN values for k
                T6.mu0_HL(cnt1)=nan;
                T6.mu0_std_HL(cnt1)=nan;
                T6.g_HL(cnt1)=nan;
                T6.g_std_HL(cnt1)=nan;
                T6.kNoN_HL(cnt1)=nan;
                T6.kNoN_std_HL(cnt1)=nan;
                T6.muN_HL(cnt1)=nan;
                T6.muN_std_HL(cnt1)=nan;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %Low Light (30% or 15% or 5% or 3% or 1%)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            k_dil=T1.k_dil_LL_no_mesh(b2);
            k_dil_nan=isnan(k_dil);
            k_wsw_NoN=T1.k_wsw_NoN_LL_no_mesh(b2);
            k_wsw_NoN_nan=isnan(k_wsw_NoN);
            k_wsw_N=T1.k_wsw_N_LL_no_mesh(b2);
            k_wsw_N_nan=isnan(k_wsw_N);

            if sum(k_dil_nan)<length(k_dil) && ...
                    (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

                [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','unequal');%Need to define the different argument to consider in the paired t-test
                %                 ResultsWSW.h_ttest(n)=h;
                %                 ResultsWSW.p_ttest(n)=p;

                if h==0%No nutrient Limitation

                    k=[k_dil;k_wsw_NoN;k_wsw_N];
                    mdl=fitlm(d1,k);


                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        mu=mean(k,'omitnan');
                        T6.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k,'omitnan');
                        T6.mu0_std_LL(cnt1)=mu_stdev;
                        g=0;
                        T6.g_LL(cnt1)=g;
                        g_stdev=0;
                        T6.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        mu=mean(k(7:18),'omitnan');%Mean value of k(1)
                        T6.mu0_LL(cnt1)=mu;
                        mu_stdev=std(k(7:18),'omitnan');%StdDev value of k(1)
                        T6.mu0_std_LL(cnt1)=mu_stdev;
                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T6.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T6.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k(1:6),'omitnan')-mean(k(7:18),'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T6.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k(7:18),'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T6.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k(7:18),'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T6.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T6.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    T6.kNoN_LL(cnt1)=nan;
                    T6.kNoN_std_LL(cnt1)=nan;
                    T6.muN_LL(cnt1)=nan;
                    T6.muN_std_LL(cnt1)=nan;

                    clear d mdl mu g mu_stdev g_stdev

                else %Nutirent Limited


                    %N ammended samples are used to compute muN and g

                    k=[k_dil;k_wsw_N];
                    mdl=fitlm(d2,k);

                    if mdl.Coefficients{2,4}>0.05%if g=0 (pvalue(slope)>0.05)
                        muN=mean(k_wsw_N,'omitnan');
                        T6.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');
                        T6.muN_std_LL(cnt1)=muN_stdev;
                        g=0;
                        T6.g_LL(cnt1)=g;
                        g_stdev=0;
                        T6.g_std_LL(cnt1)=g_stdev;

                    elseif mdl.Coefficients{2,1}>0%if g<0 (slope>0)

                        muN=mean(k_wsw_N,'omitnan');%Mean value of k(1)
                        T6.muN_LL(cnt1)=muN;
                        muN_stdev=std(k_wsw_N,'omitnan');%StdDev value of k(1)
                        T6.muN_std_LL(cnt1)=muN_stdev;
                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};%We keep the negative g value in case
                        T6.g_LL(cnt1)=g;
                        g_stdev=mdl.Coefficients{2,2};%Same for g StdDev
                        T6.g_std_LL(cnt1)=g_stdev;

                    else %if g>0 (slope<0)

                        g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));
                        %g=-mdl.Coefficients{2,1};
                        T6.g_LL(cnt1)=g;% g = -slope = [k(d)-k(1)]/(1-d)
                        mu=mean(k_wsw_N,'omitnan')+g;
                        %mu=mdl.Coefficients{1,1};
                        T6.mu0_LL(cnt1)=mu;% mu = y-intercept = g + average k(1)
                        mu_stdev=std(k_wsw_N,'omitnan');
                        %mu_stdev=mdl.Coefficients{1,2};
                        T6.mu0_std_LL(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                        g_stdev=mdl.Coefficients{2,2};
                        T6.g_std_LL(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                    end

                    clear k mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                    %Computation of muNoN (in-situ growth rate) from k(1)NoN
                    %and g calculated from N amended samples

                    kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)
                    T6.kNoN_LL(cnt1)=kNoN;
                    kNoN_stdev=std(k_wsw_NoN,'omitnan');%StdDev value of k(1)NoN (k(d) not included)
                    T6.kNoN_std_LL(cnt1)=kNoN_stdev;

                    if g<0
                        muNoN=kNoN;
                    else
                        muNoN=kNoN+g;
                    end

                    T6.mu0_LL(cnt1)=muNoN;
                    T6.mu0_std_LL(cnt1)=kNoN_stdev;%Find a better way to estimate StdDev on muNoN

                    clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

                end

            else%Nan Values if only NaN values for k
                T6.mu0_LL(cnt1)=nan;
                T6.mu0_std_LL(cnt1)=nan;
                T6.g_LL(cnt1)=nan;
                T6.g_std_LL(cnt1)=nan;
                T6.kNoN_LL(cnt1)=nan;
                T6.kNoN_std_LL(cnt1)=nan;
                T6.muN_LL(cnt1)=nan;
                T6.muN_std_LL(cnt1)=nan;

            end

        end

    end

    % Make sure cast and niskin are in text format in the table
    T2.cast=strcat(T2.cast,"'");
    T2.niskin=strcat(T2.niskin,"'");
    T3.cast=strcat(T3.cast,"'");
    T3.niskin=strcat(T3.niskin,"'");
    T4.cast=strcat(T4.cast,"'");
    T4.niskin=strcat(T4.niskin,"'");
    T5.cast=strcat(T5.cast,"'");
    T5.niskin=strcat(T5.niskin,"'");
    T6.cast=strcat(T6.cast,"'");
    T6.niskin=strcat(T6.niskin,"'");


    % Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    % cruise
    newname2=strrep(list(n1).name,'k-values','rates-GFF');%Replace raw by clean
    newtablename2=strcat(rep2,newname2);%New tablename and path
    writetable(T2,newtablename2)
    newname3=strrep(list(n1).name,'k-values','rates-u10');%Replace raw by clean
    newtablename3=strcat(rep2,newname3);%New tablename and path
    writetable(T3,newtablename3)
    newname4=strrep(list(n1).name,'k-values','rates-d10');%Replace raw by clean
    newtablename4=strcat(rep2,newname4);%New tablename and path
    writetable(T4,newtablename4)
    %only save for cruise with 10um size fractionation (>0&<10)
    if ~all(isnan(T1.k_dil_HL_10um_sf)) || ~all(isnan(T1.k_dil_LL_10um_sf))
        newname5=strrep(list(n1).name,'k-values','rates-10um-sf');%Replace raw by clean
        newtablename5=strcat(rep2,newname5);%New tablename and path
        writetable(T5,newtablename5)
    end
    %only save for cruise with no mesh experiment
    if ~all(isnan(T1.k_dil_HL_no_mesh)) || ~all(isnan(T1.k_dil_LL_no_mesh))
        newname6=strrep(list(n1).name,'k-values','rates-no-mesh');%Replace raw by clean
        newtablename6=strcat(rep2,newname6);%New tablename and path
        writetable(T6,newtablename6)
    end


end
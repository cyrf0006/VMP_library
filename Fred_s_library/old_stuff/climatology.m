% usage: climatology
% This script creates tprofiles.dat and sprofiles.dat corresponding of all
% profiles we have in a depth x no_profiles matrix. Dates corresponding to
% these profiles are saved in dateprofiles.dat
% This script also save T_climato_MM.dat and S_climato_MM.dat
% N2_climato_MM.dat corresponding to climatological profiles for each 
% months for this period (1993-2009).
% This script also save the standard deviation of these climatological
% profiles (T_climatoSTD_MM.dat, etc.)
%
% Please note that I may have commented the saving part of the script.
%
% Usually, you only need to run once this kind of script.
%
% author: F. Cyr, feb. 2010
%
% modifications: - add many output needed by CIL_stats
%                - compute climatology for model results (2010-09-01)
% ---------------------------------------------------------------------- %

clear

% some constants
g = 9.81;
rho_0 = 1.025e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
T_core = 1; %degC, threshold for the CIL core
freeze_pt = -1.8;

CIL_core = zeros(12, 8); % CIL_core = [ Tmin STD_Tmin Z(Tmin) STD_Zmin  CIL_thickness  STD_thickness Heat_cont. STD_heat]

%             ...   ...      ...     ...          ...             ...         ...       ...  
load nom.dat

s=size(nom);
no_files = s(1);

dz = 1;
P = 1:dz:300;

T(1:length(P),1:no_files) = NaN; % matrix to be saved
S(1:length(P),1:no_files) = NaN;

for i = 1:no_files
    
    
    yyyy = nom(i,1);
    mm = nom(i,2);
    dd = nom(i,3);
    hh = nom(i,4);
    mi  = nom(i,5);
    
    % transfert number into string for filename
    
    YYYY = num2str(yyyy);
    
    if mm < 10
        MM = ['0' num2str(mm)];
    else
        MM = num2str(mm);
    end
    
    if dd < 10
        DD = ['0' num2str(dd)];
    else
        DD = num2str(dd);
    end
        
    if hh < 10
        HH = ['0' num2str(hh)];
    else
        HH = num2str(hh);
    end

    if mi < 10
        MI = ['0' num2str(mi)];
    else
        MI = num2str(mi);
    end    
    
    fname = [YYYY '-' MM '-' DD '_' HH 'h' MI '.bin'];
    
    file = load(fname);
    
    if ~isempty(file)==1 % check if file empty
        
    dat(i) = datenum(yyyy, mm, dd, hh, mi, 0);
    
    pres = file(:,1);
    temp = file(:,2);
    sali = file(:,3);
    
    I = find(temp==-99);
    temp(I) = NaN;
    I = find(sali==-99);
    sali(I) = NaN;

    if max(pres) >= max(P) & min(pres)==min(P)
      
        I = find(file(:,1)>=min(P) & file(:,1)<=max(P));
        
        % fill matrix with empty value
        T(P,i)=NaN;
        S(P,i)=NaN;
        DENS(P,i)=NaN;
        N2(P,i)=NaN;
        
        % substituate empty value where it is possible
        T(pres(I),i) = temp(I);
        S(pres(I),i) = sali(I);
        
        DENS = sw_dens(sali(I), temp(I), pres(I));
        drho = diff(DENS);
        brunt = g/rho_0.*drho/dz;
        N2(pres(I(1:length(I)-1)),i) = brunt; % length(N2) == length(T) - 1
        
    else % CTD profile doesnt start at p=1m or pmax smaller than 300m
        
        I = find(file(:,1)>=max(min(pres), min(P)) & file(:,1)<=min(max(pres), max(P)));
        
        % fill matrix with empty value
        T(P,i)=NaN;
        S(P,i)=NaN;
        DENS(P,i)=NaN;
        N2(P,i)=NaN;       
        
        % substituate empty value where it is possible
        T(pres(I),i) = temp(I);
        S(pres(I),i) = sali(I);
   
        DENS = sw_dens(sali(I), temp(I), pres(I));
        drho = diff(DENS);
        brunt = g/rho_0.*drho/dz;
        N2(pres(I(1:length(I)-1)),i) = brunt; 
        
    end
  
    
       
    else %file empty
        continue
    end
    
end


% $$$ %Save all profiles in a matrix
% $$$ dlmwrite('tprofiles.dat', T,'delimiter',' ','precision',6);
% $$$ dlmwrite('sprofiles.dat', S,'delimiter',' ','precision',6);
% $$$ dlmwrite('datprofiles.dat', dat,'delimiter',' ','precision',6);

% compute CIL Heat content
[a, b, v] = find(T<T_core);

% computing CIL heat content and thickness
for i=1:no_files 
    A = find(b==i); %indices de a correspondants au profil i
    DENS = sw_dens(S(a(A),i), T(a(A),i), pres(a(A))); %density for bin into CIL core
    H(i) = cp*sum(DENS.*(T(a(A), i)-T_core))*dz; %Heat content of CIL for this profile (cf. G. Smith chap 4)
    L(i) = length(A); % CIL thickness following our definition 
end

% computing water column heat content between 0:50:300m
HH = zeros(3, no_files);
Heat_int = zeros(3, 12); %integrated heat content above freezpt
                         %(0-50, 50-200, 200-300m)  
Heat_int_std = zeros(3, 12);  

HH_T = HH; % Temp. class instead depth class
Heat_int_T = Heat_int;              
Heat_int_std_T = Heat_int_std;


for i=1:no_files 
    % 0-50 m
    I = find(P>0 & P<=50);
    %I = find(P>0 & P<=10);
    %I = find(P>0 & P<=40);
    %I = find(P>0 & P<=30);
    %I = find(P>0 & P<=20);
    %I = find(P>0 & P<=70);
    DENS = sw_dens(S(I,i), T(I,i), P(I)'); %density for bin into CIL core
    HH(1, i) = cp*nansum(DENS.*(T(I, i)-freeze_pt))*dz;
    
    % 50-200 m
    I = find(P>50 & P<=200);
    %I = find(P>10 & P<=200);
    %I = find(P>40 & P<=200);
    %I = find(P>30 & P<=200);
    %I = find(P>20 & P<=200);
    %I = find(P>60 & P<=200);
    %I = find(P>50 & P<=200);
    %I = find(P>70 & P<=200);
    
    DENS = sw_dens(S(I,i), T(I,i), P(I)'); %density for bin into CIL core
    HH(2, i) = cp*nansum(DENS.*(T(I, i)-freeze_pt))*dz;
    
    % 200-300m
    %I = find(P>200 & P<=300);
    I = find(P>200 & P<=300);
    DENS = sw_dens(S(I,i), T(I,i), P(I)'); %density for bin into CIL core

    HH(3, i) = cp*nansum(DENS.*(T(I, i)-freeze_pt))*dz;
    
    % CIL = T<1deg.C
    I = find(T(:,i)<=1);
    if ~isempty(I)==1
        A = min(I); %shallowest limit of CIL
        B = max(I); %deepest limit of CIL

        DENS = sw_dens(S(I,i), T(I,i), P(I)'); %density for bin into CIL core
 
        % - constant rho - %
% $$$         HH_T(2, i) = cp*nansum(rho_0.*(T(I, i)-freeze_pt))*dz./length(I);
        
        % - variable rho - %
        HH_T(2, i) = cp*nansum(DENS.*(T(I, i)-freeze_pt))*dz./length(I);
        
        % - not divided by d - %
% $$$         HH_T(2, i) = cp*nansum(DENS.*(T(I, i)-freeze_pt))*dz;

        DENS = sw_dens(S(A:B,i), T(A:B,i), P(A:B)'); %density for bin into CIL core        
        % above CIL
        if A>1
        DENS = sw_dens(S(1:A-1,i), T(1:A-1,i), P(1:A-1)'); %density for bin into CIL core
                                                           
        HH_T(1, i) = cp*nansum(DENS.*(T(1:A-1, i)-freeze_pt))*dz/length(1:A-1);
        %HH_T(1, i) = cp*nansum(DENS.*(T(1:A-1, i)-freeze_pt))*dz;
        else % CIL reaches surface, no surface layer
        HH_T(1, i)=0;
        end
        
        % under CIL
        DENS = sw_dens(S(B+1:300,i), T(B+1:300,i), P(B+1:300)'); %density for bin into CIL core
                                                                 
        HH_T(3, i) = cp*nansum(DENS.*(T(B+1:300, i)-freeze_pt))*dz/length(B+1:300);
        %HH_T(3, i) = cp*nansum(DENS.*(T(B+1:300, i)-freeze_pt))*dz;

    else
        HH_T(1, i) = NaN; 
        HH_T(2, i) = NaN; 
        HH_T(3, i) = NaN; 
    end
end

% $$$ % load diffus_model results
% $$$ T_model_kobs = load('T_diffus_daily_kobs.dat');
% $$$ T_model_kcst = load('T_diffus_daily_5e-5.dat');
% $$$ time_model = load('N_diffus_daily.dat');

% load diffus_model_noflux results
T_model_kobs = load('T_diffus_daily_kobs_TS.dat');
T_model_kcst = load('T_diffus_daily_5e-5_TS.dat');
S_model_kobs = load('S_diffus_daily_kobs_TS.dat');
S_model_kcst = load('S_diffus_daily_5e-5_TS.dat');
time_model = load('N_diffus_daily.dat');


% Build the climatology

for month = 4:11

    %T, S and N2
    I = find(str2num(datestr(dat,5))==month); %any profile this month
    T_climato = [P' nanmean(T(:,I), 2)];
    S_climato = [P' nanmean(S(:,I), 2)];    
    N2_climato = [P'+ 0.5 nanmean(N2(:,I), 2)];
    T_climato_STD = [P' nanstd(T(:,I), 2)];
    S_climato_STD = [P' nanstd(S(:,I), 2)];
    N2_climato_STD = [P'+ 0.5 nanstd(N2(:,I), 2)];  
    
    % same for model results
    II = find(str2num(datestr(time_model,5))==month);
    if ~isempty(II)==1
        Tobs_climato = [P' nanmean(T_model_kobs(:,II), 2)];
        Tcst_climato = [P' nanmean(T_model_kcst(:,II), 2)];
        Sobs_climato = [P' nanmean(S_model_kobs(:,II), 2)];
        Scst_climato = [P' nanmean(S_model_kcst(:,II), 2)];
    else
        Tobs_climato = [P' P'.*NaN];
        Tcst_climato = [P' P'.*NaN];
        Sobs_climato = [P' P'.*NaN];
        Scst_climato = [P' P'.*NaN];
    end
    

    %CIL core
    [Y, J] = min(T(:,I)); % position and value of Tmin
    CIL_core(month,1) = nanmean(Y);
    CIL_core(month,2) = nanstd(Y);
    CIL_core(month,3) = nanmean(P(J));
    CIL_core(month,4) = nanstd(P(J));
    CIL_core(month,5) = nanmean(L(I).*dz);
    CIL_core(month,6) = nanstd(L(I).*dz);
    CIL_core(month,7) = nanmean(H(I)); %heat content
    CIL_core(month,8) = nanstd(H(I));
    
    % Integrated heat content
    
    Heat_int(:, month) = nanmean(HH(:,I),2);
    Heat_int_std(:, month) = nanstd(HH(:,I),2);
    Heat_int_T(:, month) = nanmean(HH_T(:,I),2);
    Heat_int_std_T(:, month) = nanstd(HH_T(:,I),2);

    %Number of profile in each month (for statistics)
    No_prof(month) = length(I);
    
    
    if month <10
            Toutname  = sprintf('T_climato_0%d.dat', month);
            Soutname  = sprintf('S_climato_0%d.dat', month);
            N2outname = sprintf('N2_climato_0%d.dat', month);
            TSTDoutname = sprintf('T_climatoSTD_0%d.dat', month);
            SSTDoutname = sprintf('S_climatoSTD_0%d.dat', month);
            N2STDoutname = sprintf('N2_climatoSTD_0%d.dat', month);
            
            Tkobs_outname = sprintf('Tmodel_kobs_climato_0%d.dat', month);
            Tkcst_outname = sprintf('Tmodel_kcst_climato_0%d.dat', month);
            Skobs_outname = sprintf('Smodel_kobs_climato_0%d.dat', month);
            Skcst_outname = sprintf('Smodel_kcst_climato_0%d.dat', month);
    else
            Toutname  = sprintf('T_climato_%d.dat', month);
            Soutname  = sprintf('S_climato_%d.dat', month);
            N2outname  = sprintf('N2_climato_%d.dat', month);
            TSTDoutname = sprintf('T_climatoSTD_%d.dat', month);
            SSTDoutname = sprintf('S_climatoSTD_%d.dat', month);
            N2STDoutname = sprintf('N2_climatoSTD_%d.dat', month);
            
            Tkobs_outname = sprintf('Tmodel_kobs_climato_%d.dat', month);
            Tkcst_outname = sprintf('Tmodel_kcst_climato_%d.dat', month);
            Skobs_outname = sprintf('Smodel_kobs_climato_%d.dat', month);
            Skcst_outname = sprintf('Smodel_kcst_climato_%d.dat', month);
    end

    dlmwrite(Toutname, T_climato,'delimiter',' ','precision',6)
    dlmwrite(Soutname, S_climato,'delimiter',' ','precision',6)
    dlmwrite(N2outname, N2_climato,'delimiter',' ','precision',6)
    dlmwrite(TSTDoutname, T_climato_STD,'delimiter',' ','precision',6)    
    dlmwrite(SSTDoutname, S_climato_STD,'delimiter',' ','precision',6)      
    dlmwrite(N2STDoutname, N2_climato_STD,'delimiter',' ','precision',6) 
     
    dlmwrite(Tkobs_outname, Tobs_climato,'delimiter',' ','precision',6)      
    dlmwrite(Tkcst_outname, Tcst_climato,'delimiter',' ','precision',6) 
    dlmwrite(Skobs_outname, Sobs_climato,'delimiter',' ','precision',6)      
    dlmwrite(Skcst_outname, Scst_climato,'delimiter',' ','precision',6) 
end


    dlmwrite('CIL_core.dat', CIL_core,'delimiter',' ','precision',6) 
          
    dlmwrite('clim_no_profile.dat', No_prof,'delimiter',' ','precision',6)
    dlmwrite('HEAT_int_50-200.dat', Heat_int,'delimiter',' ','precision',6) 
    dlmwrite('HEAT_int_std_50-200.dat', Heat_int_std,'delimiter',' ','precision',6) 
    dlmwrite('HEAT_int_T.dat', Heat_int_T,'delimiter',' ','precision',6) 
    dlmwrite('HEAT_int_std_T.dat', Heat_int_std_T,'delimiter',' ','precision',6)

    
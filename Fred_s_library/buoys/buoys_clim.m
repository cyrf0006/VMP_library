% Buoys climatology
% This script has been written to build a climatology of wind/windstress
% for the IML/MPO buoys. Considering raw files from SGDO (ex: MMOB_BOUEE2002_RIMOUSKI_IML4_METOCE.ODF)
% the method is to 1st use a shell script to modify files from ODF to
% ASCII: "~/shellscripts/ODF2ASCII.sh MMOB_BOUEE2008_RIMOUSKI_IML4_METOCE.ODF meteo08 date08"
% and then just run this script after editing meteofiles and datefiles.
% 2 files will be created: wind_weekly.dat and wind_monthly.dat, containing
% [wind wind_std gust gust_std windstress windstress_std guststress guststress_std]
%  ...   ...    ...    ...        ...          ...           ...          ...     %
%
% author: F. Cyr, march 2010
%
% Modifications: - 
%
%
%
%
% ------------------------------------------------------------------------%

clear

meteofiles = ['meteo02';  'meteo03';  'meteo04';  'meteo05'; ...
              'meteo06';  'meteo07';  'meteo08'; 'meteo09'];
datefiles = ['date02';  'date03';  'date04';  'date05';  'date06'; ...
             'date07';  'date08'; 'date09'];

count = 1;

% store all data in very long vectors (all years in a row)
for i=1:size(meteofiles,1)
    
    system(['rm -rf OUT']);
    system(['sed -s "s/:/ /g" ' datefiles(i,:) ' > OUT']);
    system(['sed -s "s/-/ /g" OUT > /tmp/tmp.txt']);
    system(['mv /tmp/tmp.txt ./OUT']);
    
    meteo = load(meteofiles(i,:));
    datef = load('OUT');
    
    L = length(meteo(:,1)); %number of recorded step
    
    windsp(count:count+L-1) = meteo(:,3);
    windir(count:count+L-1) = meteo(:,7); %degrees
    gustsp(count:count+L-1) = meteo(:,5);
    pres(count:count+L-1) = meteo(:,13); %hPa
    temp(count:count+L-1) = meteo(:,9);
    sst(count:count+L-1) = meteo(:,15); 
    rh(count:count+L-1) = meteo(:,11)./100; % rel hum.
     
    yyyy(count:count+L-1) = datef(:,1);
    mm(count:count+L-1) = datef(:,2);
    dd(count:count+L-1) = datef(:,3);

    clear meteo datef

    count = count+L;
end

%vector with all data time
n = datenum(yyyy, mm, dd);


%compute wind_x wind_y
wind_x = windsp.*cosd(windir);
wind_y = windsp.*sind(windir);
gust_x = gustsp.*cosd(windir); 
gust_y = gustsp.*sind(windir);

% compute  windstress !!!

% -> Old version (from GIll)
% $$$ %compute drag coef (Gill)
% $$$ rho = pres.*100./287.1./(temp+273.15); %Kg/m^3
% $$$ tau_w = (0.61+0.063*windsp)/1000.*rho.*windsp.^2;
% $$$ tau_g = (0.61+0.063*gustsp)/1000.*rho.*gustsp.^2;
% $$$ tau_wx = (0.61+0.063*windsp)/1000.*rho.*windsp.*wind_x;
% $$$ tau_wy = (0.61+0.063*windsp)/1000.*rho.*windsp.*wind_y;
% $$$ tau_gx = (0.61+0.063*gustsp)/1000.*rho.*gustsp.*gust_x;
% $$$ tau_gy = (0.61+0.063*gustsp)/1000.*rho.*gustsp.*gust_y;
% $$$ 
% $$$ I = find(windsp<6);
% $$$ tau_w(I)=1.1e-3.*rho(I).*windsp(I).^2;
% $$$ tau_wx(I)=1.1e-3.*rho(I).*windsp(I).*wind_x(I);
% $$$ tau_wy(I)=1.1e-3.*rho(I).*windsp(I).*wind_y(I);
% $$$ I = find(gustsp<6);
% $$$ tau_g(I)=1.1e-3.*rho(I).*gustsp(I).^2;
% $$$ tau_gx(I)=1.1e-3.*rho(I).*gustsp(I).*gust_x(I);
% $$$ tau_gy(I)=1.1e-3.*rho(I).*gustsp(I).*gust_y(I);

% -> New version (Zedler et al. 2009)

%                        ----- FROM HERE -----
% calculation have been done for wind and gust, but bot for xy-dir wind
rho = pres.*100./287.1./(temp+273.15); %Kg/m^3

% parameters
kappa = 0.41;
cd_black = [1.46;1.7;1.8;1.6;1.5]; % for U10 = [15-20; 20-22; 22-25; 25-30; >30]
cd_pond = [1.2; 0.49; 0.065]; %[1.2; 0.49 + 0.065U10] for U10 = [0-11; 11-15];

% 1)  --- for average wind --- 

wind10 = windlog_convert(windsp, 2, 10); % from 2m to 10m
C_d = wind10.*0; % initialization
U = wind10;

% Drag coef calculation;
I = find(U<11);
C_d(I) = cd_pond(1);
I = find(U>=11 & U < 15);
C_d(I) = cd_pond(2) + cd_pond(3).*U(I);
I = find(U>=15 & U < 20);
C_d(I) = cd_black(1);
I = find(U>=20 & U < 22);
C_d(I) = cd_black(2);
I = find(U>=22 & U < 25);
C_d(I) = cd_black(3);
I = find(U>=25 & U < 30);
C_d(I) = cd_black(4);
I = find(U>=30);
C_d(I) = cd_black(5);

C_d = C_d./1000; % in 1e-3
tau_w = C_d.*rho.*U.^2;

clear U wind10 C_d

% 2)  --- for gusts --- 
gust10 = windlog_convert(gustsp, 2, 10); % from 2m to 10m
C_d = gust10.*0; % initialization
U = gust10;

% Drag coef calculation;
I = find(U<11);
C_d(I) = cd_pond(1);
I = find(U>=11 & U < 15);
C_d(I) = cd_pond(2) + cd_pond(3).*U(I);
I = find(U>=15 & U < 20);
C_d(I) = cd_black(1);
I = find(U>=20 & U < 22);
C_d(I) = cd_black(2);
I = find(U>=22 & U < 25);
C_d(I) = cd_black(3);
I = find(U>=25 & U < 30);
C_d(I) = cd_black(4);
I = find(U>=30);
C_d(I) = cd_black(5);

C_d = C_d./1000; % in 1e-3
tau_g = C_d.*rho.*U.^2;

clear U gust10 C_d


% 3)  --- for wind-x --- 
tau_wx = tau_w.*0;
% 4)  --- for wind-y ---
tau_wy = tau_w.*0;
% 5)  --- for gust-x --- 
tau_gx = tau_w.*0;
% 6)  --- for gust-y --- 
tau_gy = tau_w.*0;
% 3-6 not computed.

%                        ----- TO HERE -----

clear rho

%compute heat_flux (Gill p. 36)
% param (40-50deg.N)
dQdT=32; %W/m^2
dQdnc=135;
dQdr=-290;

nc=temp*0+0.3; %cloud cover

Q = (sst-temp).*dQdT + (nc-0.25).*dQdnc + (rh - 0.75).*dQdr;

% $$$ % save raw data
% $$$ n_raw = datenum(yyyy, mm, dd);
% $$$ dlmwrite('raw_wind_heat.dat', [n_raw' tau_wx' tau_wy' Q'],'delimiter',' ','precision',6);


%%%%%%%%%%%%%%%%%%
% - daily clim - %
%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty matrix for daily clim
daily_wind(1:365, 1:size(meteofiles,1))=NaN;
daily_T=daily_wind;
daily_SST = daily_wind;
daily_HR = daily_wind;
daily_pres = daily_wind;

years = 2002:2009;

for k = 1:length(years)

    
    %days in this year
    nn = datenum(years(k), 1, 1):datenum(years(k), 12, 31);
    
    %bissextile
    if length(nn)==366
        nn(60)=[];
    end
        
    
    for j = 1:length(nn)
 
        I = find(n == nn(j));

        daily_wind(j,k) = nanmean(windsp(I)); %windspeed, not tau
        daily_T(j,k) = nanmean(temp(I)); 
        daily_SST(j,k) = nanmean(sst(I));
        daily_HR(j,k) = nanmean(rh(I));
        daily_pres(j,k) = nanmean(pres(I));    
   
    end
    
    

end
clear temp sst nc rh pres

daily_wind_clim = nanmean(daily_wind, 2); %windspeed
daily_T_clim = nanmean(daily_T, 2);
daily_SST_clim = nanmean(daily_SST, 2);
daily_HR_clim = nanmean(daily_HR, 2);
daily_pres_clim = nanmean(daily_pres, 2);

daily_n_clim = [datenum(999, 1, 1):datenum(999, 12, 31)]';

% save daily clim
dlmwrite('daily_clim.dat', [daily_n_clim daily_wind_clim daily_T_clim daily_SST_clim ...
                  daily_HR_clim daily_pres_clim ],'delimiter',' ','precision',6);

% save raw daily values
dlmwrite('raw_daily_wind.dat',daily_wind ,'delimiter',' ','precision',6);
dlmwrite('raw_daily_T.dat',daily_T ,'delimiter',' ','precision',6);
dlmwrite('raw_daily_SST.dat',daily_SST ,'delimiter',' ','precision',6);
dlmwrite('raw_daily_HR.dat',daily_HR ,'delimiter',' ','precision',6);
dlmwrite('raw_daily_pres.dat',daily_pres ,'delimiter',' ','precision',6);


% $$$ for i = 1:12 % months to compute
% $$$  
% $$$     %%%%%%%%%%%%%%%%%%%
% $$$     % - weekly clim - %
% $$$     %%%%%%%%%%%%%%%%%%%
% $$$     
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     I = find(mm == i & dd<=7);
% $$$     W_w(i*4-3) = nanmean(windsp(I));
% $$$     W_wx(i*4-3) = nanmean(wind_x(I));
% $$$     W_wy(i*4-3) = nanmean(wind_y(I));
% $$$     G_w(i*4-3) = nanmean(gustsp(I));
% $$$     G_wx(i*4-3) = nanmean(gust_x(I));
% $$$     G_wy(i*4-3) = nanmean(gust_y(I));
% $$$     
% $$$     TAU_ww(i*4-3) = nanmean(tau_w(I));
% $$$     TAU_wwx(i*4-3) = nanmean(tau_wx(I));
% $$$     TAU_wwy(i*4-3) = nanmean(tau_wy(I));
% $$$     TAU_gw(i*4-3) = nanmean(tau_g(I));
% $$$     TAU_gwx(i*4-3) = nanmean(tau_gx(I));
% $$$     TAU_gwy(i*4-3) = nanmean(tau_gy(I));
% $$$     Q_w(i*4-3) = nanmean(Q(I));
% $$$     
% $$$     STD_W_w(i*4-3) = nanstd(windsp(I));
% $$$     STD_W_wx(i*4-3) = nanstd(wind_x(I));
% $$$     STD_W_wy(i*4-3) = nanstd(wind_y(I));
% $$$     STD_G_w(i*4-3) = nanstd(gustsp(I));
% $$$     STD_G_wx(i*4-3) = nanstd(gust_x(I));
% $$$     STD_G_wy(i*4-3) = nanstd(gust_y(I));
% $$$     STD_TAU_ww(i*4-3) = nanstd(tau_w(I));
% $$$     STD_TAU_wwx(i*4-3) = nanstd(tau_wx(I));
% $$$     STD_TAU_wwy(i*4-3) = nanstd(tau_wy(I));
% $$$     STD_TAU_gw(i*4-3) = nanstd(tau_g(I));
% $$$     STD_TAU_gwx(i*4-3) = nanstd(tau_gx(I));
% $$$     STD_TAU_gwy(i*4-3) = nanstd(tau_gy(I));
% $$$     STD_Q_w(i*4-3) = nanstd(Q(I));
% $$$     n_w(i*4-3) = datenum(999, i, 4, 0,0,0);
% $$$     
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     I = find(mm == i & dd<=14 & dd>7);
% $$$     W_w(i*4-2) = nanmean(windsp(I));
% $$$     W_wx(i*4-2) = nanmean(wind_x(I));
% $$$     W_wy(i*4-2) = nanmean(wind_y(I));
% $$$     G_w(i*4-2) = nanmean(gustsp(I));
% $$$     G_wx(i*4-2) = nanmean(gust_x(I));
% $$$     G_wy(i*4-2) = nanmean(gust_y(I));
% $$$     
% $$$     TAU_ww(i*4-2) = nanmean(tau_w(I));
% $$$     TAU_wwx(i*4-2) = nanmean(tau_wx(I));
% $$$     TAU_wwy(i*4-2) = nanmean(tau_wy(I));
% $$$     TAU_gw(i*4-2) = nanmean(tau_g(I));
% $$$     TAU_gwx(i*4-2) = nanmean(tau_gx(I));
% $$$     TAU_gwy(i*4-2) = nanmean(tau_gy(I));
% $$$     Q_w(i*4-2) = nanmean(Q(I));
% $$$     
% $$$     STD_W_w(i*4-2) = nanstd(windsp(I));
% $$$     STD_W_wx(i*4-2) = nanstd(wind_x(I));
% $$$     STD_W_wy(i*4-2) = nanstd(wind_y(I));
% $$$     STD_G_w(i*4-2) = nanstd(gustsp(I));
% $$$     STD_G_wx(i*4-2) = nanstd(gust_x(I));
% $$$     STD_G_wy(i*4-2) = nanstd(gust_y(I));
% $$$     STD_TAU_ww(i*4-2) = nanstd(tau_w(I));
% $$$     STD_TAU_wwx(i*4-2) = nanstd(tau_wx(I));
% $$$     STD_TAU_wwy(i*4-2) = nanstd(tau_wy(I));
% $$$     STD_TAU_gw(i*4-2) = nanstd(tau_g(I)); 
% $$$     STD_TAU_gwx(i*4-2) = nanstd(tau_gx(I));
% $$$     STD_TAU_gwy(i*4-2) = nanstd(tau_gy(I));
% $$$     STD_Q_w(i*4-2) = nanstd(Q(I));     
% $$$     n_w(i*4-2) = datenum(999, i, 11, 0,0,0);
% $$$     
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     I = find(mm == i & dd<=21 & dd>14);
% $$$     W_w(i*4-1) = nanmean(windsp(I));
% $$$     W_wx(i*4-1) = nanmean(wind_x(I));
% $$$     W_wy(i*4-1) = nanmean(wind_y(I));
% $$$     G_w(i*4-1) = nanmean(gustsp(I));
% $$$     G_wx(i*4-1) = nanmean(gust_x(I));
% $$$     G_wy(i*4-1) = nanmean(gust_y(I));
% $$$     
% $$$     TAU_ww(i*4-1) = nanmean(tau_w(I));
% $$$     TAU_wwx(i*4-1) = nanmean(tau_wx(I));
% $$$     TAU_wwy(i*4-1) = nanmean(tau_wy(I));
% $$$     TAU_gw(i*4-1) = nanmean(tau_g(I));
% $$$     TAU_gwx(i*4-1) = nanmean(tau_gx(I));
% $$$     TAU_gwy(i*4-1) = nanmean(tau_gy(I));
% $$$     Q_w(i*4-1) = nanmean(Q(I));
% $$$     
% $$$     STD_W_w(i*4-1) = nanstd(windsp(I));
% $$$     STD_W_wx(i*4-1) = nanstd(wind_x(I));
% $$$     STD_W_wy(i*4-1) = nanstd(wind_y(I));
% $$$     STD_G_w(i*4-1) = nanstd(gustsp(I));
% $$$     STD_G_wx(i*4-1) = nanstd(gust_x(I));
% $$$     STD_G_wy(i*4-1) = nanstd(gust_y(I));
% $$$     STD_TAU_ww(i*4-1) = nanstd(tau_w(I));
% $$$     STD_TAU_wwx(i*4-1) = nanstd(tau_wx(I));
% $$$     STD_TAU_wwy(i*4-1) = nanstd(tau_wy(I));
% $$$     STD_TAU_gw(i*4-1) = nanstd(tau_g(I));
% $$$     STD_TAU_gwx(i*4-1) = nanstd(tau_gx(I));
% $$$     STD_TAU_gwy(i*4-1) = nanstd(tau_gy(I));
% $$$     STD_Q_w(i*4-1) = nanstd(Q(I));
% $$$     n_w(i*4-1) = datenum(999, i, 18, 0,0,0);
% $$$     
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     I = find(mm == i & dd>21);
% $$$     W_w(i*4) = nanmean(windsp(I));
% $$$     W_wx(i*4) = nanmean(wind_x(I));
% $$$     W_wy(i*4) = nanmean(wind_y(I));
% $$$     G_w(i*4) = nanmean(gustsp(I));
% $$$     G_wx(i*4) = nanmean(gust_x(I));
% $$$     G_wy(i*4) = nanmean(gust_y(I));
% $$$     
% $$$     TAU_ww(i*4) = nanmean(tau_w(I));
% $$$     TAU_wwx(i*4) = nanmean(tau_wx(I));
% $$$     TAU_wwy(i*4) = nanmean(tau_wy(I));
% $$$     TAU_gw(i*4) = nanmean(tau_g(I));
% $$$     TAU_gwx(i*4) = nanmean(tau_gx(I));
% $$$     TAU_gwy(i*4) = nanmean(tau_gy(I));
% $$$     Q_w(i*4) = nanmean(Q(I));
% $$$     
% $$$     STD_W_w(i*4) = nanstd(windsp(I));
% $$$     STD_W_wx(i*4) = nanstd(wind_x(I));
% $$$     STD_W_wy(i*4) = nanstd(wind_y(I));
% $$$     STD_G_w(i*4) = nanstd(gustsp(I));
% $$$     STD_G_wx(i*4) = nanstd(gust_x(I));
% $$$     STD_G_wy(i*4) = nanstd(gust_y(I));
% $$$     STD_TAU_ww(i*4) = nanstd(tau_w(I));
% $$$     STD_TAU_wwx(i*4) = nanstd(tau_wx(I));
% $$$     STD_TAU_wwy(i*4) = nanstd(tau_wy(I));
% $$$     STD_TAU_gw(i*4) = nanstd(tau_g(I));
% $$$     STD_TAU_gwx(i*4) = nanstd(tau_gx(I));
% $$$     STD_TAU_gwy(i*4) = nanstd(tau_gy(I));
% $$$     STD_Q_w(i*4) = nanstd(Q(I));  
% $$$     n_w(i*4) = datenum(999, i, 25, 0,0,0);
% $$$     
% $$$     %%%%%%%%%%%%%%%%%%%%
% $$$     % - monthly clim - %
% $$$     %%%%%%%%%%%%%%%%%%%%
% $$$     I = find(mm == i);
% $$$     W_m(i) = nanmean(windsp(I));
% $$$     W_mx(i) = nanmean(wind_x(I));
% $$$     W_my(i) = nanmean(wind_y(I));
% $$$     G_m(i) = nanmean(gustsp(I));
% $$$     G_mx(i) = nanmean(gust_x(I));
% $$$     G_my(i) = nanmean(gust_y(I));
% $$$     
% $$$     TAU_wm(i) = nanmean(tau_w(I));
% $$$     TAU_wmx(i) = nanmean(tau_wx(I));
% $$$     TAU_wmy(i) = nanmean(tau_wy(I));
% $$$     TAU_gm(i) = nanmean(tau_g(I));
% $$$     TAU_gmx(i) = nanmean(tau_gx(I));
% $$$     TAU_gmy(i) = nanmean(tau_gy(I));
% $$$     Q_m(i) = nanmean(Q(I));
% $$$     
% $$$     STD_W_m(i) = nanstd(windsp(I));
% $$$     STD_W_mx(i) = nanstd(wind_x(I));
% $$$     STD_W_my(i) = nanstd(wind_y(I));
% $$$     STD_G_m(i) = nanstd(gustsp(I));
% $$$     STD_G_mx(i) = nanstd(gust_x(I));
% $$$     STD_G_my(i) = nanstd(gust_y(I));
% $$$     STD_TAU_wm(i) = nanstd(tau_w(I));
% $$$     STD_TAU_wmx(i) = nanstd(tau_wx(I));
% $$$     STD_TAU_wmy(i) = nanstd(tau_wy(I));
% $$$     STD_TAU_gm(i) = nanstd(tau_g(I));
% $$$     STD_TAU_gmx(i) = nanstd(tau_gx(I));
% $$$     STD_TAU_gmy(i) = nanstd(tau_gy(I));
% $$$     STD_Q_m(i) = nanstd(Q(I));
% $$$     n_m(i) = datenum(999, i, 15, 0,0,0);
% $$$     
% $$$ end
% $$$ 
% $$$ % save climatology [wind wind_std gust gust_std windstress windstress_std guststress guststress_std]
% $$$ %                    ...   ...    ...    ...        ...          ...
% $$$ %                    ...          ...     %
% $$$ 
% $$$ dlmwrite('wind_weekly.dat', [W_w' STD_W_w' G_w' STD_G_w' TAU_ww' ...
% $$$                     STD_TAU_ww' TAU_gw' STD_TAU_gw'],'delimiter',' ','precision',6);
% $$$ dlmwrite('windx_weekly.dat', [W_wx' STD_W_wx' G_wx' STD_G_wx' TAU_wwx' ...
% $$$                     STD_TAU_wwx' TAU_gwx' STD_TAU_gwx'],'delimiter',' ','precision',6);
% $$$ dlmwrite('windy_weekly.dat', [W_wy' STD_W_wy' G_wy' STD_G_wy' TAU_wwy' ...
% $$$                     STD_TAU_wwy' TAU_gwy' STD_TAU_gwy'],'delimiter',' ','precision',6);
% $$$ 
% $$$ dlmwrite('wind_monthly.dat', [W_m' STD_W_m' G_m' STD_G_m' TAU_wm' ...
% $$$                     STD_TAU_wm' TAU_gm' STD_TAU_gm'],'delimiter',' ','precision',6);
% $$$ dlmwrite('windx_monthly.dat', [W_mx' STD_W_mx' G_mx' STD_G_mx' TAU_wmx' ...
% $$$                     STD_TAU_wmx' TAU_gmx' STD_TAU_gmx'],'delimiter',' ','precision',6);
% $$$ dlmwrite('windy_monthly.dat', [W_my' STD_W_my' G_my' STD_G_my' TAU_wmy' ...
% $$$                     STD_TAU_wmy' TAU_gmy' STD_TAU_gmy'],'delimiter',' ','precision',6);
% $$$ 
% $$$ dlmwrite('heat_weekly.dat', [W_w' STD_Q_w'],'delimiter',' ','precision',6);
% $$$ dlmwrite('heat_monthly.dat', [Q_m' STD_Q_m'],'delimiter',' ','precision',6);
% $$$ 

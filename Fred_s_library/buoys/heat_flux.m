% Net heat flux budget computed from Wallcraft et al. 2008:
% Qnet = Qsw+Qlw+Ql+Qs. Because of missing data, Long-wave 
% and short wave are computed from NARR reanalyses. Note
% that the old method to compute LW was not that bad, but the one
% for SW was terrible...
% Net heat budget is saved as a daily climatology (365 days), but
% the winter part of the year is filled with NaN values.
% Climatology is saved under the name 'daily_heatflx_clim.dat'
% and has the format [n Qnet]
%
% Author: F. Cyr (may 2010)
%
%
% Modifications: - may 4 : compute for all years before averaging
%                  (second part of the script)
%                  ** I kept the same name to save the result **  
%
% ---------------------------------------------------- %

% $$$ 
% $$$ clear
% $$$ % $$$ 
% $$$ % $$$ dat = load('datefile_heatflux.dat');
% $$$ % $$$ n = datenum(dat(:,1), dat(:,2), dat(:,3), dat(:,4), dat(:,5), dat(:,6));
% $$$ % $$$ clear dat
% $$$ % $$$ 
% $$$ % $$$ meteo = load('IML4_meteo.dat');
% $$$ % $$$ wind = meteo(:,3);
% $$$ % $$$ temp = meteo(:,9);
% $$$ % $$$ temp_K = temp + 273.15; %temp in Kelvins
% $$$ % $$$ SST = meteo(:,15);
% $$$ % $$$ SST_K = SST + 273.15;
% $$$ % $$$ HR = meteo(:,11); % in %
% $$$ % $$$ pres = meteo(:,13); %in mB
% $$$ % $$$ pres_kPa = pres./10; %in kPa
% $$$ % $$$ %pres_HG = pres.*760./1013.25; %pressure in mmHG
% $$$ % $$$ clear meteo
% $$$ % $$$ 
% $$$ % $$$ optic = load('IML4_optic.dat');
% $$$ % $$$ radiance = [ optic(:,3) optic(:,7) optic(:,11) optic(:,15) optic(:,19) optic(:,23) optic(:,27)]; %in uW/cm^2/nm/sr
% $$$ % $$$ clear optic
% $$$ 
% $$$ %daily_clim: [n wind T SST HR pres]
% $$$ daily_clim = load('daily_clim.dat');
% $$$ 
% $$$ n = daily_clim(:,1);
% $$$ wind = daily_clim(:,2);
% $$$ temp = daily_clim(:,3);
% $$$ temp_K = temp + 273.15; %temp in Kelvins
% $$$ SST = daily_clim(:,4);
% $$$ SST_K = SST + 273.15;
% $$$ HR = daily_clim(:,5);
% $$$ pres = daily_clim(:,6);
% $$$ pres_kPa = pres./10; %in kPa
% $$$ %pres_HG = pres.*760./1013.25; %pressure in mmHG
% $$$ 
% $$$ %%%%%%%%%%%%%%
% $$$ %% Constants %
% $$$ %%%%%%%%%%%%%%
% $$$ 
% $$$ cc(1:length(temp)) = 0.3; %cloud_cover N/A, so... constant!
% $$$ cc = cc';
% $$$ 
% $$$ 
% $$$ Cp = 1004.5; %J/kg/K (air spec. heat)
% $$$ R = 287.1; %J/kg/K (ideal gas cst)
% $$$ L = 2.5e6; %J/Kg (Latent heat of vapo) 
% $$$ 
% $$$ delta = 0.95; % sea water emissivity
% $$$ sigma = 5.6698e-8; %stephen-boltzmann constant
% $$$ c1 = 0.70; % cloud coeff (smith-dobson 84)
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%
% $$$ % diagnostic vars %
% $$$ %%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ ees = 610.78*exp((17.2694.*temp)./(temp+238.3)); %saturation vapor pressure (in Pa)
% $$$ ee = HR.*ees; %vapor pressure
% $$$ qa = 0.622.*ee./(pres_kPa.*1000-ee); %air mixing ratio [adim]
% $$$ qs = 0.622.*ees./(pres_kPa.*1000-ees); %air sat.  mixing ratio [adim]
% $$$ 
% $$$ Va = max(3, min(27.5, wind)); %new wind computed for latent heat (from Kara et al. 2000)
% $$$ 
% $$$ %%%%%%%%%%%%%
% $$$ %% LongWave %
% $$$ %%%%%%%%%%%%%
% $$$ 
% $$$ % $$$ % computed from wallcraft et al. 2008
% $$$ % $$$ Qlw = -delta.*sigma.*temp_K.^3.*(temp_K.*(0.254-0.0066.*ee*760/ ...
% $$$ 
% $$$ % loaded from NARR project
% $$$ LW = load('clim_netlw.dat');
% $$$ Qlw = LW(:,2);
% $$$ 
% $$$ figure
% $$$ plot(n, Qlw)
% $$$ title('longwave irradiance')
% $$$ datetick('x', 12)
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%%%
% $$$ %% Latent heat %
% $$$ %%%%%%%%%%%%%%%%
% $$$ 
% $$$ % computed from wallcraft et al. 2008
% $$$ Cl0 = 1e-3.*(0.994 + 0.061.*Va - 0.001.*Va.^2);
% $$$ Cl1 = 1e-3.*(-0.020 + 0.691./Va - 0.817./Va.^2);
% $$$ Cl = Cl0 + Cl1.*(SST-temp);
% $$$ 
% $$$ rho_a = pres_kPa*1000./R./temp_K; %in kg/m^3, pressure in Pa
% $$$ 
% $$$ Ql = Cl.*L.*rho_a.*Va.*(qa-qs); %[Ql] = W/m^2 id [q]=[]
% $$$ 
% $$$ figure
% $$$ plot(n, Ql)
% $$$ title('latent heat')
% $$$ datetick('x', 12)
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%
% $$$ % Sensible heat %
% $$$ %%%%%%%%%%%%%%%%%
% $$$ 
% $$$ % computed from wallcraft et al. 2008
% $$$ Cs = 0.96.*Cl;
% $$$ 
% $$$ Qs = Cs.*Cp.*rho_a.*Va.*(temp_K-SST_K);
% $$$ 
% $$$ figure
% $$$ plot(n, Qs)
% $$$ title('sensible heat')
% $$$ datetick('x', 12)
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%
% $$$ %% ShortWave %
% $$$ %%%%%%%%%%%%%%
% $$$ 
% $$$ % $$$ % fisrt method
% $$$ % $$$ spec = [412;443;490;510;555;669;683]; %mesured wavelength
% $$$ % $$$ theo = 400:1:700; %wavelength of the visible spectrum
% $$$ % $$$ theo = theo';
% $$$ % $$$ % integration in the visible spectrum (400-700nm), approx with
% $$$ % $$$ % ---412---------443---------490-------510-------555-----------------669--------683--%
% $$$ % $$$ %interpolation over the sepctrum
% $$$ % $$$ YI = INTERP1(spec,radiance',theo,'nearest'); % still in uW/cm^2/nm/sr
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ % %second method
% $$$ % $$$ % spec = [300;412;443;490;510;555;669;683;1000]; %mesured wavelength
% $$$ % $$$ % theo = 300:1000; %wavelength of the visible spectrum
% $$$ % $$$ % theo = theo';
% $$$ % $$$ % 
% $$$ % $$$ % S = size(radiance);
% $$$ % $$$ % Zer(1:S(1)) = 0;
% $$$ % $$$ % rad = [Zer' radiance Zer'];
% $$$ % $$$ % 
% $$$ % $$$ % %interpolation over the sepctrum
% $$$ % $$$ % YI = INTERP1(spec,rad',theo,'cubic'); % still in uW/cm^2/nm/sr
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ clear radiance
% $$$ % $$$ 
% $$$ % $$$ %integration
% $$$ % $$$ Qsw = nansum(YI, 1); % now in uW/cm^2/sr
% $$$ % $$$ %clear YI
% $$$ % $$$ 
% $$$ % $$$ %radiance2irradiance
% $$$ % $$$ Qsw = pi.*Qsw'; % now in uW/cm^2
% $$$ % $$$ 
% $$$ % $$$ % radiance in W/m^2
% $$$ % $$$ Qsw = Qsw./100;
% $$$ 
% $$$ % loaded from NARR project
% $$$ SW = load('clim_netsw.dat');
% $$$ 
% $$$ Qsw = SW(:,2);
% $$$ 
% $$$ figure
% $$$ plot(n, Qsw)
% $$$ title('shortwave irradiance')
% $$$ datetick('x', 12)
% $$$ 
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%
% $$$ % Net heat flux %
% $$$ %%%%%%%%%%%%%%%%%
% $$$ Q = Qsw + Qlw + Ql + Qs;
% $$$ figure
% $$$ plot(n, Q)
% $$$ title('Net Heat Flux')
% $$$ datetick('x', 12)
% $$$ 
% $$$ dlmwrite('daily_heatflx_clim.dat', [n Q],'delimiter',' ','precision',6);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- Second part, compute for all years before averaging --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

wind = load('raw_daily_wind.dat');
temp = load('raw_daily_T.dat');
temp_K = temp + 273.15; %temp in Kelvins
SST = load('raw_daily_SST.dat'); 
SST_K = SST + 273.15;
HR = load('raw_daily_HR.dat'); 
pres = load('raw_daily_pres.dat'); 
pres_kPa = pres./10; %in kPa
%pres_HG = pres.*760./1013.25; %pressure in mmHG


%%%%%%%%%%%%%%
%% Constants %
%%%%%%%%%%%%%%

cc(1:length(temp)) = 0.3; %cloud_cover N/A, so... constant!
cc = cc';


Cp = 1004.5; %J/kg/K (air spec. heat)
R = 287.1; %J/kg/K (ideal gas cst)
L = 2.5e6; %J/Kg (Latent heat of vapo) 

delta = 0.95; % sea water emissivity
sigma = 5.6698e-8; %stephen-boltzmann constant
c1 = 0.70; % cloud coeff (smith-dobson 84)

%%%%%%%%%%%%%%%%%%%
% diagnostic vars %
%%%%%%%%%%%%%%%%%%%

ees = 610.78*exp((17.2694.*temp)./(temp+238.3)); %saturation vapor pressure (in Pa)
ee = HR.*ees; %vapor pressure
qa = 0.622.*ee./(pres_kPa.*1000-ee); %air mixing ratio [adim]
qs = 0.622.*ees./(pres_kPa.*1000-ees); %air sat.  mixing ratio [adim]

Va = max(3, min(27.5, wind)); %new wind computed for latent heat (from Kara et al. 2000)

%%%%%%%%%%%%%
%% LongWave %
%%%%%%%%%%%%%

% loaded from NARR project
Qlw = load('raw_netlw.dat');

%%%%%%%%%%%%%%%%
%% Latent heat %
%%%%%%%%%%%%%%%%

% computed from wallcraft et al. 2008
Cl0 = 1e-3.*(0.994 + 0.061.*Va - 0.001.*Va.^2);
Cl1 = 1e-3.*(-0.020 + 0.691./Va - 0.817./Va.^2);
Cl = Cl0 + Cl1.*(SST-temp);

rho_a = pres_kPa*1000./R./temp_K; %in kg/m^3, pressure in Pa

Ql = Cl.*L.*rho_a.*Va.*(qa-qs); %[Ql] = W/m^2 id [q]=[]



%%%%%%%%%%%%%%%%%
% Sensible heat %
%%%%%%%%%%%%%%%%%

% computed from wallcraft et al. 2008
Cs = 0.96.*Cl;

Qs = Cs.*Cp.*rho_a.*Va.*(temp_K-SST_K);


%%%%%%%%%%%%%%
%% ShortWave %
%%%%%%%%%%%%%%

% loaded from NARR project
Qsw = load('raw_netsw.dat');


%%%%%%%%%%%%%%%%%
% Net heat flux %
%%%%%%%%%%%%%%%%%
Q = Qsw + Qlw + Ql + Qs;
keyboard
daily_clim = load('daily_clim.dat');
n = daily_clim(:,1);

dlmwrite('daily_heatflx_clim.dat', [n nanmean(Q,2) nanstd(Q,2)],'delimiter',' ','precision',6);



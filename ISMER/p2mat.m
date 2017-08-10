function p2mat(data_fname);

% This script converts raw VMP data (.p files) to physical measurements
% and stored them in a Matlab file (.mat).
%
% The only parameter is the name of the .P file without the extension (ex. 'DAT000')
% Note that the setup file muste be 'SETUP.TXT'
% The function creates 'data_fname.mat' and clear other variables
%
% Naming convention:
%   Upper-case variables mean fine-scale measurements (e.g. Seabird, Pressure)
%   Lower-case variables mean micro-scale measurements (e.g. shear)
%
% Author: Daniel Bourgault - 2009/07/24
%   - modified by F. Cyr - 2009/04/28: - now a function
%                                   - returns also conductivity
%   - modified by F. Cyr - 2009/12/15: - remove the "despike" to keep raw
%  measurements, we were doing it again in further treatments
%   - modified by F. Cyr - 2010/01/?:
%      use now fix_bad_buffers in this function 
%        (!!!! ERASE ORIGINAL .P FILES !!!!)
%   - modified by F. Cyr - 2010/03/23:
%      add a condition to deconvolve variables to make sure
%      variable exist. For instance, Galbraith has no
%      micro-conductivity and this was causing problem in this
%      script. By the way I also needed to modify the variables are
%      saved that;s why I always clear unused variables
%   - modified by F. Cyr - Dec. 2011:
%       - Correct conductivity for thermal lag using (garau_etal2011)
%       - Compute practical and absolute salinity and conservative
%       temperature from TEOS-10. New variables (SP, SA, CT) have been created
%       such that old variable (SBS, SBT, SBC) can still be used,
%       and they have been corrected for thermal lag only. 
%       ** NOTE: gsw_SA_from_SP.m needs LAT,LON coordinates of CTD
%       profile. These variable are passed to the function via
%       CALIB.DAT. If they are absent, SA, SP, CT won't be
%       calculated. If lat-lon is there, they will be saved with
%       the profile to be used by var_profile_cal.m
%       - this finction now returns a struct with variables to
%       var_profile_cal.m, so this function doesnt need to load the
%       .mat file each time it works with a profile. (to speed up)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
%FC data_fname = 'DAT000';    % Without the .P extension 
setup_fname = 'SETUP.TXT';

%% Calibration parameters 
% Read calibration parameters from file 'calib.dat'
% This file should look like this:
%
% G1 = 1.045;  % Shear Channel 1 differentiator gain.
% G2 = 1.010;  % Shear Channel 2 differentiator gain.
% S1 = 0.1240; % Shear probe sensitivity (M484)
% S2 = 0.1269; % Shear probe sensitivity (M485)
%
% % Diff. gain for pressure deconvolution (from calibration report)
% diff_gain_P  = 20.1;
% diff_gain_t1 = 0.980;
% diff_gain_t2 = 0.995;
% diff_gain_c1 = 1.000;
%
% % Calibration coefficients for fluorescence and transmissometer
% A_trans = -3.195060;
% B_trans = 5.508724e-2;
% A_fluoro = -6.41017;
% B_fluoro = 8.90302e-2;

% % Lat Lon of the cast (needed for GSW tool kit)
% profile_lat = 48.6667;
% profile_lon = -69.5;


calib_fname = 'CALIB.DAT';
fid = fopen(calib_fname);
tline = fgetl(fid);
while ischar(tline)
  eval(tline);
  tline = fgetl(fid);
end
fclose(fid);

% Correct file for "bad_buffer" problem...  F. Cyr modif.
% WILL ERASE THE ORIGINAL .P FILE!!!
[bad_records, chop, truncate] = fix_bad_buffers([data_fname, '.P']);
if length(bad_records)>0
    disp('fix_bad_buffer.m had been applied!');
    disp(sprintf('number of points chopped at the begining of the file =%d', chop));
    disp(sprintf('number of points truncated at the end of the file =%d', truncate));
    disp(sprintf('number of points fixed by the method =%d', bad_records));
    disp('Press any key to continue');
    pause
end
clear bad_records chop truncate


% Read the .p file and write the data to a mat-file.
% From this step the data are not yet in physical units. 
read_odas([data_fname]);

% Load the data from the .mat file just created above.
load([data_fname,'.mat']);

% Determine from the header the number of fast and slow column
% and other information from the address matrix.
matrix_rows  = header(2,31);
fast_columns = header(2,29);
slow_columns = header(2,30);
columns = slow_columns + fast_columns;

clock_rate = header(2,21) + header(2,22)/1000; 
fs = clock_rate / columns; % Calculate the exact fast sampling frequency.
FS = fs / matrix_rows;     % Calculate the exact slow sampling frequency.
fs_units = 'Hz'; FS_units = 'Hz';
clear clock_rate column slow_columns fast_columns matrix_rows

% Make time vectors
nmax_slow = length(P); 
nmax_fast = length(sh1);
t_slow = (0:nmax_slow-1)'/FS; 
t_fast = (0:nmax_fast-1)'/fs;
clear nmax_slow nmax_fast 

% Pressure deconvolution and conversion to physical units.
Ptmp = deconvolve('pressure', P, P_dP, FS, diff_gain_P);
clear P P_dP
[P, P_units] = convert_odas(Ptmp,'pres','file',setup_fname);
clear Ptmp;

% Interpolate the pressure P (slow sampling) to the fast sampling rate.
p =  interp1(t_slow, P, t_fast, 'linear','extrap'); % There<s a lag bet. P and p (see correction further)
p_units = 'dBar';

% Make mtime vectors
hh =        str2num(time(:,1:2)); 
mm =        str2num(time(:,4:5)); 
ss =        str2num(time(:,7:12)); 
the_day =   str2num(date(:,1:2));
the_month = str2num(date(:,4:5));
the_year =  str2num(date(:,7:10));

MTIME_tmp = datenum(the_year,the_month,the_day,hh,mm,ss);

IMTIME_tmp = linspace(MTIME_tmp(1),MTIME_tmp(end),length(MTIME_tmp));
IMTIME     = linspace(MTIME_tmp(1),MTIME_tmp(end),length(MTIME_tmp)*round(FS));
iMTIME     = linspace(MTIME_tmp(1),MTIME_tmp(end),length(MTIME_tmp)*round(fs));
MTIME = interp1(IMTIME_tmp,MTIME_tmp,IMTIME','linear','extrap');
mtime = interp1(IMTIME_tmp,MTIME_tmp,iMTIME','linear','extrap');
clear MTIME_tmp IMTIME_tmp IMTIME iMTIME;
clear hh mm ss the_day the_month the_year;

% Calculate fall speed W from P and interpolate at the fast sampling rate. 
% W is smoothed at 1.5 Hz.
Nyquist = FS/2; 
Wn = 1.5/Nyquist;
W = gradient(P, 1/FS);
[b,a] = butter(4,Wn);
W = filtfilt(b,a,W);
w = interp1(t_slow, W, t_fast,'linear','extrap');
W_units = 'm/s'; w_units = 'm/s';
clear a b Nyquist Wn

% Convert accelerometer data to physical units.
[ax, ax_units] = convert_odas(Ax,'Pitch','file',setup_fname);
[ay, ay_units] = convert_odas(Ay,'Roll','file',setup_fname);
[az, az_units] = convert_odas(Az,'Az','file',setup_fname);
g = 9.81;            % Gravitational acceleration.
pitch = asind(ax/g); % Convert ax in m/s^2 to degree.
roll = asind(ay/g);  % Convert ay in m/s^2 to degree.
pitch_units = 'degree'; 
roll_units = 'degree';
clear Ax Ay Az g ax ay

% Correction of the pressure shift between p and P
p = p + dp; %corrected pressure
clear dp;

% Convert shear probe signals to physical units (s^{-1}).
if exist('sh1')~=0
    shear1 = byte2shear(sh1,S1,G1,w);
    shear1_units = 's^{-1}'; 
    clear sh1
end
if exist('sh2')~=0
    shear2 = byte2shear(sh2,S2,G2,w);
    shear2_units = 's^{-1}';
    clear sh2
end

% Convert SeaBird T and C signals to physical units and C to salinity
if exist('sbt')~=0 & exist('sbc')~=0 
    [SBT, SBT_units] = convert_odas(sbt,'SBT1E','file',setup_fname);
    [SBC, SBC_units] = convert_odas(sbc,'SBC1E','file',setup_fname);

    % correct Conductivity for thermal Lag
    basicProfileData.ptime = MTIME;
    basicProfileData.depth = P;
    basicProfileData.temp = SBT;
    basicProfileData.cond = SBC;
    prof_corr = correctThermalLag(basicProfileData);    
    SBC = prof_corr.cond;
    SBT = prof_corr.temp;
    
    SBS = gsw_SP_from_C(SBC, SBT, P); %SP
    SBS_units = 'psu'; 

    if exist('profile_lat') ~= 0 &  exist('profile_lon') ~= 0
        % pratical salinity (SP), Absolute salinity (SA) and conserv. T (CT)
        SP = SBS;
        [SA, in_ocean] = gsw_SA_from_SP(SP,P,profile_lon, profile_lat); %SA
        CT = gsw_CT_from_t(SA, prof_corr.temp, P); %CT
        % units
        SP_units = 'psu';
        SA_units = 'g/kg';
        CT_units = 'degC';        
        clear in_ocean prof_corr basicProfileData
    else
        disp('Lat Lon coordinates of the casts are not provided in CALIB.DAT')
        disp(['SP, SA and CT would not be calculated!! (press any ' ...
              'key to continue)'])
    end
    clear sbt sbc 
end


% Deconvolve temperature microstructure
if exist('T1_dT1') ~= 0 % variable exists?
    t1_nonphys = deconvolve('temperature', [], T1_dT1, fs, diff_gain_t1);
    t1_nonphys_units = 'non physical temperature';
    clear T1 T1_dT1
end
if exist('T2_dT2') ~= 0
    t2_nonphys = deconvolve('temperature', [], T2_dT2, fs, diff_gain_t2);
    t2_nonphys_units = 'non physical temperature';
    clear T2 T2_dT
end


% Deconvolve conductivity microstructure
if exist('C1_dC1') ~= 0
    c1_nonphys = deconvolve('conductivity', [], C1_dC1, fs, diff_gain_c1);
    c1_nonphys_units = 'non physical conductivity';
    clear C1_dC1
end
   
% Convert fluorescence and turbidity to volts and then to physical units.
if exist('ch14') ~= 0
    [fluoro, fluoro_units] = convert_odas(ch14,'fluoro','file',setup_fname);
    %    plot(fluoro);
% $$$     fluoro = A_fluoro + B_fluoro*fluoro*1000; % x 1000 to convert to mV.
    fluoro = A_fluoro + B_fluoro*fluoro; % x 1000 to convert to mV.
    fluoro_units = 'ppb';
    clear ch14
end
if exist('ch15') ~= 0
    [trans, trans_units]   = convert_odas(ch15,'trans','file',setup_fname);
% $$$     trans  = A_trans  + B_trans*trans*1000;   % x 1000 to convert to mV.
    trans  = A_trans  + B_trans*trans;   % x 1000 to convert to mV.
    trans_units = 'FTU';
    clear ch15
end



%clear rest of variables
clear tline sp_char setup_fname header gnd fid calib_fname
if exist('ans') == 0
    clear ans
end

% save the whole workspace after clearing unnecessary variables
save(data_fname) 

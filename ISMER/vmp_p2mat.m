function vmp_p2mat(data_fname);

% New version of older function 'p2mat.m'.
%
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
% ---------------------------------------------------------------------------- %


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
%       .maty file each time it works with a profile. (to speed up)

%  June 2013: F. Cyr introduced vmp_p2mat.m
%     Few changes are included:
%       - Now the function should deal with both header_version,
%       i.e., with setup.cfg in *.P header (ODAS V3). When new
%       odas V3 are used, CALIB.DAT and SETUP.TXT are not needed
%       anymore. If old file type are used (prior to 2013 in our
%       case), CALIB.DAT and SETUP.TXT should be present in the
%       current folder.
%       - added patch_odas prior to to fix_bad_buffer
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ----  fix_bad_buffers.m ---- %%%%
[bad_records, fix_manually] = patch_odas([data_fname]);

if isempty(bad_records) ~= 1
    [bad_records, chop, truncate] = fix_bad_buffers([data_fname, '.P']);
    if length(bad_records)>0
        disp('fix_bad_buffer.m had been applied!');
        disp(sprintf('number of points chopped at the begining of the file =%d', chop));
        disp(sprintf('number of points truncated at the end of the file =%d', truncate));
        disp(sprintf('number of points fixed by the method =%d', bad_records(end)));
        disp('Press any key to continue');
        pause
    end
    clear bad_records chop truncate
end

%%%% ---- read_odas.m (Version 3) ---- %%%%%
read_odas(data_fname); % convert binary data to .mat
load([data_fname,'.mat']); % load .mat just created 


%%%%% ---- Make time vectors ---- %%%%
FS = fs_slow;
fs = fs_fast;
fs_units = 'Hz'; FS_units = 'Hz';
nmax_slow = length(P); 
nmax_fast = length(sh1);
t_slow = (0:nmax_slow-1)'/FS; 
t_fast = (0:nmax_fast-1)'/fs;
clear nmax_slow nmax_fast fs_slow fs_fast

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


%%%% ---- Deconvolve and convert to phyisical units ---- %%%%

if header_version < 6 % Use legacy ODAS Library
                      
    % Look for SETUP.TXT file
    setup_fname = 'SETUP.TXT';
    fid = fopen(setup_fname);
    if fid == -1
        disp('No setup file was found, please check this! [Quit]')
        return
    end
    
    % get CALIB.DAT (Return if not present)
    calib_fname = 'CALIB.DAT';
    fid = fopen(calib_fname);
    if fid == -1
        disp('No CALIB.DAT was found, please check this! [Quit]')
        return
    else
        tline = fgetl(fid);
        while ischar(tline)
            eval(tline);
            tline = fgetl(fid);
        end
        fclose(fid);
    end

    % slow channels conversion
    Ptmp = deconvolve('pressure', P, P_dP, FS, diff_gain_P);
    [P, P_units] = convert_odas(Ptmp,'pres','file',setup_fname);
    
    if exist('sbt')~=0 & exist('sbc')~=0 
        [SBT, SBT_units] = convert_odas(sbt,'SBT1E','file',setup_fname);
        [SBC, SBC_units] = convert_odas(sbc,'SBC1E','file',setup_fname);
    end
    
    % vertical velocity (for shear)
    Nyquist = FS/2; 
    Wn      = 1.5/Nyquist; % W is smoothed at 1.5 Hz.
    W       = gradient(P, 1/FS);
    [b,a]   = butter(4,Wn);
    W       = filtfilt(b,a,W);
    w       = interp1(t_slow, W, t_fast,'linear','extrap');
    W_units = 'm/s'; w_units = 'm/s';
    clear a b Nyquist Wn
    
    % fast channels conversion
    [ax, ax_units] = convert_odas(Ax,'Pitch','file',setup_fname);
    [ay, ay_units] = convert_odas(Ay,'Roll','file',setup_fname);
    [az, az_units] = convert_odas(Az,'Az','file',setup_fname);

    % shear
    if exist('sh1')~=0
        shear1 = byte2shear(sh1,S1,G1,w);
        shear1_units = 's^{-1}'; 
    end
    
    if exist('sh2')~=0
        shear2 = byte2shear(sh2,S2,G2,w);
        shear2_units = 's^{-1}';
    end
    
    if exist('T1_dT1') ~= 0 
        t1_nonphys = deconvolve('temperature', [], T1_dT1, fs, diff_gain_t1);
        %[t1 t1_units] = convert_odas(t1_nonphys, 'therm', 'file', setup_fname);
        t1_nonphys_units = 'non physical temperature';
    end
    
    if exist('T2_dT2') ~= 0
        t2_nonphys = deconvolve('temperature', [], T2_dT2, fs, diff_gain_t2);
        t2_nonphys_units = 'non physical temperature';
    end

    if exist('C1_dC1') ~= 0
        c1_nonphys = deconvolve('conductivity', [], C1_dC1, fs, diff_gain_c1);
        c1_nonphys_units = 'non physical conductivity';
    end

    if exist('ch14') ~= 0
        [fluoro, fluoro_units] = convert_odas(ch14,'fluoro','file',setup_fname);
        fluoro = A_fluoro + B_fluoro*fluoro*1000; % *1000 to convert to mV
        fluoro_units = 'ppb';
    end
    if exist('ch15') ~= 0
        [trans, trans_units]   = convert_odas(ch15,'trans','file',setup_fname);
        trans  = A_trans  + B_trans*trans*1000;   
        trans_units = 'FTU';
    end

else % Use ODAS V3 
    
    % Get lat lon
    %profile_lat = str2num(char(inifile_with_instring(setupfilestr, 'readstring', {'cruise info','','profile_lat'})));
    %profile_lon = str2num(char(inifile_with_instring(setupfilestr, 'readstring', {'cruise info','','profile_lon'})));
    cfg = setupstr(setupfilestr);
    profile_lat = setupstr(cfg,'cruise info','profile_lat');
    profile_lon = setupstr(cfg,'cruise info','profile_lon');
    profile_lat = str2num(profile_lat{1});
    profile_lon = str2num(profile_lon{1});
    

    % convert slow channels
    P_hres       = deconvolve('pres', P, P_dP, FS, setupfilestr, header_version);
    [P, P_units] = convert_odas(P_hres, 'pres', 'string', setupfilestr, header_version);
    
    %if exist('sbt')~=0 & exist('sbc')~=0 
    if exist('SBT')~=0 & exist('SBC')~=0 
        [SBT, SBT_units] = convert_odas(SBT, 'sbt', 'string', setupfilestr, header_version);
        [SBC, SBC_units] = convert_odas(SBC, 'sbc', 'string', setupfilestr, header_version);
    end    
    
    % vertical velocity (for shear)
    Nyquist = FS/2; 
    Wn      = 1.5/Nyquist; % W is smoothed at 1.5 Hz.
    W       = gradient(P, 1/FS);
    [b,a]   = butter(4,Wn);
    W       = filtfilt(b,a,W);
    w       = interp1(t_slow, W, t_fast,'linear','extrap');
    W_units = 'm/s'; w_units = 'm/s';
    clear a b Nyquist Wn
    
    % convert fast channels
    [az, az_units] = convert_odas(az, 'Az', 'string', setupfilestr, header_version);
    [ax, ax_units] = convert_odas(ax, 'Ax', 'string', setupfilestr, header_version);
    [ay, ay_units] = convert_odas(ay, 'Ay', 'string', setupfilestr, header_version);
    
 
    if exist('sh1') ~= 0 
        [shear1, shear1_units] = convert_odas(sh1, 'shear1', 'string', setupfilestr, header_version);
    end
    
    if exist('sh2') ~= 0     
        [shear2, shear2_units] = convert_odas(sh2, 'shear2', 'string', setupfilestr, header_version);
    end
    
    if exist('t1_dt1') ~= 0 
        %t1_nonphys  = deconvolve('dtherm1', [], T1_dT1, fs, setupfilestr, header_version); 
        t1_nonphys  = deconvolve('dtherm1', [], t1_dt1, fs, setupfilestr, header_version); 
        %t1 = convert_odas(t1_nonphys, 'dtherm1', 'string', setupfilestr, header_version);
    end
    
    if exist('t2_dt2') ~= 0 
        %t2_nonphys  = deconvolve('dtherm2', [], T2_dT2, fs, setupfilestr, header_version);   
        t2_nonphys  = deconvolve('dtherm2', [], t2_dt2, fs, setupfilestr, header_version);   
        %t2 = convert_odas(t2_nonphys, 'dtherm2', 'string', setupfilestr, header_version)
    end
       
    if exist('c1_dc1') ~= 0
        %c1_nonphys  = deconvolve('ucond', [], C1_dC1, fs, setupfilestr, header_version); 
        c1_nonphys  = deconvolve('ucond', [], c1_dc1, fs, setupfilestr, header_version); 
        %c1 = convert_odas(c1_nonphys, 'ucond', 'string', setupfilestr, header_version);
    end
    
    %if exist('trans') ~= 0
    if exist('turb') ~= 0
        [trans, trans_units] = convert_odas(turb, 'turb', 'string', setupfilestr, header_version);
    end       
    
    if exist('fluoro') ~= 0
        [fluoro, fluoro_units] = convert_odas(fluoro, 'fluoro', 'string', setupfilestr, header_version);
    end
end
% Clean unnecessary variables
%clear P_dP Ptmp sh1 sh2 T1 T1_dT1 T2 T2_dT2 C1_dC1 ch14 ch15
clear P_dP Ptmp sh1 sh2 t1 t1_dt1 t2 t2_dt2 c1_dc1 ch14 ch15


%%%% ----- post-convert treatment ---- %%%%
% Interpolate the pressure P (slow sampling) to the fast sampling rate.
p =  interp1(t_slow, P, t_fast, 'linear','extrap'); % There<s a lag bet. P and p (see correction further)
p_units = 'dBar';

% Accel to degree
g = 9.81;             % Gravitational acceleration.
pitch = asind(ax/g);  % Convert ax in m/s^2 to degree.
roll  = asind(ay/g);  % Convert ay in m/s^2 to degree.
pitch_units = 'degree'; 
roll_units  = 'degree';

% Fine-scale conductivity, temperature and salinity
if exist('SBT')~=0 & exist('SBC')~=0 

    % correct Conductivity for thermal Lag
    basicProfileData.ptime = MTIME;
    basicProfileData.depth = P;
    basicProfileData.temp  = SBT;
    basicProfileData.cond  = SBC;
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
end


%clear rest of variables
clear tline sp_char setup_fname header gnd fid calib_fname
if exist('ans') == 0
    clear ans
end


%%%%% ----  save the whole workspace ---- %%%%
save(data_fname) 

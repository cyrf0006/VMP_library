%function p2mat_ADCP(data_fname,setup_fname,fname_ADCP,year,month,day,useADCP,sensitivity1,sensitivity2)
function p2mat_ADCP(data_fname,setup_fname,fname_ADCP,year,month,day,useADCP,sensitivity1,sensitivity2,Wp,std_Wp)

% Diff. gain for pressure deconvolution (from calibration report)
diff_gain_P=20.3;
diff_gain_T=1.0;

% Parameters for despiking
thresh_spike=7;
smooth_spike=0.5;
N_spike=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the .p file and write to a mat-file and load the data
read_odas([data_fname,'.p'],setup_fname);
load([data_fname,'.mat']);

matrix_rows  = header(2,31);
fast_columns = header(2,29);
slow_columns = header(2,30);

columns = slow_columns + fast_columns;
clock_rate = header(2,21) + header(2,22)/1000; 
fs = clock_rate / columns; 
FS = fs / matrix_rows;

% Make time vectors
nmax_slow=length(Pres); 
nmax_fast=length(az);
t_slow = (0:nmax_slow-1)'/FS; 
t_fast = (0:nmax_fast-1)'/fs;

% Pressure deconvolution
Pres = deconvolve('pressure', Pres, ch11, round(FS), diff_gain_P);
clear ch11;
save(data_fname,'Pres','-append'); 

% Convert to physical units and save
[var_list,units]=conv_odas(data_fname, setup_fname);
load([data_fname,'.mat']);

% Change variable names
P = Pres;
pitch = Pitch; clear Pitch;
roll = Roll; clear Roll;
pitch=asin(pitch/9.81)*180/pi;
roll=asin(roll/9.81)*180/pi;

% Change variable name
SBT = SBT1; clear SBT1;
SBC = SBC1; clear SBC1;

% Make complete time vector
hh = str2num(time(:,1:2)); 
mm = str2num(time(:,4:5)); 
ss = str2num(time(:,7:12)); 
MTIME_tmp = datenum(year,month,day,hh,mm,ss);
P_tmp = P(1:64:length(P));

% Interpolate at high sampling rate the pressure and the time
p =  interp1(t_slow, P, t_fast, 'linear','extrap');
mtime = interp1((P_tmp-mean(P_tmp))/std(P_tmp),(MTIME_tmp-mean(MTIME_tmp))/std(MTIME_tmp),(p-mean(P_tmp))/std(P_tmp),'linear','extrap');
mtime = mtime*std(MTIME_tmp)+mean(MTIME_tmp);
MTIME = interp1((P_tmp-mean(P_tmp))/std(P_tmp),(MTIME_tmp-mean(MTIME_tmp))/std(MTIME_tmp),(P-mean(P_tmp))/std(P_tmp),'linear','extrap');
MTIME = MTIME*std(MTIME_tmp)+mean(MTIME_tmp);

% Calculate the falling speed
W = gradient(P, 1/FS);

%% Adjusting the flow past the sensor using the ADCP measurements.
if useADCP == 1
    
  load(fname_ADCP); 
  w_ADCP = ADCP.vert_vel;
  ' WARNING - Something specific to 5 July 2007'
  w_ADCP = w_ADCP(:,16900:17800);
  
  Iw = find(isnan(w_ADCP) == 1);
  w_ADCP(Iw) = 0;

  [nmax_ADCP kmax_ADCP]=size(w_ADCP);

  % Filters specifications
  order = 4;
  low_pass_period = 15; 
  SamplingPeriod = (ADCP.mtime(2)-ADCP.mtime(1))*86400;
  Fs = 1/SamplingPeriod;
  Fn_low = 1/low_pass_period;

  % Vertical
  low_pass_dz = 2;
  Ks = 1/ADCP.config.cell_size;  % Sampling inmterval dz
  Kn = 1/low_pass_dz;

  [B1,A1] = butter(order,Fn_low/(Fs/2),'low');
  [B3,A3] = butter(order,Kn/(Ks/2),'low');
  
  % Remove the noise
  for k=1:kmax_ADCP   
    w_ADCP(:,k)=filtfilt(B1,A1,w_ADCP(:,k));
  end
  
  % Vertical smoothing
  for i=1:nmax_ADCP
    w_ADCP(i,:)=filtfilt(B3,A3,w_ADCP(i,:));
  end

  ' WARNING - Something specific to 5 July 2007'
  z0=1.02;  % Depth of the ADCP 
  ADCPmtime = ADCP.mtime(16900:17800);
  
  %WI=interp2((ADCP.mtime-mean(ADCP.mtime))/std(ADCP.mtime),ADCP.config.ranges+z0,w_ADCP,(MTIME-mean(ADCP.mtime))/std(ADCP.mtime),P,'linear');
  WI=interp2((ADCPmtime-mean(ADCPmtime))/std(ADCPmtime),ADCP.config.ranges+z0,w_ADCP,(MTIME-mean(ADCPmtime))/std(ADCPmtime),P,'linear');
  III=find(isnan(WI)==1); WI(III)=0.0;
  Wc = W+WI;

  %z0=1.02;  % Depth of the ADCP 
  %III=find(isnan(ADCP.vert_vel)==1);
  %ADCP.vert_vel(III)=0.0; clear III;
  %WI=interp2(ADCP.mtime,ADCP.config.ranges+z0,ADCP.vert_vel,MTIME,P,'linear');
  %III=find(isnan(WI)==1); WI(III)=0.0;
  %Wc=W+WI;
  
else

  Wc=W;
  k = find(abs(W-Wp) > 1.5*std_Wp);
  Wc(k) = Wp;
  
end
%%

[b,a] = butter(4,0.5/FS/2);
W = filtfilt(b,a,W);
Wc = filtfilt(b,a,Wc);

% Interpolate at high sampling rate
w =  interp1(t_slow, W, t_fast,'linear');
wc = interp1(t_slow, Wc, t_fast,'linear');
p =  interp1(t_slow, P, t_fast,'linear');

% ** SHEAR **
% Convert shear signal to physical units
sh1 = she2du('odas', ch8, wc, sensitivity1);
sh2 = she2du('odas', ch9, wc, sensitivity2);
%sh1c = she2du('odas', ch8, wc, sensitivity1);
%sh2c = she2du('odas', ch9, wc, sensitivity2);
clear ch8 ch9

% ** TEMPERATURE MICROSTRUCTURE **
ch_t1 = deconvolve('temperature', [], ch5, fs, diff_gain_T);
ch_t2 = deconvolve('temperature', [], ch7, fs, diff_gain_T);
clear ch5 ch7

% ** CONDUCTIVITY AND SALINITY MICROSTRUCTURE **
ch_c1 = deconvolve('conductivity', [], ch12, fs, diff_gain_T);
clear ch12

% ** SEA-BIRD SALINITY
SBS = salinity(P, SBT, SBC);

save(data_fname,'fs','FS','time','mtime','MTIME','P','p','W','Wc','w','wc','pitch','roll','az','sh1','sh2','ch_t1','ch_t2','ch_c1','SBC','SBS','SBT'); 
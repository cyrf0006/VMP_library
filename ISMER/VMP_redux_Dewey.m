
fname_in =  'CFL_03aDEC2007_000.mat';
fname_out = 'CFL_03aDEC2007_000_redux.mat';

year = 2007; month = 12; day = 03;

% Options for Dewey's scripts
zbin = 5.0;
iplt = [1 2 3 5.0];
  
% Parameters for despiking
thresh_spike = 5; smooth_spike = 0.5; N_spike = 5;

% Load the data
load(fname_in);

FS = round(FS);
fs = round(fs);
  
% Despike the shear signals
[sh1 spike_index1] = despike(sh1, thresh_spike, smooth_spike, fs, N_spike);
[sh2 spike_index2] = despike(sh2, thresh_spike, smooth_spike, fs, N_spike);

% Flag data for which pitch or roll > 10
ax = ang2acc(pitch); ay = ang2acc(roll);
J = find(pitch > 10 | roll > 10);
sh1(J)= NaN; sh2(J) = NaN;
    
[eps1,p_eps1] = epsilon_vmp(sh1,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin);
[eps2,p_eps2] = epsilon_vmp(sh2,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin);
    
mtime_eps  = interp1(P,MTIME,p_eps1,'linear');
    
  
save(fname_out,'mtime_eps','p_eps1','p_eps2','eps1','eps2');
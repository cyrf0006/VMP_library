% $$$ % ---- For Buoy's ADCP --- %
% $$$ 
% $$$ %load IML4_ADCP2010
% $$$ load IML4_ADCP2009
% $$$ 
% $$$ east_vel  = ADCP.east_vel;
% $$$ north_vel  = ADCP.north_vel;
% $$$ time_adcp =  ADCP.mtime;
% $$$ z = ADCP.config.ranges;
% $$$ w = ADCP.vert_vel;
% $$$ T_adcp = ADCP.temperature;
% $$$ P_adcp = ADCP.pressure;
% $$$ 
% $$$ 
% $$$ 
% $$$ % for 2009, remove last time step
% $$$ east_vel = east_vel(:,1:end-1);
% $$$ north_vel = north_vel(:,1:end-1);
% $$$ w = w(:,1:end-1);
% $$$ time_adcp(end) = [];
% $$$ T_adcp(end) = [];
% $$$ P_adcp(end) = [];
% $$$ 
% $$$ 
% $$$ clear ADCP
% $$$ 
% $$$ 
% $$$ % $$$ save ADCP_RIKI2010_vel.mat east_vel north_vel w
% $$$ % $$$ save ADCP_RIKI2010_PTzt.mat P_adcp T_adcp z time_adcp
% $$$ 
% $$$ save ADCP_RIKI2009_vel.mat east_vel north_vel w
% $$$ save ADCP_RIKI2009_PTzt.mat P_adcp T_adcp z time_adcp 
% $$$ 
% $$$ % ---------------------------- %



% ----- For moorings' ADCP ----- %
%  -- First file -- %
%load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_RIKI/MRIKIa.mat')
load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/N_M080a.mat')
east_vel1  = ADCP.east_vel;
north_vel1  = ADCP.north_vel;
mtime1 =  ADCP.mtime;
z = ADCP.config.ranges;
w1 = ADCP.vert_vel;
T_adcp1 = ADCP.temperature;
P_adcp1 = ADCP.pressure;

% Backscatter intensity
I1 = ADCP.intens(:,1,:); I1 = squeeze(I1);
I2 = ADCP.intens(:,2,:); I2 = squeeze(I2);
I3 = ADCP.intens(:,3,:); I3 = squeeze(I3);
I4 = ADCP.intens(:,4,:); I4 = squeeze(I4);
Imean1 = (I1 + I2 + I3 + I4)/4;

clear ADCP

% -- Second file -- %
%load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_RIKI/MRIKIb.mat')
load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/N_M080b.mat')

east_vel2  = ADCP.east_vel;
north_vel2  = ADCP.north_vel;
mtime2 =  ADCP.mtime;
w2 = ADCP.vert_vel;
T_adcp2 = ADCP.temperature;
P_adcp2 = ADCP.pressure;

% Backscatter intensity
I1 = ADCP.intens(:,1,:); I1 = squeeze(I1);
I2 = ADCP.intens(:,2,:); I2 = squeeze(I2);
I3 = ADCP.intens(:,3,:); I3 = squeeze(I3);
I4 = ADCP.intens(:,4,:); I4 = squeeze(I4);
Imean2 = (I1 + I2 + I3 + I4)/4;

clear ADCP

% -- merge 2 files -- %
east_vel  = [east_vel1 east_vel2];
north_vel  = [north_vel1 north_vel2];
w = [w1 w2];
Imean = [Imean1 Imean2];
time_adcp = [mtime1 mtime2];
T_adcp = [T_adcp1 T_adcp2];
P_adcp = [P_adcp1 P_adcp2];

% Correct for magnetic declination!
%[east_vel, north_vel] = rotate_vecd(east_vel, north_vel, 18.5);


save M_N080_vel east_vel north_vel w 
save M_N080_PTzt P_adcp T_adcp z time_adcp
save M_N080_int Imean I1 I2 I3 I4
%save M_RIKI_vel east_vel north_vel w 
%save M_RIKI_PTzt P_adcp T_adcp z time_adcp


clear east_vel1 east_vel2 north_vel1 north_vel2 w1 w2 Imean1 Imean2 ...
    mtime1 mtime2 T_adcp1 T_adcp2 P_adcp1 P_adcp2

% ----------------------------------------------- %
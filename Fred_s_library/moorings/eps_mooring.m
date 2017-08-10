function eps_mooring(vel_file, time_file, eps_files, dt, hab)
    
% Function that will compute Epsilon from U^3 parametrization. Will
% need Cd and Z0 from Cd_fit.m
% usage ex:
% eps_mooring('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/M_N080_vel.mat','~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/M_N080_PTzt.mat','hit_bottom_epsprofiles_slope', 5/60/24, 20)
%eps_mooring('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/M_N080_vel.mat','~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/M_N080_PTzt.mat','hit_bottom_epsprofiles_slope', 1/60/24, 20)%
% Will return: - timeVec (1xn)
%              - time2Vec (1xn)
%              - zhab (1xm)
%              - epsMat (mxn)
%  For manuscript:
% eps_mooring('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/M_N080_vel.mat','~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/M_N080_PTzt.mat','hit_bottom_epsprofiles_slope', 15/60/24, 20)

     
% few info:
kappa = .41;
% at z = 1m
C_d = 1.8e-3;
Z_0 = 0.0007;
% at z = 5m`
C_d = 1.7e-3;
Z_0 = 0.0003;
zUb = 5;

% Soulsby...
C_d = 3e-3;
zUb = 2;

%C_d = 3e-3; %Soulsby1997

% ------ ADCP VELOCITIES ------ %
load(vel_file);
load(time_file);

% find bottom ()
for i = 1:size(w, 1)
    I = find(isnan(w(i, :))==1);
    if length(I)>size(w, 2)./4
        bot = i;
        break
    end
end

% reduce file to desired depth range
% We calculate u,v just in case, but only U^2+V^2 will be used
u = east_vel(1:bot-1, :);
v = north_vel(1:bot-1, :);
[u,v] = rotate_vecd(v, u, 33.5);

%w = w(1:bot-1, :);
z = z(1:bot-1);
dz = z(2)-z(1);
zhab = z(end)-z+dz;


% Vector norm
U = sqrt(u.^2 + v.^2);


% reduce time range
tmin = datenum(2011, 09, 20, 0, 0, 0);
tmax = datenum(2011, 10, 12, 0, 0, 0);
I = find(time_adcp>tmin & time_adcp<tmax);
J = find(zhab<=hab);
time_adcp = time_adcp(I);
U = U(J,I);
u = u(J,I);
v = v(J,I);
zhab = zhab(J);

% Original binning presscribed in function input
% 5 min., 15 min., etc. bins
t1 = round(time_adcp(1)*24)/24;
t2 = round(time_adcp(end)*24)/24;
t_bin = t1:dt:t2;
for i = 1:length(t_bin);
    I = find(time_adcp >= t_bin(i)-dt/2 & time_adcp < t_bin(i)+dt/2);
    U_bin(:,i) = nanmean(U(:,I), 2);
    u_bin(:,i) = nanmean(u(:,I), 2);
    v_bin(:,i) = nanmean(v(:,I), 2);
end

% -------- MEAN SHEAR (Method 1) ------- %
% $$$ % 30 minutes bin for Shear
% $$$ t1 = round(time_adcp(1)*24)/24;
% $$$ t2 = round(time_adcp(end)*24)/24;
% $$$ t_bin = t1:dt:t2;
% $$$ for i = 1:length(t_bin);
% $$$     I = find(time_adcp >= t_bin(i)-dt/2 & time_adcp < t_bin(i)+dt/2);
% $$$     u_bin_tmp(:,i) = nanmean(u(:,I), 2);
% $$$     v_bin_tmp(:,i) = nanmean(v(:,I), 2);
% $$$ end
% vertical bins
% $$$ 
% $$$ dzS2 = 1;
% $$$ zhab_S2bin = dzS2/2:dzS2:max(zhab);
% $$$ for i = 1:length(zhab_S2bin);
% $$$     I = find(zhab >= zhab_S2bin(i)-dzS2/2 & zhab < zhab_S2bin(i)+dzS2/2);
% $$$     u_bin_S2(i,:) = nanmean(u_bin_tmp(I,:), 1);
% $$$     v_bin_S2(i,:) = nanmean(v_bin_tmp(I,:), 1);
% $$$ end

dzS2 = 1;
zhab_S2bin = dzS2/2:dzS2:max(zhab);
for i = 1:length(zhab_S2bin);
    I = find(zhab >= zhab_S2bin(i)-dzS2/2 & zhab < zhab_S2bin(i)+dzS2/2);
    u_bin_S2(i,:) = nanmean(u_bin(I,:), 1);
    v_bin_S2(i,:) = nanmean(v_bin(I,:), 1);
end

du = diff(u_bin_S2,1,1);
dv = diff(v_bin_S2,1,1);
% $$$ du = diff(u_unfilt,1,1);
% $$$ dv = diff(v_unfilt,1,1);
dz = dzS2;

S2Mat = (du/dz).^2 + (dv/dz).^2;
zhab_S2 = zhab_S2bin(1:end-1)+dz/2;
% ----------------------------- %


% $$$ % Low pass velocity (30 minutes) 
% $$$ fs = 1/(dt*86400);
% $$$ freq_low = 1/(30*60); %Hz
% $$$ Wn_low = freq_low/(fs/2);
% $$$ [b,a] = butter(4, Wn_low);
% $$$ 
% $$$ U_filt = nan(size(U_bin));
% $$$ u_filt = nan(size(U_bin));
% $$$ v_filt = nan(size(U_bin));
% $$$ 
% $$$ for i = 1:size(U_bin, 1)
% $$$     % remove NANs
% $$$     Uvec = U_bin(i,:);
% $$$     I = find(isnan(Uvec)==1);
% $$$     xitp = 1:length(Uvec);
% $$$     x = xitp;
% $$$     Uvec(I) = [];
% $$$     x(I) = [];
% $$$     Uvec = interp1(x, Uvec, xitp);    
% $$$     U_filt(i, :) = filtfilt(b, a, Uvec);
% $$$     
% $$$     uvec = u_bin(i,:);
% $$$     I = find(isnan(uvec)==1);
% $$$     xitp = 1:length(uvec);
% $$$     x = xitp;
% $$$     uvec(I) = [];
% $$$     x(I) = [];
% $$$     uvec = interp1(x, uvec, xitp);    
% $$$     u_filt(i, :) = filtfilt(b, a, uvec);
% $$$     
% $$$     vvec = v_bin(i,:);
% $$$     I = find(isnan(vvec)==1);
% $$$     xitp = 1:length(vvec);
% $$$     x = xitp;
% $$$     vvec(I) = [];
% $$$     x(I) = [];
% $$$     vvec = interp1(x, vvec, xitp);    
% $$$     v_filt(i, :) = filtfilt(b, a, vvec);
% $$$ end
% $$$ 
% $$$ % Rename variables
% $$$ u_unfilt = u_bin;
% $$$ v_unfilt = u_bin;

% $$$ U_bin = U_filt;
% $$$ u_bin = u_filt;
% $$$ v_bin = v_filt;

% Bottom velocity
[Y, I] = min(abs(zhab-zUb));
Ub = U_bin(I,:); %Ub

% $$$ % Velocity defect law?
% $$$ Ub= nanmean(U_bin(1,:),1)-nanmean(U_bin(I,:),1);

% ------- EPSILON FROM U^3 ------- %
epsMat = nan(size(U_bin));
for i = 1:length(Ub)
    %epsMat(:,i) = C_d.^(3/2)*Ub(i).^(3)./(kappa*zhab);
    epsMat(:,i) = C_d.^(3/2)*Ub(i).^(3)./(kappa*zhab);
end
%save epsilon_U3.mat eps zhab t_bin

% $$$ % ------- MEAN SHEAR (method 2) ------- %
% $$$ % vertical bins
% $$$ dzS2 = 1;
% $$$ zhab_S2bin = dzS2/2:dzS2:max(zhab);
% $$$ for i = 1:length(zhab_S2bin);
% $$$     I = find(zhab >= zhab_S2bin(i)-dzS2/2 & zhab < zhab_S2bin(i)+dzS2/2);
% $$$     u_bin_S2(i,:) = nanmean(u_ubin(I,:), 1);
% $$$     v_bin_S2(i,:) = nanmean(v_ubin(I,:), 1);
% $$$ end
% $$$ %remove NaNs
% $$$ 
% $$$ 
% $$$ du = diff(u_bin_S2,1,1);
% $$$ dv = diff(v_bin_S2,1,1);
% $$$ % $$$ du = diff(u_unfilt,1,1);
% $$$ % $$$ dv = diff(v_unfilt,1,1);
% $$$ 
% $$$ S2Mat = (du/dz).^2 + (dv/dz).^2;
% $$$ zhab_S2 = zhab_S2bin(1:end-1)+dz/2;
%S2Mat = S2_bin;
% ---------------------------------------- %


% -------------- TIME TO HIGHTIDE ------------ %
tide  = load('~/WINDEX/data_processing/BMix_study/tide_2009-2012.dat');
    
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);

% find high tide time 
count = 1;
for i = 2:length(mtime)-1

    if level(i)>level(i-1) & level(i)>level(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        count = count+1;
    end
end
clear mtime

for i = 1:length(t_bin)
    [Y, I] = min(abs(t_bin(i)-T));
    time2(i) = (t_bin(i)-T(I))*24;
end

timeVec = t_bin;


save epsMat_fromU3.mat epsMat S2Mat zhab zhab_S2 timeVec time2

% $$$ zhab_U_bin = zhab_S2bin;
% $$$ save mooringVel.mat U_bin_S2,  zhab_S2bin timeVec time2

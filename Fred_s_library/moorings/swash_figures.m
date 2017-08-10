function swash_figures(adcp_vel, adcp_sensors, tempFile, adcp_depth, drange,  tide_file, timeVec)

% 
% function swash_figures(adcp_vel, adcp_sensors, tempFile, adcp_depth, drange,  tide_file, timeVec)
%
% Written to have a idea of each M2 cycle that occured during the
% deployment of mooring N080. Pannel plot TUVS and their average
% profiles. Was called in
% /home/cyrf0006/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080
%    
% ex:
% swash_figures('M_N080_vel.mat', 'M_N080_PTzt.mat', 'Tfield_M_N080.mat',59, [59 83], '../tide_shear/tide_2009-2012.dat', datenum(2011,09,20,0,0,0):10/(24*60):datenum(2011,10,12,0,0,0)) 
%


theta =33.5;
M2Period = 12.42/24;

% ----------- temperature ------------ %
load(tempFile);
Thab = hab;



% ----------- ADCP ------------ %
load(adcp_vel)
load(adcp_sensors)

z_adcp = z + adcp_depth;
if z_adcp(1) < 0 % upward looking
    z_adcp = abs(z_adcp);
end

z_adcp_orig = z_adcp;
% vertical reduction
I = find(z_adcp > drange(1) & z_adcp < drange(2));
J = find(time_adcp > timeVec(1) & time_adcp < timeVec(end));
E = east_vel(I, J);
N = north_vel(I, J);
z_adcp = z_adcp(I);
timeADCP = time_adcp(J);

% time average
dt_adcp = timeADCP(2)-timeADCP(1);
dt = timeVec(2)-timeVec(1); % in days...

% Compute shear...
dz = z_adcp(2)-z_adcp(1);
du = diff(E, 1);
dv = diff(N, 1);
z_adcp_S = z_adcp(1:end-1)+dz/2;

S2_raw = (du./dz).^2; + (dv./dz).^2;

hab = abs(z_adcp-drange(end)); % we consider drange(end) is bottom
habS = abs(z_adcp_S-drange(end));

% -------------------------------- %


% ------------------ compute time relative to high tide  --------------------- % 
springtide = 3.5; %m
    
tide  = load(tide_file);
    
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
% restrict Hightide to mooring deployement
I = find(T>=timeVec(1) & T<=timeVec(end));
HT = T(I);


% Here's the main loop that will generate a figure for each M2 cycle
for i = 1: length(HT)
    
    I = find(time_vec_N080>= HT(i) - M2Period/2  & time_vec_N080 <= HT(i) + M2Period/2 +dt);
    T_tide = Tgrid(:,I);
    timeT = time_vec_N080(I);
    
    I = find(timeADCP >= HT(i) - M2Period/2 - dt & timeADCP <= HT(i) + M2Period/2 +dt);
    E_tmp = E(:,I);
    N_tmp = N(:,I);
    S_tmp = S2_raw(:,I);
    time_tmp = timeADCP(I);

    I = find(timeADCP >= HT(i) - M2Period/2 & timeADCP <= HT(i) + M2Period/2);
    time_tide = timeADCP(I(1)):dt:timeADCP(I(end));
    
    % bin (in time) this portion of the time serie
    E_tide = nan(size(E_tmp, 1), length(time_tide));
    N_tide = nan(size(E_tmp, 1), length(time_tide));
    S_tide = nan(size(S_tmp, 1), length(time_tide));
    for j = 1:length(time_tide)
        I = find(time_tmp>=time_tide(j)-dt/2 & time_tmp<=time_tide(j)+dt/2);
        E_tide(:,j) = nanmean(E_tmp(:,I), 2);
        N_tide(:,j) = nanmean(N_tmp(:,I), 2);
        S_tide(:,j) = nanmean(S_tmp(:,I), 2);
    end
    
    plot_swash
    figureName = ['HighTide_' datestr(HT(i), 29) '_' datestr(HT(i), 13)];
    disp(['saving ' figureName])
    print('-dpng', '-r300', figureName)
    
    
end

% $$$ clear time2_level  level2 LEV level2
% $$$ 
% $$$ % This compute the time and level versus closest high tide
% $$$ for i = 1:length(mtime)
% $$$     [Y, I] = min(abs(mtime(i) - T));
% $$$     time2_level(i) = (mtime(i)-T(I))*24; % time relative to the closest hightide
% $$$                                          % level2(i) = L(I); %level of the closest hightide
% $$$ end
% $$$ 
% $$$ % This compute the mean tide cycle using preceeding
% $$$ for i = 1: length(reg_tide) 
% $$$     I = find(time2_level > reg_tide(i) - dtide & time2_level < reg_tide(i) + dtide);
% $$$     size(I)
% $$$     LEV(i) = nanmean(level(I));
% $$$ end

% ---------- Plot Mean M2 cycle ------------- %
% $$$ t1 = datenum(2011, 09, 28, 9, 0, 0);
% $$$ t2 = datenum(2011, 09, 29, 0, 0, 0);
% $$$ t1 = datenum(2011, 09, 20, 0, 0, 0);
% $$$ %t2 = datenum(2011, 10, 12, 0, 0, 0);
% $$$ %t2 = datenum(2011, 10, 7, 0, 0, 0);
% $$$ t2 = datenum(2011, 10, 7, 0, 0, 0);


% $$$ timeVec = reg_tide;
% $$$ %level = LEV;
% $$$ mtime = reg_tide;
% $$$ S2_ave = log10(S_tide);
% $$$ E_ave = E_tide;
% $$$ N_ave = N_tide;
% $$$ t1 = timeVec(1);
% $$$ t2 = timeVec(end);
% $$$ 
% $$$ 
% $$$ 
% $$$ keyboard
% $$$ 
% $$$ 
% $$$ plot_shear_M2ave
% $$$ 
% $$$ print('-dpng', '-r300', 'TuvS_M2ave.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc', 'TuvS_M2ave.eps')

% ----------------------------------------------------- %






















function tide_shear_M2ave(adcp_vel, adcp_sensors, adcp_depth, drange,  tide_file, time_vec)

% 
% function tide_shear(adcp_vel, adcp_sensors, adcp_depth, drange,  tide_file, time_vec)
%
% ex:
% tide_shear_M2ave('../M_RIKI/M_RIKI_vel.mat','../M_RIKI/M_RIKI_PTzt.mat',-149, [70 149], 'tide_2009-2011.dat', datenum(2011,09,20,0,0,0):1/84:datenum(2011,10,08,0,0,0))
% tide_shear_M2ave('../M_RIKI/M_RIKI_vel.mat','../M_RIKI/M_RIKI_PTzt.mat',-149, [50 149], 'tide_2009-2011.dat', datenum(2011,09,20,0,0,0):1/84:datenum(2011,10,08,0,0,0))
% tide_shear_M2ave('../M_N080/M_N080_vel.mat','../M_N080/M_N080_PTzt.mat',59, [59 83], 'tide_2009-2012.dat', datenum(2011,09,20,0,0,0):15/1440:datenum(2011,10,12,0,0,0)) 
%
% -> Manuscript:
% tide_shear_M2ave('../M_N080/M_N080_vel.mat','../M_N080/M_N080_PTzt.mat',59, [59 83], 'tide_2009-2012.dat', datenum(2011,09,20,0,0,0):30/1440:datenum(2011,10,12,0,0,0)) 

% tide_shear_M2ave('/home/cyrf0006/PhD/IML4_ADCP/2010/ADCP_RIKI2010_vel.mat','/home/cyrf0006/PhD/IML4_ADCP/2010/ADCP_RIKI2010_PTzt.mat',2, [6 100], 'tide_2009-2011.dat', datenum(2010,04,21,0,0,0):1/48:datenum(2010,06,04,0,0,0)) 
% tide_shear_M2ave('/home/cyrf0006/PhD/IML4_ADCP/2009/ADCP_RIKI2009_vel.mat','/home/cyrf0006/PhD/IML4_ADCP/2009/ADCP_RIKI2009_PTzt.mat',2, [6 100], 'tide_2009-2011.dat', datenum(2009,05,01,0,0,0):1/24:datenum(2009,11,12,0,0,0)) 



% *_vel.mat et *_PTzt.mat come from 


theta =33.5;


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
J = find(time_adcp > time_vec(1) & time_adcp < time_vec(end));
E = east_vel(I, J);
N = north_vel(I, J);
z_adcp = z_adcp(I);
time_adcp = time_adcp(J);


% time average
dt_adcp = time_adcp(2)-time_adcp(1);
dt = time_vec(2)-time_vec(1); % in days...

% This WAS an error
% $$$ % Compute shear...
% $$$ dz = z_adcp(2)-z_adcp(1);
% $$$ du = diff(E, 1);
% $$$ dv = diff(N, 1);
% $$$ z_adcp_S = z_adcp(1:end-1)+dz;
% $$$ S2_raw = (du./dz).^2; + (dv./dz).^2;


clear S2_ave E_ave N_ave
% time average

for i = 1:length(time_vec)
    I = find(time_adcp > time_vec(i)-dt/2 &  time_adcp < time_vec(i)+dt/2);
    if ~isempty(I)==1
        %S2_ave(:,i) = nanmean(S2_raw(:,I), 2);
        E_ave(:,i) = nanmean(E(:,I), 2);
        N_ave(:,i) = nanmean(N(:,I), 2);
    else
        %S2_ave(:,i) = nan(size(S2_raw,1),1);
        E_ave(:,i) = nan(size(E,1),1);
        N_ave(:,i) = nan(size(N,1),1);
    end 
end


% Compute shear...
dz = z_adcp(2)-z_adcp(1);
du = diff(E_ave, 1, 1);
dv = diff(N_ave, 1, 1);
z_adcp_S = z_adcp(1:end-1)+dz;
S2_ave = (du./dz).^2; + (dv./dz).^2;

mean_S2 = nanmean(S2_ave, 1);
% -------------------------------- %


% ------------------ compute time relative to high tide  --------------------- % 
springtide = 3.5; %m
    
tide  = load(tide_file);
    
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);

count = 1;
for i = 2:length(mtime)-1

    if level(i)>level(i-1) & level(i)>=level(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        count = count+1;
    end
end


count2=1;
count3=1;
clear spring_count neap_count
for i = 1:length(time_vec)

    [Y, I] = min(abs(time_vec(i) - T));
    A(i) = (time_vec(i)-T(I))*24; % time relative to the closest hightide
    B(i) = L(I); %level of the closest hightide

    if L(I) >= springtide
 
        spring_count(count2) = i; % indices in spring
        count2 = count2+1;
    else
        neap_count(count3) = i; % indices in neap
        count3 = count3+1;
    end
    
end

time2 = A;

% $$$ % to consider just spring
% $$$ time2 = A(spring_count);
% $$$ mean_S2 = mean_S2(spring_count);

% $$$ % to consider just neap
% $$$ time2 = A(neap_count);
% $$$ mean_S2 = mean_S2(neap_count);


dtide = .25;
reg_tide = -6:dtide:6;

E_tide = nan(length(z_adcp), length(reg_tide));
N_tide = nan(length(z_adcp), length(reg_tide));
S_tide = nan(length(z_adcp)-1, length(reg_tide));

for i = 1: length(reg_tide)
% $$$     if i == 1
% $$$         I = find(time2 < (reg_tide(i) + dtide/2));  
% $$$     elseif i == length(reg_tide)
% $$$         I = find(time2 > (reg_tide(i) - dtide/2));
% $$$     else            
% $$$         I = find(time2 > (reg_tide(i) - dtide/2) & time2 < (reg_tide(i) + dtide/2));  
% $$$     end
    I = find(time2 > (reg_tide(i) - dtide/2) & time2 < (reg_tide(i) + dtide/2));  

    E_tide(:,i) = nanmean(E_ave(:,I), 2);
    N_tide(:,i) = nanmean(N_ave(:,I), 2);
    S_tide(:,i) = nanmean(S2_ave(:,I), 2);
end

% ---------- Plot Mean M2 cycle ------------- %
% $$$ t1 = datenum(2011, 09, 28, 9, 0, 0);
% $$$ t2 = datenum(2011, 09, 29, 0, 0, 0);
% $$$ t1 = datenum(2011, 09, 20, 0, 0, 0);
% $$$ %t2 = datenum(2011, 10, 12, 0, 0, 0);
% $$$ %t2 = datenum(2011, 10, 7, 0, 0, 0);
% $$$ t2 = datenum(2011, 10, 7, 0, 0, 0);


time_vec = reg_tide;
%level = LEV;
mtime = reg_tide;
S2_ave = log10(S_tide);
E_ave = E_tide;
N_ave = N_tide;
t1 = time_vec(1);
t2 = time_vec(end);






plot_shear_M2ave

print('-dpng', '-r300', 'TuvS_M2ave.png')
set(gcf, 'renderer', 'painters')
print('-depsc2', 'TuvS_M2ave.eps')

% ----------------------------------------------------- %


keyboard

function tide_shear(adcp_vel, adcp_sensors, adcp_depth, drange,  tide_file, time_vec)

% 
% function tide_shear(adcp_vel, adcp_sensors, adcp_depth, drange,  tide_file, time_vec)
%
% ex:
% tide_shear('../M_RIKI/M_RIKI_vel.mat','../M_RIKI/M_RIKI_PTzt.mat',-149, [70 149], 'tide_2009-2011.dat', datenum(2011,09,20,0,0,0):1/84:datenum(2011,10,08,0,0,0))
% tide_shear('../M_RIKI/M_RIKI_vel.mat','../M_RIKI/M_RIKI_PTzt.mat',-149, [50 149], 'tide_2009-2011.dat', datenum(2011,09,20,0,0,0):1/84:datenum(2011,10,08,0,0,0))
% tide_shear('../M_N080/M_N080_vel.mat','../M_N080/M_N080_PTzt.mat',59, [59 83], 'tide_2009-2011.dat', datenum(2011,09,20,0,0,0):1/84:datenum(2011,10,12,0,0,0)) 
%
% tide_shear('/home/cyrf0006/PhD/IML4_ADCP/2010/ADCP_RIKI2010_vel.mat','/home/cyrf0006/PhD/IML4_ADCP/2010/ADCP_RIKI2010_PTzt.mat',2, [6 100], 'tide_2009-2011.dat', datenum(2010,04,21,0,0,0):1/48:datenum(2010,06,04,0,0,0)) 
% tide_shear('/home/cyrf0006/PhD/IML4_ADCP/2009/ADCP_RIKI2009_vel.mat','/home/cyrf0006/PhD/IML4_ADCP/2009/ADCP_RIKI2009_PTzt.mat',2, [6 100], 'tide_2009-2011.dat', datenum(2009,05,01,0,0,0):1/24:datenum(2009,11,12,0,0,0)) 



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

% Compute shear...
dz = z_adcp(2)-z_adcp(1);
du = diff(E, 1);
dv = diff(N, 1);
z_adcp_S = z_adcp(1:end-1)+dz;

S2_raw = ((du./dz).^2 + (dv./dz).^2);


clear S2_ave E_ave N_ave
% time average

for i = 1:length(time_vec)
    I = find(time_adcp > time_vec(i)-dt/2 &  time_adcp < time_vec(i)+dt/2);
    if ~isempty(I)==1
        S2_ave(:,i) = nanmean(S2_raw(:,I), 2);
        E_ave(:,i) = nanmean(E(:,I), 2);
        N_ave(:,i) = nanmean(N(:,I), 2);
    else
        S2_ave(:,i) = nan(size(S2_raw,1),1);
        E_ave(:,i) = nan(size(E,1),1);
        N_ave(:,i) = nan(size(N,1),1);
    end 
end
mean_S2 = nanmean(S2_ave, 1);
% -------------------------------- %


% $$$ % Little test with IML4_ADCP_2009
% $$$ load S2_mat
% $$$ mean_S2= S2_ave;
% $$$ time_vec = TT;

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

% to consider just neap
time2 = A(neap_count);
mean_S2 = mean_S2(neap_count);

% regular tide vector (-6, 6)
reg_tide = unique(round(time2));

reg_S2 = nan(round(length(mean_S2)/2), length(reg_tide)); % /2 its just approx

for i = 1: length(time2)
    [Y, I] = min(abs(time2(i)-reg_tide));
    J = find(isnan(reg_S2(:, I))==1);
    reg_S2(J(1), I)=mean_S2(i);
end
% ---------------------------------------------------------------------------- % 


% --------------------- LOESS Statistical test ------------------------ %
II = find(isnan(mean_S2)==1);
mean_S2(II)=[]; time2(II)=[];
%smooth_shear2 = loess(time2, mean_S2, sort(time2), 0.3, 2);
smooth_shear1 = loess(time2, mean_S2, sort(time2), 0.3, 1);
time_shear = sort(time2);

%load eps_tide_raw_bndry.mat
%load eps_tide_raw_bndry_spring.mat
load eps_tide_raw_bndry_neap.mat
%load eps_tide_raw_riki.mat
%load eps_tide_raw_riki_neap.mat
%load eps_tide_raw_riki_spring.mat

%smooth_eps2 = loess(time2, EPS,sort(time2), 0.3, 2);
smooth_eps1 = loess(time2, EPS,sort(time2), 0.3, 1);
time_eps = sort(time2);

%save neap time_shear smooth_shear1 time_eps smooth_eps1
% --------------------------------------------------------------------- %

% ------------------ Bootstrap --------------------- %
nboot = 500;
for i = 1:length(reg_tide)
    I = find(~isnan(reg_S2(:,i)));
    samp = reg_S2(I,i);
    NN = length(samp);
    clear S2_boot_b
    % create random sampling
    for b = 1:nboot
        r = rand(NN,1);
        r = ceil(r*NN/1);
        S2_boot_b(b) = nanmean(samp(r));
    end

    % mean of random sampling
    S2_boot_dot(i) = nanmean(S2_boot_b);

    % compute error (EFRON & GONG 83)
    S2_error(i) = sqrt(sum((diff([S2_boot_b S2_boot_dot(i)],1, 2)).^2, 2)./(nboot-1));

end

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 8])   

subplot('position', [.15 .15 .8 .35])
errorbar(reg_tide, S2_boot_dot, S2_error, '.k')
hold on
plot(time_shear, smooth_shear1, 'k', 'linewidth', 2)
%plot(time_shear, smooth_shear2, 'b', 'linewidth', 2)
hold off
xlim([-6.5 6.5])
ylabel('S^2 (m^2 s^{-2})')
xlabel('time versus high tide (hour)')
%set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'fontsize', 10)
%set(gca, 'yscale', 'log')
%ylim([.082 .09])

% $$$ keyboard
% $$$ load('N2_tide_bndry_spring.mat')

subplot('position', [.15 .6 .8 .35])
%load('eps_tide_riki.mat')
%load('eps_tide_bndry_spring.mat')
load('eps_tide_bndry_neap.mat')
%load('eps_tide_bndry.mat')
%load('eps_tide_riki_spring.mat')
%load('eps_tide_riki_neap.mat')

errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
set(gca, 'yscale', 'log')
hold on
plot(time_eps, smooth_eps1, 'k', 'linewidth', 2)
%plot(time_eps, smooth_eps2, 'b', 'linewidth', 2)
hold off
ylabel('\epsilon (W kg^{-1})')
xlim([-6.5 6.5])
ylim([1e-9 1e-6])
set(gca, 'xticklabel', [])
%xlabel('time versus high tide (hour)')
set(gca, 'xgrid', 'on')
set(gca, 'fontsize', 10)

print('-dpng', '-r300', 'S2_vs_M2.png')
set(gcf, 'renderer', 'painters')
print('-depsc', '-r300', 'S2_vs_M2.png')


keyboard
save subplot_bndry.mat reg_tide S2_boot_dot S2_error time_shear smooth_shear1 eps_boot_dot eps_error time_eps smooth_eps1
save subplot_bndry_neap.mat reg_tide S2_boot_dot S2_error time_shear smooth_shear1 eps_boot_dot eps_error time_eps smooth_eps1
save subplot_bndry_spring.mat reg_tide S2_boot_dot S2_error time_shear smooth_shear1 eps_boot_dot eps_error time_eps smooth_eps1
save subplot_riki_spring.mat reg_tide S2_boot_dot S2_error time_shear smooth_shear1 eps_boot_dot eps_error time_eps smooth_eps1
% ---------------------------------------------------- %




% ---------- Shear at particular moment ------------- %
% $$$ t1 = datenum(2011, 09, 28, 9, 0, 0);
% $$$ t2 = datenum(2011, 09, 29, 0, 0, 0);
t1 = datenum(2011, 09, 20, 0, 0, 0);
%t2 = datenum(2011, 10, 12, 0, 0, 0);
%t2 = datenum(2011, 10, 7, 0, 0, 0);
t2 = datenum(2011, 10, 7, 0, 0, 0);

%plot_shear_N080
plot_shear_RIKI

% ----------------------------------------------------- %



keyboard

% ------------ Residual circulation ------------ %
% $$$ t1 = datenum(2011, 09, 20, 0, 0, 0);
% $$$ t2 = datenum(2011, 10, 12, 0, 0, 0);
% $$$ I = find(time_adcp >= t1 & time_adcp <= t2);
% $$$ [valong, vcross]=rotate_vecd(east_vel(:,I), north_vel(:,I), theta);

% positive is downstream....
[valong, vcross]=rotate_vecd(E, N, theta);

res_along = nanmean(valong, 2);
res_cross = nanmean(vcross, 2);

figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 7 10])   
plot(res_along*100, z_adcp, '--b')
hold on
plot(res_cross*100, z_adcp, '--r')
set(gca, 'ydir', 'reverse')
xlabel('u,v (cm s^{-1})')
ylabel('z (m)')
legend('u', 'v', 'location', 'southwest')
%print('-dpng', '-r300', 'uv_MRIKI.png')
print('-dpng', '-r300', 'uv_MN080.png')

% ----------------------------------------------- %

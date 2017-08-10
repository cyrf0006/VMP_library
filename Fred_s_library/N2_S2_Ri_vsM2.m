clear

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 5; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

% load N2 info
load(['/home/cyrf0006/WINDEX/data_processing/batchelor/' ...
      'chi_hit_bottom/N2_M2.mat'])

% load S2 info
load(['~/WINDEX/data_processing/sept_2011_mission/Mouillages/' ...
      'tide_shear/subplot_bndry.mat'])


% load chi info
load(['/home/cyrf0006/WINDEX/data_processing/batchelor/' ...
      'chi_hit_bottom/chi_M2.mat']);

% load IWs info
load(['/home/cyrf0006/WINDEX/data_processing/sept_2011_mission/' ...
      'Mouillages/M_N080/Ea_M2.mat'])


% Compute Stability index
Ri = N2_ave(1:13)./S2_boot_dot ;
Ri_error = abs(Ri - (N2_ave(1:13)+N2_error(1:13))/(S2_boot_dot+S2_error));

N2_itp = interp1(time_N2, smooth_N2, time_shear);
smooth_Ri = N2_itp./smooth_shear1;



figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 15])

%subplot(311)
subplot(411)
errorbar(reg_tide, N2_ave(1:13), N2_error(1:13), '.k')
ylabel('N^2 (s^{-2})')
set(gca, 'xticklabel', [])
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
xlim([-6.5 6.5]) 
ylim([7e-5 3e-4]) 
hold on
plot(time_N2, smooth_N2, 'k', 'linewidth', 2)
adjust_space


%subplot(312)
subplot(512)
errorbar(reg_tide, S2_boot_dot, S2_error, '.k')
ylabel('S^2 (s^{-2})')
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
xlim([-6.5 6.5]) 
%ylim([3e-4 4.25e-4]) 
hold on
plot(time_shear, smooth_shear1, 'k', 'linewidth', 2)
adjust_space

%subplot(313)
subplot(513)
errorbar(reg_tide, Ri, Ri_error, '.k')
ylabel('N^2/S^2')
%xlabel('time relative to high tide')
hold on
plot(time_shear, smooth_Ri, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
%ylim([0.1 0.9]) 
adjust_space

subplot(514)
errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
ylabel('\epsilon (W kg^{-1})')
%xlabel('time relative to high tide')
hold on
plot(time_eps, smooth_eps1, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
set(gca, 'xticklabel', [])
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
adjust_space

subplot(515)
errorbar(reg_tide_chi, chi_boot_dot, chi_error, '.k')
ylabel('\chi (^{\circ}C^2 s^{-1})')
xlabel('time relative to high tide')
hold on
plot(time_chi, smooth_chi, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
adjust_space

print('-dpng', '-r300',  'N2_S2_Ri.png') % no colorbar



figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
subplot(311)
errorbar(reg_tide, N2_ave(1:13), N2_error(1:13), '.k')
ylabel('N^2 (s^{-2})')
set(gca, 'xticklabel', [])
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
xlim([-6.5 6.5]) 
ylim([7e-5 3e-4]) 
hold on
plot(time_N2, smooth_N2, 'k', 'linewidth', 2)
adjust_space


subplot(312)
errorbar(reg_tide, S2_boot_dot, S2_error, '.k')
ylabel('S^2 (s^{-2})')
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
xlim([-6.5 6.5]) 
%ylim([3e-4 4.25e-4]) 
hold on
plot(time_shear, smooth_shear1, 'k', 'linewidth', 2)
adjust_space

subplot(313)
errorbar(reg_tide, Ri, Ri_error, '.k')
ylabel('N^2/S^2')
%xlabel('time relative to high tide')
hold on
plot(time_shear, smooth_Ri, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 

set(gca, 'xgrid', 'on')
xlabel('time relative to high tide (hour)')
%ylim([0.1 0.9]) 
adjust_space
print('-dpng', '-r300',  'sub1.png') % no colorbar




figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
subplot(311)
errorbar(reg_tide, Ri, Ri_error, '.k')
ylabel('N^2/S^2')
%xlabel('time relative to high tide')
hold on
plot(time_shear, smooth_Ri, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
%ylim([0.1 0.9]) 
adjust_space

subplot(312)
errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
ylabel('\epsilon (W kg^{-1})')
%xlabel('time relative to high tide')
hold on
plot(time_eps, smooth_eps1, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
set(gca, 'xticklabel', [])
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
adjust_space

subplot(313)
errorbar(reg_tide_chi, chi_boot_dot, chi_error, '.k')
ylabel('\chi (^{\circ}C^2 s^{-1})')
xlabel('time relative to high tide (hour)')
hold on
plot(time_chi, smooth_chi, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
adjust_space

print('-dpng', '-r300',  'sub2.png') % no colorbar


% Epsilon vs Ea
figure(4)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 6])

plot(time_eps, smooth_eps1, 'k', 'linewidth', 2)
% $$$ hold on
% $$$ % load S2 info
% $$$ load(['~/WINDEX/data_processing/sept_2011_mission/Mouillages/' ...
% $$$       'tide_shear/subplot_bndry_spring.mat'])
% $$$ plot(time_eps, smooth_eps1, 'r', 'linewidth', 2)
% $$$ load(['~/WINDEX/data_processing/sept_2011_mission/Mouillages/' ...
% $$$       'tide_shear/subplot_bndry_neap.mat'])
% $$$ plot(time_eps, smooth_eps1, 'b', 'linewidth', 2)
% $$$ hold off
%legend('all', 'spring', 'neap', 'location', 'southeast')
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
ylabel('\epsilon (W kg^{-1})')
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
xlabel('time relative to high tide (hour)')

print('-dpng', '-r300',  'eps_Ea1.png') % no colorbar

figure(5)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 6])

plot(time_Ea_spring, Ea_smooth_spring, 'k', 'linewidth', 2)
ylabel('E_a (J m^{-2})')
% $$$ hold on
% $$$ plot(time_Ea_spring, Ea_smooth_spring, 'r', 'linewidth', 2)
% $$$ plot(time_Ea_neap, Ea_smooth_neap, 'b', 'linewidth', 2)
% $$$ hold off
% $$$ legend('all', 'spring', 'neap', 'location', 'southwest')
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
%set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
%set(gca, 'xticklabel', [])
xlabel('time relative to high tide (hour)')

print('-dpng', '-r300',  'eps_Ea2.png') % no colorbar




% PECS PLOT
figure(6)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
subplot(311)
errorbar(reg_tide, N2_ave(1:13), N2_error(1:13), '.k')
ylabel('N^2 (s^{-2})')
set(gca, 'xticklabel', [])
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
xlim([-6.5 6.5]) 
ylim([7e-5 3e-4]) 
hold on
plot(time_N2, smooth_N2, 'k', 'linewidth', 2)
adjust_space

subplot(312)
errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
ylabel('\epsilon (W kg^{-1})')
%xlabel('time relative to high tide')
hold on
plot(time_eps, smooth_eps1, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
set(gca, 'xticklabel', [])
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
adjust_space

subplot(313)
errorbar(reg_tide_chi, chi_boot_dot, chi_error, '.k')
ylabel('\chi (^{\circ}C^2 s^{-1})')
xlabel('time relative to high tide (hour)')
hold on
plot(time_chi, smooth_chi, 'k', 'linewidth', 2)
xlim([-6.5 6.5]) 
%ylim([0.1 0.9]) 
set(gca, 'yscale', 'log') 
set(gca, 'xgrid', 'on')
adjust_space

print('-dpng', '-r300',  'pecs_sub.png') % no colorbar

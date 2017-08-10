% head_profiles.m:
%
% Script that aims to hughlught the difference between profiles
% over the sill in the deep portion of the Head, before the sill.
% Will also plot Nox_S relationship.
% 
% NOTE: Will need that 
% NOx_relationship('datFiles.list', 50, 300, 'NOx', 'S') is run
% before calling the script so plot_NOxRel_info.mat is created.
%
% to be run in anywhere
%
% F. Cyr - Feb. 2014

load('/home/cyrf0006/PhD/Nitrates/allStat/plot_NOxRel_info.mat')


prof1 = '~/WINDEX/data_processing/Mission_tadoussac_2009/profile_2009-10-01_001.mat';
prof2 = '~/WINDEX/data_processing/Mission_tadoussac_2009/profile_2009-10-01_068.mat';

load(prof1)
SA1 = SA;
CT1 = CT;
P1 = P;
load(prof2)
SA2 = SA;
CT2 = CT;
P2 = P;
rho1 = gsw_rho(SA1,CT1,P1);
rho2 = gsw_rho(SA2,CT2,P2);
[Y1, I1] = sort(rho1);
[Y2, I2] = sort(rho2);

% test sort salinity
[Y1, I1] = sort(SA1);
[Y2, I2] = sort(SA2);

% $$$ 
% $$$ 
% $$$ plot(rho1, P1, 'k', 'linewidth', 2)
% $$$ hold on
% $$$ plot(rho2, P2, 'r', 'linewidth', 2)
% $$$ set(gca, 'ydir', 'reverse')
% $$$ adjust_space
% $$$ 
% $$$ subplot(2,2,4)
% $$$ 
% $$$ plot(CT1, P1, 'k', 'linewidth', 2)
% $$$ hold on
% $$$ plot(CT2, P2, 'linewidth', 2)
% $$$ set(gca, 'ydir', 'reverse')
% $$$ adjust_space
% $$$ 
% $$$ print(gcf, '-dpng', 'head_profiles1.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print(gcf, '-depsc2', 'head_profiles1.eps')




figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 10])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.1 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.13; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.11; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


subplot(1,2,1)

plot(xVec, yVec, '.r')
hold on
plot(xVec3, yVec3, '.k')
plot(xVec2, yVec2, '.m')
plot(xVec3, yVec3, '.k')
hold off
legend('0-50m', '50-300m', '>300m', 'location', 'NorthWest')
legend('0-32 g kg^{-1}', '32-34.3 g kg^{-1}', '>34.3 g kg^{-1}', 'location', 'NorthWest')


[xVecRaw, I] = sort(xVec3);
yVecRaw = yVec3(I);


p = polyfit(xVecRaw, yVecRaw, 1);

%xVecItp = min(xVecRaw):max(xVecRaw);
yVecItp = p(2)+p(1).*xVecRaw;
hold on
plot(xVecRaw, yVecItp, 'color', [1 1 1]*.5)

[r,p] = corrcoef(yVecRaw,yVecItp);

xlim([min(xVec) 35])
ylim([min(yVec) max(yVec)])
xlabel('S_A (g kg^{-1})')
ylabel('C_{NO_3} (mmol m^{-3})')
text(34, 1.5, 'a', 'fontsize', 12, 'fontweight', 'bold')
adjust_space


subplot(1, 2, 2)
% $$$ plot(rho1, P1, 'k', 'linewidth', 2)
% $$$ hold on
% $$$ plot(rho2, P2, '--k', 'linewidth', 2)
NO31 = 7.3*SA1(I1) - 226.5;%32-34.3 SA
I = find(SA1(I1)<32 | SA1(I1)>34.3);
NO31(I) = NaN;

NO32 = 7.3*SA2(I2) - 226.5; %32-34.3 SA
I = find(SA2(I2)<32 | SA2(I2)>34.3);
NO32(I) = NaN;

plot(SA1(I1), P1, 'k', 'linewidth', 2)
hold on
plot(SA2(I2), P2, '-k', 'linewidth', 1)
set(gca, 'ydir', 'reverse')
ylabel('Depth (m)')
%xlabel('\rho (kg m^{-3})')
%xlim([1019 1029])
xlabel('S_A (g kg^{-1})')
xlim([25 34])
ylim([0 125])
adjust_space

% 2nd axis system
hAxes = gca;
hAxes_pos = get(hAxes,'Position');
hAxes_pos(4) = hAxes_pos(4);
hAxes2 = axes('Position',hAxes_pos);
% $$$ plot(CT1, P1, 'color', [0 .7 .7], 'linewidth', 2) % temp
% $$$ hold on
% $$$ plot(CT2, P2, 'linewidth', 2, 'color', [0 .7 .7], 'linestyle', '--')

plot(NO31, P1, 'color', [0 .7 .7], 'linewidth', 2) % NO3
hold on
plot(NO32, P2, 'linewidth', 1, 'color', [0 .7 .7], 'linestyle', '-')
set(hAxes2,'XAxisLocation','top', 'Color','none', 'YTickLabel',[])
h1_ylim = get(hAxes,'YLim'); % store x-axis limits of first axes
set(hAxes2,'YLim',h1_ylim) % specify x-axis limits of second axes
set(hAxes2,'xcolor',[0,.7,0.7])
%set(hAxes2,'XLim',[0 6.5]) % specify y-axis limits of second axes
set(hAxes2,'XLim',[6 17]) % specify y-axis limits of second axes

set(hAxes2, 'box', 'off')
set(hAxes2, 'ydir', 'reverse')
set(hAxes2, 'tickdir', 'out')
xpos = xlabel('C_{NO_3} (mmol m^{-3})');
XPOS = get(xpos, 'pos');
XPOS(2) = -7;
set(xpos, 'pos', XPOS) 
text(0.25, 265, 'b', 'fontsize', 12, 'fontweight', 'bold')




% $$$ subplot(2,2,4)
% $$$ 
% $$$ plot(CT1, P1, 'k', 'linewidth', 2)
% $$$ hold on
% $$$ plot(CT2, P2, 'linewidth', 2)
% $$$ set(gca, 'ydir', 'reverse')
% $$$ adjust_space
% $$$ 
print(gcf, '-dpng', '-r300', 'head_profiles.png')
set(gcf, 'renderer', 'painters')
print(gcf, '-depsc2', 'head_profiles.eps')



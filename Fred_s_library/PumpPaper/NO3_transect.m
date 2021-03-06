% script sels_transect.m
% run in /home/cyrf0006/PhD/Nitrates/transect
% see also ~/matlab/inPath/VMP_Library/Fred_s_library/o2/boot_NO3_all_stat.m
%    for .dat file generation
% See ecnlosed function below
clear 
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.05; % top of figure
bots = 0.13; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

[NO3field, xbin, zbin, dist] = getTransect('station_file_NO3_scatter');
%[NO2field] = getTransect('station_file_NO2');
%[PO4field] = getTransect('station_file_PO4');
stats = ([16:25]);

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 20 10])

V=[0:.5:25];
contourf(xbin, zbin, NO3field, V, 'linestyle', 'none')
%pcolor(xbin, zbin, NO3field)
%shading interp
hold on
contour(xbin, zbin, NO3field, [5:5:20 22 24], 'color', [1 1 1]*1, 'ShowText','on');

ylabel('Depth (m)')
ylim([0 425])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
    text(dist(i), -5, num2str(stats(i)), 'vertical', 'bot', 'horizontal', 'center');
end

% Add topography 
load ~/data/GSL_transect.mat
x_trans = x_trans-212;

x = [x_trans x_trans(end) x_trans(1) x_trans(1)];
y = [H_trans max(zbin) max(zbin) H_trans(1)];

patch(x, y, [1 1 1]*.0, 'edgecolor', 'none');
xlim([-25 700])

LSLE = [-20 190];
GSL = [190 800];
Yvert = 415;
plot(LSLE, [Yvert Yvert], '-', 'color', [1 1 1], 'linewidth', 2)
plot(GSL, [Yvert Yvert], '-', 'color', [1 1 1], 'linewidth', 2)
text(190, Yvert, '|', 'color', [1 1 1], 'fontSize', 15, 'verticalAlignment', ...
     'middle', 'fontWeight', 'bold')
text(-20, Yvert, '|', 'color', [1 1 1], 'fontSize', 15, 'verticalAlignment', ...
     'middle', 'horizontalAlignment', 'right', 'fontWeight', 'bold')
text(75, Yvert, 'LSLE', 'color', [1 1 1], 'fontSize', 12, 'verticalAlignment', ...
     'bottom', 'horizontalAlignment', 'center')
text(450, Yvert, 'GSL', 'color', [1 1 1], 'fontSize', 12, 'verticalAlignment', ...
     'bottom', 'horizontalAlignment', 'center')
% $$$ cb = colorbar;
% $$$ ti = ylabel(cb, '[NO_x] (mmol m^{-3})');

% $$$ cbPos = get(cb, 'pos');
% $$$ cbPos(1) = .01+cbPos(1);
% $$$ cbPos(2) = cbPos(2)+.1;
% $$$ cbPos(3) = .75*cbPos(3);
% $$$ cbPos(4) = cbPos(4)-.2;
% $$$ set(cb, 'pos', cbPos)
xlabel('Distance from Stat. 25 (km)');
adjust_space

% save transect info (for sketch)
save transect_info.mat xbin zbin NO3field

%save figure
print('-dpng', 'NOxTransect.png')
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'NOxTransect.eps')




%%% Second figure, With plot
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 4; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.1; % top of figure
bots = 0.13; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 20 18])


subplot(4,1,[1 3])
V=[0:.5:25];
contourf(xbin, zbin, NO3field, V, 'linestyle', 'none')
%pcolor(xbin, zbin, NO3field)
%shading interp
hold on
contour(xbin, zbin, NO3field, [5:5:20 22 24], 'color', [1 1 1]*1, 'ShowText','on');

ylabel('Depth (m)')
ylim([0 425])
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'tickdir', 'out')
for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
    text(dist(i), -5, num2str(stats(i)), 'vertical', 'bot', 'horizontal', 'center');
end

% Add topography 
load ~/data/GSL_transect.mat
x_trans = x_trans-212;

x = [x_trans x_trans(end) x_trans(1) x_trans(1)];
y = [H_trans max(zbin) max(zbin) H_trans(1)];

patch(x, y, [1 1 1]*.0, 'edgecolor', 'none');
xlim([-25 700])

LSLE = [-20 190];
GSL = [190 800];
Yvert = 415;
plot(LSLE, [Yvert Yvert], '-', 'color', [1 1 1], 'linewidth', 2)
plot(GSL, [Yvert Yvert], '-', 'color', [1 1 1], 'linewidth', 2)
text(190, Yvert, '|', 'color', [1 1 1], 'fontSize', 15, 'verticalAlignment', ...
     'middle', 'fontWeight', 'bold')
text(-20, Yvert, '|', 'color', [1 1 1], 'fontSize', 15, 'verticalAlignment', ...
     'middle', 'horizontalAlignment', 'right', 'fontWeight', 'bold')
text(75, Yvert, 'LSLE', 'color', [1 1 1], 'fontSize', 12, 'verticalAlignment', ...
     'bottom', 'horizontalAlignment', 'center')
text(450, Yvert, 'GSL', 'color', [1 1 1], 'fontSize', 12, 'verticalAlignment', ...
     'bottom', 'horizontalAlignment', 'center')
% $$$ cb = colorbar;
% $$$ ti = ylabel(cb, '[NO_x] (mmol m^{-3})');

% $$$ cbPos = get(cb, 'pos');
% $$$ cbPos(1) = .01+cbPos(1);
% $$$ cbPos(2) = cbPos(2)+.1;
% $$$ cbPos(3) = .75*cbPos(3);
% $$$ cbPos(4) = cbPos(4)-.2;
% $$$ set(cb, 'pos', cbPos)

adjust_space
pos1 = get(gca, 'pos');
vertExt = pos1(2)+pos1(4);
adjust_space
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = vertExt - pos2(2);
set(gca, 'pos', pos2)


subplot(414)

plot(xbin, nanmean(NO3field, 1), 'k', 'linewidth', 2)

xlabel('Distance from Stat. 25 (km)');
ylabel('[NO_x] (mmol m^{-3})')
xlim([-25 700])
ylim([15.5 19])
set(gca, 'tickdir', 'out')

adjust_space

%save figure
print('-dpng', 'NOxTransect_withplot.png')
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'NOxTransect_withplot.eps')
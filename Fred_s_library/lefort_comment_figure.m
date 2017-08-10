clear

load H_of_x

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.01 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

subplotheight = 5.2; %cm
totalfigheight = subplotheight./figh; %cm
totalfigwidth = 20; %cm

% --- customizing figure --- %
CONTOUR_SPACE = 10000;
CONTOUR_FONT = 8;
CONTOUR_NOLINES = 200;
V2 = [60:10:90];
V1 = [100 125 150];


figure(1)
clf
%set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 40 15])
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 totalfigwidth totalfigheight])
%imagesc(x_c, sigma, c_tp1);

load cmap_modif.mat
colormap(cmap)

% DANIEL FORTRAN
% Load la solution
load O2_FORTRAN_1km.mat

% Niveau de référence du rigid lid.
z0 = 150;

xx = (chi(1)-chi(2));

subplot(3,1,1)
contourf((chi-5*xx)/1000,z+z0,O', CONTOUR_NOLINES, 'linestyle', 'none');
hold on
[cs1, h1] = contour(chi/1000,z+z0,O', V1, 'k');
[cs2, h2] = contour(chi/1000,z+z0,O', V2, 'color', [1 1 1]);

clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, ...
       'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
axis([0 824 150 525])
xlim([0 825])
ylim([150 525])
caxis([0 350])
set(gca, 'xticklabel', []);
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
adjust_space


% DANIEL MATLAB
% Load la solution
load O2_dx_1km.mat

% Niveau de référence du rigid lid.
z0 = 150;

subplot(3,1,2)
contourf(chi/1000,z+z0,O', CONTOUR_NOLINES, 'linestyle', 'none')
hold on
[cs1, h1] = contour(chi/1000,z+z0,O', V1, 'k');
[cs2, h2] = contour(chi/1000,z+z0,O', V2, 'color',[1 1 1]);
clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, ...
       'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
axis([0 824 150 525])
xlim([0 825])
ylim([150 525])
caxis([0 350])
ylab = ylabel('Depth (m)');
set(gca, 'xticklabel', []);
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
cb = colorbar;
adjust_space

% ylabel position
ylabPos = get(ylab, 'pos');
ylabPos(2) = 140;
ylabPos(1) = -68;
set(ylab, 'pos', ylabPos);

set(cb,'YTick',[0:50:350])


% FRED MATLAB
%load O2_fred_Kx1.mat
%load O2_Fred_1km_20121016.mat
%load O2_1km_2012-10-09.mat
%load O2_final_slow.mat
load O2_final_20121016.mat

subplot(3,1,3)
contourf(X/1000,Z+150,c_t, CONTOUR_NOLINES, 'linestyle', 'none')
hold on
[cs1, h1] = contour(X/1000,Z+150,c_t, V1, 'k');
[cs2, h2] = contour(X/1000,Z+150,c_t, V2, 'color',[1 1 1]);
clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, ...
       'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
axis([0 824 150 525])
xlim([0 825])
ylim([150 525])
caxis([0 350])
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
xlabel('Seaward distance from STN 25 (km)')
adjust_space



set(gcf, 'renderer', 'painters')
print('-dpng','comment_fig.png')
print('-depsc','comment_fig.eps')


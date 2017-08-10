clear

load H_of_x


% DANIEL FORTRAN
% Load la solution
load O2_FORTRAN_1km.mat

% Lefort solut. (interpolate)
load Oxygen_Output.mat
xLefort = 1:size(Ox_ref,2);
zLefort = 1:size(Ox_ref,1);
zLefort = fliplr(zLefort)+150;


% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
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


% Niveau de référence du rigid lid.
z0 = 150;

xx = (chi(1)-chi(2));

contourf((chi-5*xx)/1000,z+z0,O', CONTOUR_NOLINES, 'linestyle', 'none');
hold on
[cs1, h1] = contour(chi/1000,z+z0,O', V1, 'k');
[cs2, h2] = contour(chi/1000,z+z0,O', V2, 'color', [1 1 1]);

clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, ...
       'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
xlim([0 825])
ylim([150 525])
caxis([0 350])
ylab = ylabel('Depth (m)', 'fontWeight', 'bold');
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
cb = colorbar;

set(cb,'YTick',[0:50:350])

xlabel('Seaward distance from STN 25 (km)')
adjust_space

% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-dpng','comment_fig.png')
% $$$ print('-depsc','comment_fig.eps')





% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.01 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
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


figure(2)
clf
%set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 40 15])
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 totalfigwidth totalfigheight])
%imagesc(x_c, sigma, c_tp1);

%load cmap_modif.mat
load cmap_modif2.mat
colormap(cmap)

% DANIEL FORTRAN
% Load la solution
load O2_FORTRAN_1km.mat


subplot(211)
contourf(xLefort, zLefort, Ox_ref, CONTOUR_NOLINES, 'linestyle', 'none');
hold on
[cs1, h1] = contour(xLefort, zLefort, Ox_ref, V1, 'k');
[cs2, h2] = contour(xLefort, zLefort, Ox_ref, V2, 'color', [1 1 1]);

clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, 'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
xlim([0 825])
ylim([150 525])
caxis([0 350])
ylab = ylabel('Depth (m)', 'fontSize', 11,'fontWeight', 'bold');
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
set(gca, 'xticklabel', [])
set(gca, 'fontSize', 11)
cb = colorbar;
set(cb,'YTick',[0:50:350])
adjust_space
pause(1)
colorbar('off')
text(10, 510, 'a) Lefort et al. [2012]', 'verticalAlignment', 'bottom', ...
     'horizontalAlignment', 'left', 'fontSize', 12, 'fontWeight', ...
     'bold', 'color', [1 1 1])

subplot(212)
% Niveau de référence du rigid lid.
z0 = 150;

xx = (chi(1)-chi(2));
contourf((chi-5*xx)/1000,z+z0,O', CONTOUR_NOLINES, 'linestyle', 'none');
hold on
[cs1, h1] = contour(chi/1000,z+z0,O', V1, 'k');
[cs2, h2] = contour(chi/1000,z+z0,O', V2, 'color', [1 1 1]);

clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, ...
       'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
xlim([0 825])
ylim([150 525])
caxis([0 350])
ylab = ylabel('Depth (m)', 'fontSize', 11, 'fontWeight', 'bold');
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
cb = colorbar;
set(cb,'YTick',[0:50:350])
adjust_space
pause(1)
cb_pos = get(cb, 'pos');
cb_pos(1) = cb_pos(1)-.02;
cb_pos(3) = cb_pos(3)*.8;
cb_pos(2) = cb_pos(2)+.05;
cb_pos(4) = cb_pos(4)*2-.1;
set(cb,'pos',cb_pos)
set(gca, 'fontSize', 11)
xlabel('Seaward distance from STN 25 (km)', 'fontSize', 11, 'fontWeight', 'bold')
text(10, 500, 'b) FORTRAN code adapted from Bourgault and Kelley [2004]', 'verticalAlignment', 'bottom', ...
     'horizontalAlignment', 'left', 'fontSize', 12, 'fontWeight', 'bold', 'color', [1 1 1])


set(gcf, 'renderer', 'painters')
print('-dpng','comment_fig2.png')
print('-depsc','comment_fig2.eps')

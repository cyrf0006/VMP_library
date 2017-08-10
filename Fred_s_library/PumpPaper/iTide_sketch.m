clear all

% --- Get topography  --- %
load ~/data/GSL_transect.mat
x_trans = x_trans-190;
x = [x_trans x_trans(end) x_trans(1) x_trans(1)];
y = [H_trans 200 200 H_trans(1)];

% --- Get No3 field --- %
load transect/transect_info.mat
xbin = (xbin+10);

ox = xbin;
oz = zbin(2:end);

% Extrapolate data
for i = 1:length(xbin)
    I = find(~isnan(NO3field(:,i)));
    NO3field(:,i) = interp1(zbin(I), NO3field(I,i), zbin, 'nearest', 'extrap');
end
NO3field = [NO3field NO3field(:,end)];
xbin = [xbin -30];

% --- Build displacements matrix --- %
%eta = load('~/PhD/Nitrates/eta_struct.mat');
zz = [7.5 25:25:150];
xeta = -100:200;
% fake structure
eta0 = [-10 0 10 25 35 25 10];
etaMat = nan(length(zz), length(xeta));
wavelen = 60; %km
for i = 1:length(zz)
% $$$     [Y, I] = min(abs(zz(i)-eta.z));
% $$$     eta0 = eta.eta(I)*3;
    etaMat1(i, :) = eta0(i).*cos(2*pi*xeta./wavelen)+zz(i);
    etaMat2(i, :) = eta0(i).*cos(2*pi*xeta./wavelen+pi)+zz(i);
end
    



figure(2)
clf
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.5 ; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.07; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

XMAX = 125;
XMIN = -30;
ZMAX = 150;
ZMIN = 0;

load mycm
colormap(mycm)

subplot(211)
V=[0:.1:25];
contourf(xbin, zbin, NO3field, V, 'linestyle', 'none')
caxis([3 21])
hold on
plot([ox(end) ox(end) XMAX-eps  XMAX-eps ox(end)], [ZMAX-eps, oz(1), oz(1), ...
                    ZMAX-eps, ZMAX-eps], '--k', 'linewidth', 2)
plot(xeta, etaMat2, 'color', [.3 .3 .3], 'linewidth', 0.25)
patch(x, y, [1 1 1]*.7, 'edgecolor', 'none');
hold off
%xlim([0 125])
xlim([-30 XMAX])    
ylim([0 150])    
set(gca, 'ydir', 'reverse')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'ytick', [0:25:150])
set(gca, 'ygrid', 'on')
ylabel('Depth (m)')
set(gca, 'xticklabel', [])
text(122, 140, 'a', 'fontsize', 12, 'fontWeight', 'bold','BackgroundColor',[1 1 1])
adjust_space


subplot(212)
contourf(xbin, zbin, NO3field, V, 'linestyle', 'none')
caxis([3 21])
hold on
plot([ox(end) ox(end) XMAX-eps  XMAX-eps ox(end)], [ZMAX-eps, oz(1), oz(1), ...
                    ZMAX-eps, ZMAX-eps], '--k', 'linewidth', 2)
plot(xeta, etaMat1, 'color', [.3 .3 .3], 'linewidth', 0.25)
patch(x, y, [1 1 1]*.7, 'edgecolor', 'none');
hold off
%xlim([0 125])
xlim([-30 XMAX])    
ylim([0 150])    
set(gca, 'ydir', 'reverse')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'ytick', [0:25:150])
set(gca, 'ygrid', 'on')
xlabel('Distance from sill (km)')
text(122, 140, 'b', 'fontsize', 12, 'fontWeight', 'bold','BackgroundColor',[1 1 1])
text(-29, 137.5, 'Head of the', 'verticalAlignment', 'bottom', ...
     'horizontalAlignment', 'left')
text(-29, 137.5, 'Laurentian Channel', 'verticalAlignment', 'top', ...
     'horizontalAlignment', 'left')

adjust_space
print('-dpng', '-r300', 'sketch_topo2.png')
set(gcf, 'renderer', 'painters')
print('-depsc2', 'sketch_topo2.eps')





% $$$ figure(1)
% $$$ clf
% $$$ set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 20 15])
% $$$ patch(x, y, [1 1 1]*.7, 'edgecolor', 'none');
% $$$ %xlim([0 125])
% $$$ xlim([-30 125])    
% $$$ ylim([0 150])    
% $$$ set(gca, 'ydir', 'reverse')
% $$$ set(gca, 'box', 'on')
% $$$ set(gca, 'tickdir', 'out')
% $$$ set(gca, 'ytick', [0:25:150])
% $$$ set(gca, 'ygrid', 'on')
% $$$ print('-dpng', '-r300', 'sketch_topo.png')
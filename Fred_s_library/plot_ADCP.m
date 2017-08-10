function plot_ADCP(fname)

% function plot_ADCP(fname)
%
% plot couple of fields from the ADCP. fname should be the .mat
% ADCP file, i.e., after rdradcp.m have been performed. 


% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.07; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

umin = -1.0;
umax =  1.0;
wmin = -0.15;
wmax =  0.1;

load(fname);


% Set the z axis
z0 = 0.8; % Depth of the transducer. This info can come from the ADCP.depth data
z = ADCP.config.ranges + z0;

% Set the time axis
mtime = ADCP.mtime;
dvec = datevec(mtime(1));
mtime0 = datenum(dvec(1),dvec(2),dvec(3),0,0,0);


% Velocity
U = ADCP.east_vel;
V = ADCP.north_vel;
W = ADCP.vert_vel;

% Backscatter intensity
I1 = ADCP.intens(:,1,:);
I2 = ADCP.intens(:,2,:);
I3 = ADCP.intens(:,3,:);
I4 = ADCP.intens(:,4,:);
 
I1 = squeeze(I1);
I2 = squeeze(I2);
I3 = squeeze(I3);
I4 = squeeze(I4);
Imean = (I1 + I2 + I3 + I4)/4;

% figure(1);
% pcolor(mtime-mtime0,z,I1);
% shading('interp');
% set(gca,'ydir','reverse');
% set(gca,'XTickLabel',[]);
% ylabel('Depth (m)');
% datetick('x',15);
% xlabel(['Time on ',datestr(mtime0)]);
% colorbar;
% 
% figure(2);
% pcolor(mtime-mtime0,z,I2);
% shading('interp');
% set(gca,'ydir','reverse');
% set(gca,'XTickLabel',[]);
% ylabel('Depth (m)');
% datetick('x',15);
% xlabel(['Time on ',datestr(mtime0)]);
% colorbar;
% 
% figure(3);
% pcolor(mtime-mtime0,z,I3);
% shading('interp');
% set(gca,'ydir','reverse');
% set(gca,'XTickLabel',[]);
% ylabel('Depth (m)');
% datetick('x',15);
% xlabel(['Time on ',datestr(mtime0)]);
% colorbar;
% 
% figure(4);
% pcolor(mtime-mtime0,z,I4);
% shading('interp');
% set(gca,'ydir','reverse');
% set(gca,'XTickLabel',[]);
% ylabel('Depth (m)');
% datetick('x',15);
% xlabel(['Time on ',datestr(mtime0)]);
% colorbar;



figure(5)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 20 14])


set(gcf, 'renderer', 'painters')
subplot(3,1,1);
pcolor(mtime-mtime0,z,U);
shading('interp');
caxis([umin umax]);
set(gca,'ydir','reverse');
%set(gca,'XTickLabel',[]);
ylabel('Depth (m)', 'FontSize', 10);
datetick('x',15);
%xlabel(['Time on ',datestr(mtime0)], 'FontSize', 10);
title('Meridional velocity (m/s)', 'FontSize', 10)
set(gca, 'FontSize', 10)
c = colorbar('FontSize', 10);
ylabel(c,'Velocity (m/s)', 'FontSize', 10)
set(gca, 'xticklabel', [])

adjust_space

subplot(3,1,2);
pcolor(mtime-mtime0,z,V);
shading('interp');
caxis([umin umax]);
set(gca,'ydir','reverse');
%set(gca,'XTickLabel',[]);
ylabel('Depth (m)', 'FontSize', 10);
datetick('x',15);
%xlabel(['Time on ',datestr(mtime0)], 'FontSize', 10);
title('Zonal velocity (m/s)', 'FontSize', 10)
set(gca, 'FontSize', 10)
c = colorbar('FontSize', 10);
ylabel(c,'Velocity (m/s)', 'FontSize', 10)
set(gca, 'xticklabel', [])

adjust_space

subplot(3,1,3);
pcolor(mtime-mtime0,z,W);
shading('interp');
caxis([wmin wmax]);
set(gca,'ydir','reverse');
ylabel('Depth (m)', 'FontSize', 10);
datetick('x',15);
xlabel(['Time on ',datestr(mtime0)], 'FontSize', 10);
title('Vertical velocity (m/s)', 'FontSize', 10)
set(gca, 'FontSize', 10)
c = colorbar('FontSize', 10);
ylabel(c,'Velocity (m/s)', 'FontSize', 10)

adjust_space

% Save figure
print('-dpng', '-r300','uvw_adcp.png')


figure(8);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 20 14])

set(gcf, 'renderer', 'painters')


pcolor(mtime-mtime0,z,Imean);
shading('interp');
set(gca,'ydir','reverse');
set(gca,'XTickLabel',[]);
ylabel('Depth (m)', 'FontSize', 10);
datetick('x',15);
xlabel(['Time on ',datestr(mtime0)], 'FontSize', 10);
title('Mean backscatter intensity', 'FontSize', 10)
set(gca, 'FontSize', 10)
colorbar;
c = colorbar('FontSize', 10)   ;
ylabel(c,'Backscatter', 'FontSize', 10)

% Save figure
print('-dpng', '-r300','intens_adcp.png')

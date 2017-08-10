clear
% $$$ % *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.05; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.05; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

% --- Customizable info -- %
paperwidth = 26;%cm
paperheight = 20;%cm

SEABED_COLOR = [.7 .7 .7];
LAND_COLOR = [0.9529 0.8706 0.7333];
WATER_COLOR = [0.0600 0.4820 0.6750];



% --- Get bathym
load bathySTLE.mat % Thx to Simon Senn.

Z = myBathym(1:65, 7:end-3);
x = [1:size(Z,2)]*400/1000;
y = [1:size(Z,1)]*400/1000;
x = fliplr(x);

% get coastline
X = X(1:65, 7:end-3); % to draw coastline
Y = Y(1:65, 7:end-3);
myY = X(1:end,1); % ok if swap
myX = Y(1,1:end);
X = myX;
Y = myY;

load stle.xpyp; % Thx to Simon Senn.
id = find(stle(:,1)==99999.99);
stle(id,1)=NaN;
stle(id,2)=NaN;

I = find((stle(:,1)>=Y(1) & stle(:,1)<=Y(end)) | isnan(stle(:,1))==1);
stle = stle(I,:);
I = find((stle(:,2)>=X(1) & stle(:,2)<=X(end)) | isnan(stle(:,2))==1);
stle = stle(I,:);
stle(:,1) = (stle(:,1)-Y(1))*.4+.4;
stle(:,2) = abs(stle(:,2)-X(end))*.4+.4;


% top mask
% $$$ maskTop = ones(size(Z));
% $$$ I = find(isnan(Z)==1);
% $$$ maskTop(I) = 0;
% $$$ for i = 1:length(x)
% $$$     %I = find(mask)
% $$$ end
xstle = stle(:,2);
ystle = stle(:,1);
I = find(isnan(xstle) == 1);
xstle(I) = [];
ystle(I) = [];
% This is tweak!
xstle(end-6:end) = [];
ystle(end-6:end) = [];


% side-x mask
zMax = Z(end,:);
zpatchx = [zMax 350 350 zMax(1)];
xpatch = [x x(end) x(1) x(1)];
zpatchx_blue = [zMax 0 0 zMax(1)];

% side-y mask
zMax = Z(:,end);
I = find(isnan(zMax)==1);
zMax(I) = 0;
zpatchy = [zMax' 350 350 zMax(1)];
zpatchy_blue = [zMax' 0 0 zMax(1)];
ypatch = [y y(end) y(1) y(1)];

% $$$ % Side temperature % Thx STOS (failed attemp)
% $$$ load '~/PhD/CTD_MPO/transect_O2/transect_lin_interp.mat'
% $$$ I  = find(xbin<100);
% $$$ xbin = xbin(I);
% $$$ T = T_transect(:,I);
% $$$ % $$$ pcol3=pcolor(T), shading('flat');
% $$$ % $$$ load BR_symetric
% $$$ % $$$ % $$$ contourf(T, 100, 'linestyle', 'none')
% $$$ % $$$ % $$$ hold on
% $$$ % $$$ caxis([0 5])
% $$$ [c, h] = contour(xbin, zbin, T, [1.5 1.5]);
% $$$ 
% $$$ % I = find(diff(c(1,:))<0); 
% $$$ % This is also full tweak!!
% $$$ xVec = [c(1,2:121) c(1,121) fliplr(c(1,123:end)) c(1,123)];
% $$$ xVec1 = max(xVec)-xVec;
% $$$ %zVec = [c(1,2:121) c(1,121) fliplr(c(1,123:end)) c(1,123)]; [c(1,2:121) c(1,121) fliplr(c(1,123:end)) c(1,123)];
% $$$ yVec1 = [c(2,2:121) c(2,121) fliplr(c(2,123:end)) c(2,1)];    
% $$$ 
% $$$ xVec2 = [26 15 10 5 0 0 5 15 26 26];
% $$$ yVec2 = [111.5998 113 114 118 120 46 48 52 53.5525 111.5998];

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])

%surf(x, y, Z, 'linestyle', 'none')
%set(gca, 'zdir', 'reverse')
%shading interp

% Add masks
%hold on
%pcolor(x, y, maskTop)
%shading interp
%caxis([0 1])

fill3(xpatch, xpatch*0+y(end), zpatchx_blue, WATER_COLOR)
hold on
fill3(ypatch*0+x(end), ypatch, zpatchy_blue, WATER_COLOR)

% CIL attempt
% $$$ fill3(xVec1, xVec1*0+y(end), yVec1, WATER_COLOR*1.2)
% $$$ fill3(xVec2*0+x(end), xVec2, yVec2, WATER_COLOR*1.2)

fill3(xpatch, xpatch*0+y(end), zpatchx, SEABED_COLOR)
set(gca, 'zdir', 'reverse')
%hold on
fill3(xpatch, xpatch*0+y(end), zpatchx, SEABED_COLOR)
fill3(ypatch*0+x(end), ypatch, zpatchy, SEABED_COLOR)

% $$$ fill3(xpatch, xpatch*0+y(end), zpatchx_blue, WATER_COLOR)
% $$$ fill3(ypatch*0+x(end), ypatch, zpatchy_blue, WATER_COLOR)

% Top mask
xpatch = [xstle; xstle(end); xstle(1); xstle(1)];
ypatch = [ystle; y(1); y(1); ystle(1)];
patch(xpatch, ypatch, LAND_COLOR);

xpatch = [xstle; xstle(end); xstle(1); xstle(1)];
ypatch = [ystle; y(end); y(end); ystle(1)];
patch(xpatch, ypatch, WATER_COLOR);

ylim([y(1) y(end)])
xlim([x(end) x(1)])
zlim([0 350])

% Customize axis
%set(gca, 'xtick', [10:90])
set(gca, 'xtickLabel', [90:-10:10])
%set(gca, 'ytick', [5:5:25])
%pause(.1)
set(gca, 'ytickLabel', [35:-5:15])


% PLot box display
set(gca, 'box', 'on')
set(gca, 'PlotBoxAspectRatio', [4 2 1])
set(gca, 'DataAspectRatio', [0.9492 0.5000 13.6719])
set(gca, 'CameraViewAngle', [10])
set(gca, 'CameraPosition', [-386.013 189.031 -1217.6])
set(gca, 'CameraTarget', [49 13.2 175])
set(gca, 'CameraUpVector', [0 0 -1])
adjust_space


% $$$ zlabel('Depth (m)')
% $$$ ylabel('Dist from SouthShore (km)')
% $$$ xlabel('Dist from Tadoussac (km)')

keyboard

print(gcf, '-dpng', '-r300', 'recapFigure.png')
set(gcf, 'renderer', 'painters')
print(gcf, '-depsc2', 'recapFigure.eps')

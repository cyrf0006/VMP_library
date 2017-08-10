function rotateLSLE(field, middleTransect, latLims, lonLims)
    

    
% usage ex in /home/cyrf0006/PhD/Nitrates/SST:
% rotateLSLE('SST_averaged.mat', 'middleTransect.dat', [46.8 50], [-71 -66.8])

    
    
data = load(field);
names = fieldnames(data);
eval(['x = data.' names{1} ';']);
eval(['y = data.' names{2} ';']);
eval(['z = data.' names{3} ';']);


transect = load(middleTransect);
transLon = transect(1:end-1,1);
transLat = transect(1:end-1,2);

dlat = .05;
dlon = .05;
transLonItp = transLon(1):dlon:transLon(end);
transLatItp = interp1(transLon, transLat, transLonItp, 'cubic');

% Remove pts not in LSLE
z = flag_SST(z, x, y, [46.8 50],[-71 -66.8]);

    
% loop on transect pts
% rename tranect since it is long!
xt = transLonItp;
yt = transLatItp;
slat = nan(size(xt));
slon = nan(size(xt));
nlat = nan(size(xt));
nlon = nan(size(xt));
d = 35000; %m
for i = 2:length(transLonItp)-1
    X1 = xt(i-1);
    X2 = xt(i+1);
    Y1 = yt(i-1);
    Y2 = yt(i+1);
    
    ang = atand((Y2-Y1)/(X2-X1));
    
    
    [slon(i), slat(i), a21] = m_fdist(xt(i), yt(i), 180-ang, d);
    [nlon(i), nlat(i), a21] = m_fdist(xt(i), yt(i), 360-ang, d);
end
    
nlon = nlon-360;
slon = slon-360;

 % Find pixel along the track
[X,Y] = meshgrid(x,y); 
T = [];
xT = [];
W = [];
Tint = [];
for i = 2:length(nlat)-1    
    
    % Across shore temperature transect
    [range, lons, lats] = m_lldist([nlon(i) slon(i)], [nlat(i) slat(i)], 80);    
    Titp = interp2(X,Y,z,lons,lats);
    T = [T nanmean(Titp)]; 
    Tint = [Tint nansum(Titp)]; 

    % distance from Ile d'Orlean
    range = m_lldist([transLonItp(1) transLonItp(i)], [transLatItp(1) transLatItp(i)], 80);  
    xT = [xT range];
    
    % channel width
    I = find(~isnan(Titp)==1);
    W = [W length(I)];
end



% PLots!
    
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 20])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.01; % top of figure
bots = 0.02; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

m_proj('mercator', 'long',[lonLims(1) lonLims(2)],'lat',[latLims(1) latLims(2)]);
m_grid('box','fancy')
hold on
%m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
m_pcolor(x,y,z)     
shading flat
m_gshhs_h('patch',[.5 .5 .5]); %coastlines
m_plot(transLonItp, transLatItp, 'k', 'linewidth', 2)
m_plot(transLonItp, transLatItp, 'color', [1 1 1], 'linestyle','.', 'linewidth', 2)
hold off
colorbar

% add transects
hold on
for i = 1:length(nlat)
    m_plot([nlon(i) slon(i)], [nlat(i) slat(i)], 'color', [1 1 1])
end
hold off

% add km marks
marks = 50:50:350;
for i = 1:length(marks)
    [Y, I] = min(abs(xT- marks(i)));
    m_text(nlon(I), nlat(I), num2str(marks(i)),'vertical', 'bottom', ...
           'horizontal', 'right', 'color', [1 1 1],'FontSize',10)
end

adjust_space






print(gcf, '-dpng', 'rotate_LSLE_map.png')




figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 20])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


subplot(311)
plot(xT, T)
ylabel('T ^{\circ}C')
set(gca, 'xticklabel', [])
adjust_space

subplot(312)       
plot(xT, W*1.1*1.1)
ylabel('pixel area (km^2)')            
set(gca, 'xticklabel', [])
adjust_space

subplot(313)
plot(xT, Tint./(W*1100*1100));
ylabel('H (^{\circ}C m^{-2})')
xlabel('dist. from Ile d''Orleans (km)')
adjust_space

print(gcf, '-dpng', 'HeatFlux_subplot.png')

keyboard
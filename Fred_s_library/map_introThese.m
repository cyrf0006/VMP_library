clear
% to be run in ~/WINDEX/data_processing/Mission_tadoussac_2009/

% *********************** Adjust_space.m ************************ %
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
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
% ------------- LSLE MAP ---------------- %

%load '/home/cyrf0006/data/Topo_Gulf_St_Lawrence.mat'
load '/home/cyrf0006/data/SHC/500m/CHS_500m_gulf.mat'


%figure dimension
paperwidth = 22;%cm
paperheight = 18;%cm

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; % offset between heigth of colorbar and heigth of figure

lon_min=-70;
lon_max=-67;
lat_min=48;
lat_max=49.5;

%Points for transect (x,y = lon,lat)
%A = [-68-50/60 48+52/60];
%B = [-68-25/60 48+32/60];
A = [-68.422966 48.537659];
B = [-68.837865 48.850145];

I=find(lat<lat_max & lat>lat_min);
latitude=lat(I);
longitude=lon(I);
bathy=z(I);

I=find(longitude<lon_max & longitude>lon_min);
latitude2=latitude(I);
longitude2=longitude(I);
bathy2=bathy(I);

clear latitude longitude I bathy lat lon z

% on the new grid
y = lat_min:(lat_max-lat_min)/110/5:lat_max;
x = lon_min:(lon_max-lon_min)/80/5:lon_max;
[X,Y] = meshgrid(x,y);
[XI,YI,Z] = griddata(longitude2,latitude2,bathy2,X,Y);

%b=gaussfir(0.001, 5, 3);  %Gaussian filter designer ()
b=gaussfir(10);  %Gaussian filter designer ()

Z2=filter2(b,Z);  %Application of the filter

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])


m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
V=[0:10:300];
caxis([V(1) V(end)]) ;

% -- custom colormap! -- %
load ~/PhD/bathym/gebco64
colormap(gebco);

%pack
hold on
[H, H] = m_contourf(X,Y,Z2,V, 'linestyle', 'none');
[HH, HH] = m_contour(X,Y,Z2, [20 120 200 300], 'color', 'k');


xlabel('Longitude', 'FontSize', 10)
ylabel('Latitude', 'FontSize', 10)
m_gshhs_h('patch',[.7 .7 .7]); %coastlines                               
%m_coast('patch',[.5 .5 .5]); %coastlines

m_line(-69.716606, 48.156604,'marker','o','MarkerFaceColor','k','markersize',6,'color','k');
m_text(-69.8, 48.156604, 'Tadoussac', 'vertical', 'bottom', 'horizontal', 'center', 'color', 'k','FontSize',10, 'FontWeight', 'bold');
m_line(-68.522186, 48.482571, 'marker','o','MarkerFaceColor','k','markersize',6,'color','k');
m_text(-68.522186, 48.482571, 'Rimouski', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',10, 'FontWeight', 'bold');
m_line(-67.386367, 49.320459, 'marker','o','MarkerFaceColor','k','markersize',6,'color','k');
m_text(-67.386367, 49.320459, 'Pointe-des-Monts', 'vertical', 'bottom', 'horizontal', 'right', 'color', 'k','FontSize',10, 'FontWeight', 'bold');
m_line(-69.55, 48.15, 'marker','o','markersize',30, 'color', 'k', 'linewi', 2, 'LineStyle', '--');
m_line([-69.55 -69],[48.15 48.15],'color','k','linewi',2, 'linestyle', '--');    
m_text(-69, 48.15, 'tete du', 'vertical', 'bottom', 'horizontal', 'left', 'color', 'k','FontSize',10, 'FontWeight', 'bold');
m_text(-69, 48.15, 'chenal Laurentien', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',10, 'FontWeight', 'bold');

set(gca, 'fontsize', 10)


adjust_space

m_grid('box','fancy')

% colorbar
cbar_offset = 0.04; 
offset2 = 0.02; % offset between heigth of colorbar 
c = colorbar('location','south');
cb_pos = get(c, 'position');
cb_back = cb_pos;
cb_pos(4) = cb_pos(4)*.3;
cb_pos(1) = .5;
cb_pos(3) = cb_pos(3)*.4;
cb_pos(2) = cb_pos(2)+cbar_offset;
set(c, 'pos', cb_pos)
set(c, 'fontsize', 9)

ti = ylabel(c,'Profondeur (m)', 'FontSize', 9);
ti_pos = get(ti, 'position');
set(ti, 'Rotation',0.0);     
clim = get(gca, 'clim');
set(ti, 'position', [(clim(2)-clim(1))./2  4 ti_pos(3)]); 






% ------------------ Gulf Map ------------------------ %

load '/home/cyrf0006/data/Topo_Gulf_St_Lawrence.mat'

LSLE_lon = [-69.341126, -67.006531];
LSLE_lat = [47.916342, 48.951366];

lon_min=-71;
lon_max=-56;
lat_min=45;
lat_max=52;

I=find(lat<lat_max & lat>lat_min);
latitude=lat(I);
longitude=lon(I);
bathy=z(I);

I=find(longitude<lon_max & longitude>lon_min);
latitude2=latitude(I);
longitude2=longitude(I);
bathy2=bathy(I);

clear latitude longitude I bathy lat lon z

% on the new grid
y = lat_min:(lat_max-lat_min)/110/5:lat_max;
x = lon_min:(lon_max-lon_min)/80/5:lon_max;
[X,Y] = meshgrid(x,y);
[XI,YI,Z] = griddata(longitude2,latitude2,bathy2,X,Y);

%b=gaussfir(0.001, 5, 3);  %Gaussian filter designer ()
b=gaussfir(10); 
Z2=filter2(b,Z);  %Application of the filter


%a2 = axes('position',[0.15 0.35 0.28 0.80]) ; % inset
a2 = axes('position',[0.08 0.36 0.35 0.82]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

V=[0:25:300];
caxis([V(1) V(end)]) ;

%pack
hold on
[H, H] = m_contourf(X,Y,Z2,V); %filtered

%[H, H] = m_contourf(X,Y,Z,V); %not filtered
set(H, 'LineStyle', 'none')

m_gshhs_h('patch',[.7 .7 .7]); %coastlines
m_line([-56.43 -59],[51.55 51.55],'color','k','linewi',1, 'linestyle', '-');    
m_text(-59, 51.55, 'detroit de', 'vertical', 'bottom', 'horizontal', 'right', 'color', 'k','FontSize',8, 'FontWeight', 'bold')             
m_text(-59, 51.55, 'Belle-Isle', 'vertical', 'top', 'horizontal',  'right', 'color', 'k','FontSize',8, 'FontWeight', 'bold') 
% $$$ m_text(-58.5, 50.9, 'detroit de', 'vertical', 'bottom', 'horizontal', 'left', 'color', 'k','FontSize',8, 'FontWeight', 'bold')             
% $$$ m_text(-58.7, 50.9, 'Belle-Isle', 'vertical', 'top', 'horizontal',  'left', 'color', 'k','FontSize',8, 'FontWeight', 'bold') 
% $$$ m_text(-57.5, 51.4, '\rightarrow', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',10, 'rotation', 35)
m_text(-60.5, 47.35, 'detroit', 'vertical', 'bottom', 'horizontal', 'left', 'color', 'k','FontSize',8, 'FontWeight', 'bold')
m_text(-60.5, 47.35, 'de Cabot', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',8, 'FontWeight', 'bold')
m_text(-59, 46.8, 'chenal', 'vertical', 'bottom', ...
       'horizontal', 'left', 'color', 'k','FontSize',8, 'rotation', -38, 'FontWeight', 'bold')
m_text(-59, 46.8, 'Laurentien', 'vertical', 'top', ...
       'horizontal', 'left', 'color', 'k','FontSize',8, 'rotation', -38, 'FontWeight', 'bold')
m_text(-63.364105, 48.25, 'golfe du', 'vertical', 'bottom', ...
       'horizontal', 'left', 'color', 'k','FontSize',8, 'rotation', 0, 'FontWeight', 'bold')
m_text(-63.364105, 48.25, 'Saint-Laurent', 'vertical', 'top', ...
       'horizontal', 'left', 'color', 'k','FontSize',8, 'rotation', 0, 'FontWeight', 'bold')
m_text(-68.522186, 48.482571, 'EMSL', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',8, 'FontWeight', 'bold');
m_text(-64.17, 49.3, 'det. d''Honguedo', 'vertical', 'middle', ...
       'horizontal', 'center', 'color', 'k','FontSize',8, 'rotation', -30, 'FontWeight', 'bold')


% LSLE restangle
lon_min=-70;
lon_max=-67;
lat_min=48;
lat_max=50;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', ...
          1, 'edgecolor', 'k', 'linestyle', '--') 

%includeStatsGulf


set(gca, 'fontsize', 8)


%save figure
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'map_introThese.eps')


keyboard



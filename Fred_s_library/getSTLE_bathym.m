% in /home/cyrf0006/PhD/RecapFigure
clear 

krom = nc_varget('STLE400m_GRID.nc','KROM');
irom = nc_varget('STLE400m_GRID.nc','IROM');
grid = nc_varget('STLE400m_GRID.nc','GRID');
J = nc_varget('STLE400m_GRID.nc','INDEXVALID2D');
lonrom = nc_varget('STLE400m_GRID.nc','LONGITUDES');
latrom = nc_varget('STLE400m_GRID.nc','LATITUDES');

[X, Y] = meshgrid(irom, krom);

bathym = nan(size(X));
bathym(J) = grid;

myBathym = bathym';
X = X';
Y = Y';
%lonrom = lonrom';
%latrom = latrom';
myBathym = myBathym(225:end, 700:end);
lonrom = lonrom(225:end, 700:end);
latrom = latrom(225:end, 700:end);
X = X(225:end, 700:end);
Y = Y(225:end, 700:end);


pcolor(lonrom, latrom, myBathym)
shading interp


save bathySTLE.mat myBathym lonrom latrom X Y


% $$$ % Find the corresponding coastline
% $$$ 
% $$$ latMin = min(min(latrom));
% $$$ lonMin = min(min(lonrom));
% $$$ latMax = max(max(latrom));
% $$$ lonMax = max(max(lonrom));
% $$$ 
% $$$ 
% $$$ m_proj('mercator','long',[lonMin lonMax],'lat',[latMin latMax]);
% $$$ m_gshhs_h('patch',[.7 .7 .7])
% $$$ m_gshhs('h','save','myCoast') 
% $$$ 
% $$$ [xx, yy] = m_ll2xy(lon, lat);

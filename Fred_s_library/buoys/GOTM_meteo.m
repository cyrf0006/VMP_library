% build meteo file (script not complete), used march 10



% 10m wind
meteo = load('meteo_raw.dat');
wind = load('wind_IML4_2008.dat');

wind2 = windlog_convert(wind(:,1), 2, 10);

uwind = wind2.*cosd(wind(:,2));
vwind = wind2.*sind(wind(:,2));

meteo(:,1) = uwind;
meteo(:,2) = vwind;

dlmwrite('meteo_w10m.dat', meteo,'delimiter',' ','precision',6)


% Gust wind speed
clear
meteo = load('meteo_raw.dat');
IML4 = load('IML4_meteo.dat');
wind = load('wind_IML4_2008.dat'); % needed for orientation

gust = IML4(:,5);
gust2 = windlog_convert(gust, 2, 10);

uwind = gust2.*cosd(wind(:,2));
vwind = gust2.*sind(wind(:,2));

meteo(:,1) = uwind;
meteo(:,2) = vwind;

dlmwrite('meteo_gust10m.dat', meteo,'delimiter',' ','precision',6)




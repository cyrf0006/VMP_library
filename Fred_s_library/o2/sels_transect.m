% script sels_transect.m
%
% See ecnlosed function below
clear 


[NO3field, xbin, zbin, dist] = getTransect('station_file_NO3');
[NO2field] = getTransect('station_file_NO2');
[PO4field] = getTransect('station_file_PO4');


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 25])

subplot(311)
V=[0:.01:.2];
contourf(xbin, zbin, NO2field, V, 'linestyle', 'none')
hold on
%contour(xbin, zbin, O_transect, [62.5 75 150 250 350], 'color', 'k', 'ShowText','on');
title('Nitrite concentration (mmol m^{-3})')
ylabel('depth (m)')
ylim([0 400])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
end
colorbar

subplot(312)
V=[0:.5:25];
contourf(xbin, zbin, NO3field, V, 'linestyle', 'none')
hold on
%contour(xbin, zbin, O_transect, [62.5 75 150 250 350], 'color', 'k', 'ShowText','on');
title('Nitrate + nitrite concentration (mmol m^{-3})')
ylabel('depth (m)')
ylim([0 400])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
end
colorbar

subplot(313)
V=[0:.05:2.5];
contourf(xbin, zbin, PO4field, V, 'linestyle', 'none')
hold on
%contour(xbin, zbin, O_transect, [62.5 75 150 250 350], 'color', 'k', 'ShowText','on');
title('Phosphate concentration (mmol m^{-3})')
ylabel('depth (m)')
ylim([0 400])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
end
xlabel('distance from Tadoussac (km)')
colorbar


figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 25])

subplot(311)
V=[0:.5:25];
contourf(xbin, zbin, NO3field, V, 'linestyle', 'none')
hold on
%contour(xbin, zbin, O_transect, [62.5 75 150 250 350], 'color', 'k', 'ShowText','on');
title('Nitrates  + nitrites concentration (mmol m^{-3})')
ylabel('depth (m)')
ylim([0 400])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
end
colorbar

subplot(312)
V=[0:.05:2.5];
contourf(xbin, zbin, PO4field, V, 'linestyle', 'none')
hold on
%contour(xbin, zbin, O_transect, [62.5 75 150 250 350], 'color', 'k', 'ShowText','on');
title('Phosphates concentration (mmol m^{-3})')
ylabel('depth (m)')
ylim([0 400])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
end
colorbar

subplot(313)
V=[0:.5:20];
contourf(xbin, zbin, NO3field./PO4field, V, 'linestyle', 'none')
hold on
contour(xbin, zbin, NO3field./PO4field, [16 16], 'color', 'k', 'ShowText','on');
title('N/P (mmol m^{-3})')
ylabel('depth (m)')
ylim([0 400])
set(gca, 'ydir', 'reverse')

for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 400], '--k')
end
xlabel('distance from Tadoussac (km)')
colorbar


%save transect_griddata_lin.mat xbin zbin T_transect O_transect



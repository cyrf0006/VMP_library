clear
% --- customizing figure --- %
CONTOUR_SPACE = 10000;
CONTOUR_FONT = 8;
CONTOUR_NOLINES = 200;
V2 = [40:10:90];
V1 = [100 125 150];



totalfigheight = 6;
totalfigwidth = 20; %cm

figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 totalfigwidth totalfigheight])
load O2.mat
load H_of_x
load cmap_modif.mat
colormap(cmap)

contourf(X/1000,Z+150,c_t, CONTOUR_NOLINES, 'linestyle', 'none')
hold on
[cs1, h1] = contour(X/1000,Z+150,c_t, V1, 'k');
[cs2, h2] = contour(X/1000,Z+150,c_t, V2, 'color',[1 1 1]);
clabel(cs1, h1, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE)
clabel(cs2, h2, 'fontsize', CONTOUR_FONT, 'labelspacing', CONTOUR_SPACE, ...
       'color', [1 1 1])
fill([x0/1000 x0(1)/10000], [H0 H0(end)], [40 32 33]/255)
hold off
axis([0 825 150 525])
xlim([0 825])
ylim([150 525])
caxis([0 350])
set(gca, 'ytick', [200 300 400 500]);
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out')
xlabel('Seaward distance from STN 25 (km)')
colorbar
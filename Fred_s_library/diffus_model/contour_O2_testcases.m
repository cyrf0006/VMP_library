clear

% Output idea
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 16 10])
O2min = 60;
O2max = 200;
ucurrent = 0.01; %m/s



% s1
s1 = subplot('position', [.1 .9 .72 .08]);
load SOD
SOD = SOD/10000*365*24*3600;
X = abs(t-t(end))*ucurrent/1000; % t is in sec.
plot(X,SOD, 'k')
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
ylabel('F_b')

% s2
s2 = subplot('position', [.1 .7 .72 .15]);
load test_a.mat
V = O2min:5:O2max;
contourf(X, P, Tmat)
hold on
[c, h] = contour(X, P, Tmat, V, 'color', 'k');
set(h,'ShowText','on', 'textstep', get(h,'LevelStep')*2);
hold off
set(gca, 'ydir', 'reverse')
cmap = lbmap(128,'RedBlue');
colormap(flipud(cmap));
caxis([O2min O2max])
ylim([200 300])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])

% s2.2
s22 = subplot('position', [.85 .7 .1 .15]);
plot(Tmat(:,1), P, '--k')
hold on
plot(Tmat(:,end), P, 'k')
hold off
set(gca, 'ydir', 'reverse')
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
xlim([0 220])
ylim([200 300])


% s3
s3 = subplot('position', [.1 .4 .72 .25]);
load test_b.mat
V = O2min:10:O2max;
contourf(X, P, Tmat)
hold on
[c, h] = contour(X, P, Tmat, V, 'color', 'k');
set(h,'ShowText','on', 'textstep', get(h,'LevelStep')*2);
hold off
set(gca, 'ydir', 'reverse')
cmap = lbmap(128,'RedBlue');
colormap(flipud(cmap));
caxis([O2min O2max])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])
ylabel('Depth (m)')
ylim([150 300])


% s3.3
s33 = subplot('position', [.85 .4 .1 .25]);
plot(Tmat(:,1), P, '--k')
hold on
plot(Tmat(:,end), P, 'k')
hold off
set(gca, 'ydir', 'reverse')
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
xlim([0 220])
ylim([150 300])

% s4
s4 = subplot('position', [.1 .1 .72 .25]);
load test_c.mat
V = O2min:20:O2max;
contourf(X, P, Tmat)
hold on
[c, h] = contour(X, P, Tmat, V, 'color', 'k');
set(h,'ShowText','on', 'textstep', get(h,'LevelStep')*2);
hold off
set(gca, 'ydir', 'reverse')
cmap = lbmap(128,'RedBlue');
colormap(flipud(cmap));
caxis([O2min O2max])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
xlabel('Distance (km)')
ylim([150 300])


% s4.4
s44 = subplot('position', [.85 .1 .1 .25]);
plot(Tmat(:,1), P, '--k')
hold on
plot(Tmat(:,end), P, 'k')
hold off
set(gca, 'ydir', 'reverse')
set(gca, 'yticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
xlim([0 220])
ylim([150 300])

xlabel('[O_2 (umol L^{-1})]')

print('-dpng', '-r300','O2_testcases.png')

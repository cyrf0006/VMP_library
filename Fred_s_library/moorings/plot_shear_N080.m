figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 14])   

subplot('position', [.1 .8 .78 .15])
I = find(mtime>t1 & mtime<t2);
plot(mtime(I), level(I))
set(gca, 'xgrid', 'on')
set(gca, 'xtick', [t1:3/24:t2])
set(gca, 'xticklabel', [])
xlim([t1 t2])
ylabel('\eta (m)')

subplot('position', [.1 .45 .78 .3])
I = find(time_vec >= t1 & time_vec <= t2);
imagesc(time_vec(I), abs(z_adcp_S-83), S2_ave(:,I));
set(gca, 'ydir', 'normal')
set(gca, 'xgrid', 'on')
set(gca, 'xtick', [t1:3/24:t2])
set(gca, 'xticklabel', []);
ylabel('hab(m)')
xlim([t1 t2])
ylim([0 50])
c = colorbar('position', [.89 .5 .025 .2]);
caxis([0.05 .35])
ti = ylabel(c,'S^2 (m^2 s^{-1})', 'FontSize', 10);

subplot('position', [.1 .1 .78 .3])
load '../M_N080/Tfield_M_N080.mat'  % this will erase time_vec!!!
I = find(time_vec_N080 >= t1 & time_vec_N080 <= t2);
imagesc(time_vec_N080(I), hab, Tgrid(:, I))
set(gca, 'ydir', 'normal')
%datetick('x', 7)
datetick('x', 15)
set(gca, 'xgrid', 'on')
set(gca, 'xtick', [t1:1:t2])
caxis([0 3])
xlim([t1 t2])
ylim([0 50])
%xlabel('sept-2011')
xlabel('28 sept')
ylabel('hab (m)')
c = colorbar('position', [.89 .15 .025 .2]);
ti = ylabel(c,'{T(^{\circ}C)}', 'FontSize', 10);

print('-dpng', '-r300', 'plot_shear_N080.png')


[valong, vcross]=rotate_vecd(E_ave, N_ave, theta);


figure(5)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 20])   

subplot('position', [.22 .85 .65 .12])
I = find(mtime>t1 & mtime<t2);
plot(mtime(I), level(I))
set(gca, 'xgrid', 'on')
set(gca, 'xtick', [t1+1:1:t2-1])
%set(gca, 'yticklabel', []);
%set(gca, 'xtick', [t1:3/24:t2])
set(gca, 'xticklabel', [])
xlim([t1 t2])
ylabel('\eta (m)')

subplot('position', [.22 .60 .65 .21])
I = find(time_vec >= t1 & time_vec <= t2);
imagesc(time_vec(I), abs(z_adcp-83), valong(:,I)*100);
set(gca, 'ydir', 'normal')
set(gca, 'xgrid', 'on')
%set(gca, 'xtick', [t1:3/24:t2])
set(gca, 'xtick', [t1+1:1:t2-1])
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
%ylabel('hab(m)')
xlim([t1 t2])
ylim([0 23])
caxis([-40 40])


subplot('position', [.22 .35 .65 .21])
imagesc(time_vec(I), abs(z_adcp-83), vcross(:,I)*100);
set(gca, 'ydir', 'normal')
set(gca, 'xgrid', 'on')
%set(gca, 'xtick', [t1:3/24:
set(gca, 'xtick', [t1+1:1:t2-1])
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
%ylabel('hab(m)')
xlim([t1 t2])
ylim([0 23])
c = colorbar('position', [.89 .5 .025 .2]);
caxis([-40 40])
ti = ylabel(c,'u,v (cm s^{-1})', 'FontSize', 10);

subplot('position', [.22 .1 .65 .21])
imagesc(time_vec(I), abs(z_adcp_S-83), S2_ave(:,I));
set(gca, 'ydir', 'normal')
set(gca, 'xtick', [t1:1:t2])
set(gca, 'xgrid', 'on')
datetick('x', 7)
%ylabel('hab(m)')
set(gca, 'yticklabel', []);
xlabel('sept. / oct.')
xlim([t1 t2])
ylim([0 23])
c = colorbar('position', [.89 .125 .025 .15]);
caxis([0.05 .25])
ti = ylabel(c,'S^2 (m^2 s^{-2})', 'FontSize', 10);


% Residual circul
res_along = nanmean(valong(:,I), 2);
res_cross = nanmean(vcross(:,I), 2);
res_shear = nanmean(S2_ave(:,I), 2);

subplot('position', [.07 .6 .12 .21])
plot(res_along*100, abs(z_adcp-83), 'k')
set(gca, 'ydir', 'normal')
% $$$ xlab = xlabel('u(cm s^{-1})');
% $$$ posx = get(xlab, 'pos');
% $$$ posx(2) = posx(2)+1;
% $$$ set(xlab, 'pos', posx)
% $$$ ylabel('hab (m)')
ylim([0 23])
xlim([-10 0])

subplot('position', [.07 .35 .12 .21])
plot(res_cross*100, abs(z_adcp-83), 'k')
set(gca, 'ydir', 'normal')
% $$$ xlab = xlabel('v(cm s^{-1})');
% $$$ posx = get(xlab, 'pos');
% $$$ posx(2) = posx(2)+1;
% $$$ set(xlab, 'pos', posx)
ylabel('hab (m)')
ylim([0 23])
xlim([-5 0])

subplot('position', [.07 .1 .12 .21])
plot(res_shear, abs(z_adcp_S-83), 'k')
set(gca, 'ydir', 'normal')
% $$$ xlab = xlabel('S^2 (m^2 s^{-2})');
% $$$ ylabel('hab (m)')
ylim([0 23])
xlim([0.05 0.15])

print('-dpng', '-r300', 'UVS_N080_residual.png')

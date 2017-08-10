[valong, vcross]=rotate_vecd(E_ave, N_ave, theta);
xtik = [-6:2:6];

FONTSIZE = 10;

load ~/WINDEX/data_processing/sept_2011_mission/Mouillages/T_regtide.mat


figure(5)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 16 15])   

subplot('position', [.22 .76 .65 .18])
imagesc(reg_tide, hab, T_tide)
set(gca, 'ydir', 'normal')
set(gca, 'xgrid', 'on')
set(gca, 'xtick', xtik)
%set(gca, 'yticklabel', []);
set(gca, 'xticklabel', [])
xlim([t1 t2])
%ylabel('hab (m)')
c = colorbar('position', [.89 .78 .025 .14]);
%caxis([-40 40])
ti = ylabel(c,'T(^{\circ}C)', 'FontSize', 10);
set(gca, 'fontsize',FONTSIZE)

subplot('position', [.22 .54 .65 .18])
I = find(time_vec >= t1 & time_vec <= t2);
imagesc(time_vec(I), abs(z_adcp-83), valong(:,I)*100);
set(gca, 'ydir', 'normal')
set(gca, 'xgrid', 'on')
set(gca, 'xtick', xtik)
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
%ylabel('hab(m)')
c = colorbar('position', [.89 .56 .025 .14]);
caxis([-30 30])
ti = ylabel(c,'u (cm s^{-1})', 'FontSize', 10);
xlim([t1 t2])
ylim([0 23])
set(gca, 'fontsize',FONTSIZE)


subplot('position', [.22 .32 .65 .18])
imagesc(time_vec(I), abs(z_adcp-83), vcross(:,I)*100);
set(gca, 'ydir', 'normal')
set(gca, 'xgrid', 'on')
set(gca, 'xtick', xtik)
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
%ylabel('hab(m)')
xlim([t1 t2])
ylim([0 23])
c = colorbar('position', [.89 .34 .025 .14]);
caxis([-5 5])
ti = ylabel(c,'v (cm s^{-1})', 'FontSize', 10);
set(gca, 'fontsize',FONTSIZE)

subplot('position', [.22 .1 .65 .18])
%imagesc(time_vec(I), abs(z_adcp_S-83), S2_ave(:,I));
contourf(time_vec(I), abs(z_adcp_S-83), S2_ave(:,I), 'linestyle', 'none');
set(gca, 'ydir', 'normal')
set(gca, 'xtick', xtik)
set(gca, 'xgrid', 'on')
%datetick('x', 7)
%ylabel('hab(m)')
set(gca, 'yticklabel', []);
xlabel('time to hightide (hour)')
xlim([t1 t2])
ylim([0 23])
c = colorbar('position', [.89 .12 .025 .14]);
%caxis([0.05 .25])
ti = ylabel(c,'log(S^2) (s^{-2})', 'FontSize', 10);
tiPos = get(ti, 'pos');
tiPos(1) = 9;
set(ti, 'pos', tiPos)
set(gca, 'fontsize',FONTSIZE)


% Residual circul
res_along = nanmean(valong(:,I), 2);
res_cross = nanmean(vcross(:,I), 2);
res_shear = nanmean(S2_ave(:,I), 2);

subplot('position', [.07 .54 .12 .18])
plot(res_along*100, abs(z_adcp-83), 'k')
set(gca, 'ydir', 'normal')
% $$$ xlab = xlabel('u(cm s^{-1})');
% $$$ posx = get(xlab, 'pos');
% $$$ posx(2) = posx(2)+1;
% $$$ set(xlab, 'pos', posx)
% $$$ ylabel('hab (m)')
ylim([0 23])
xlim([-10 0])
set(gca, 'fontsize',FONTSIZE)

subplot('position', [.07 .32 .12 .18])
plot(res_cross*100, abs(z_adcp-83), 'k')
set(gca, 'ydir', 'normal')
% $$$ xlab = xlabel('v(cm s^{-1})');
% $$$ posx = get(xlab, 'pos');
% $$$ posx(2) = posx(2)+1;
% $$$ set(xlab, 'pos', posx)
ylabel('hab (m)')
ylim([0 23])
xlim([-5 0])
set(gca, 'fontsize',FONTSIZE)

subplot('position', [.07 .1 .12 .18])
plot(res_shear, abs(z_adcp_S-83), 'k')
set(gca, 'ydir', 'normal')
% $$$ xlab = xlabel('S^2 (m^2 s^{-2})');
% $$$ ylabel('hab (m)')
ylim([0 23])
%xlim([0.05 0.15])
set(gca, 'fontsize',FONTSIZE)

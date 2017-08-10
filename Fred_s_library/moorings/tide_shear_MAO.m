% just search and replace "riki" for "n080"
% and run in /home/cyrf0006/WINDEX/data_processing/sept_2011_mission/Mouillages/tide_shear
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 8])   

subplot('position', [.15 .15 .8 .35])
load subplot_riki.mat
errorbar(reg_tide, S2_boot_dot, S2_error, '.k')
hold on
plot(time_shear, smooth_shear1, 'k', 'linewidth', 2)
hold off
xlim([-6.5 6.5])
ylabel('S^2 (m^2 s^{-2})')
xlabel('time versus high tide (hour)')
set(gca, 'xgrid', 'on')
set(gca, 'fontsize', 10)
%set(gca, 'yscale', 'log')
%ylim([.082 .09])
text(6,0.0415,'b', 'fontsize',10, 'fontweight', 'bold')
text(6,1.425e-4,'b', 'fontsize',10, 'fontweight', 'bold')


subplot('position', [.15 .6 .8 .35])
load('eps_tide_riki.mat')

errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
set(gca, 'yscale', 'log')
hold on
plot(time_eps, smooth_eps1, 'k', 'linewidth', 2)
load subplot_riki_spring.mat
plot(time_eps, smooth_eps1, 'k', 'linewidth', 1)
load subplot_riki_neap.mat
plot(time_eps, smooth_eps1, '--k', 'linewidth', 1)
hold off

ylabel('\epsilon (W kg^{-1})')
xlim([-6.5 6.5])
ylim([1e-9 1e-6])
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'fontsize', 10)
text(6,2e-9,'a', 'fontsize',10, 'fontweight', 'bold')


print('-dpng', '-r300', 'S2_vs_M2_riki3.png')
set(gcf, 'renderer', 'painters')
print('-depsc', 'S2_vs_M2_riki3.eps')

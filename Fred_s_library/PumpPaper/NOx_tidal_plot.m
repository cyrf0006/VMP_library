figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 13])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.02; % very right of figuer
tops = 0.03; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
FS = 14;

subplot(311)
semilogy(time2, KVec, '.k')
hold on
patch([-6 6 6 -6 -6], [1.1 1.1 13 13 1.1]*10^(-2), [1 1 1]*.7, 'linestyle', 'none')
plot([-6, 6], [4.5e-2 4.5e-2], '--k')
semilogy(time2, KVec, '.k', 'markerSize', 10)
xlim([-6 6])
ylim([1e-7 1e0])
set(gca, 'ytick', [1e-8 1e-6 1e-4 1e-2 1e0])
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])
ylabel('$K (\rm m^2 s^{-~~1})$', 'interpreter', 'latex')
text(-5.7, 5e-7, 'a', 'fontSize', FS, 'fontWeight', 'bold')
adjust_space

subplot(312)
semilogy(time2, dNVec, '.k', 'markerSize', 10)
hold on
patch([-6 6 6 -6 -6], [.2 .2 .43 .43 .2], [1 1 1]*.7, 'linestyle', 'none')
plot([-6, 6], [.3 .3], '--k')
semilogy(time2, dNVec, '.k', 'markerSize', 10)
ylim([1e-3 1e4])
xlim([-6 6])
set(gca, 'ytick', [1e-4 1e-2 1e0 1e2 1e4])
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])
ylabel('$\frac{dC_{NO_3}}{dz} (\rm mmol~m^{-~~4})$', 'interpreter', 'latex')
text(-5.7, 5e-3, 'b', 'fontSize', FS, 'fontWeight', 'bold')
adjust_space


subplot(313)
semilogy(time2, FVec*86400, '.k')
hold on
patch([-6 6 6 -6 -6], [24 24 400 400 24], [1 1 1]*.7, 'linestyle', 'none')
plot([-6, 6], [121 121], '--k')
semilogy(time2, FVec*86400, '.k', 'markerSize', 10)
ylim([1e-3 1e4])
xlim([-6 6])
set(gca, 'ytick', [1e-4 1e-2 1e0 1e2 1e4])
set(gca, 'ygrid', 'on')
ylabel('$F_{NO_3} (\rm mmol m^{-~~2} d^{-~~1})$', 'interpreter', 'latex')
xlabel('t_{HT} (h)')
text(-5.7, 5e-3, 'c', 'fontSize', FS, 'fontWeight', 'bold')
adjust_space


set(gcf, 'renderer', 'painters')
print('-depsc2', 'CNOx_tidal.eps')

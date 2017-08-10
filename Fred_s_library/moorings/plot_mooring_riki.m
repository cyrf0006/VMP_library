
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 4; % no. subplot row
dx = 0.001 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 25 25])
   
xtickmark = timelim(1):timelim(2);

subplot(4,1,1)
plot(time_tide, tide)
ylabel('\eta (m)')
xlim([timelim])
set(gca, 'xtick', xtickmark)
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
title('Predicted tide heights at Rimouski')
adjust_space

subplot(4,1,2)
plot(time_wind, wind)
ylabel('U (m/s)')
xlim([timelim])
set(gca, 'xtick', xtickmark)
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
title('Wind speed at Ile Bicquette (6-h filtered)')
adjust_space

subplot(4,1,3)
[AX,H1,H2] = plotyy(time_rbr, T_rbr, time_rbr, P_rbr);
%plot(time_wind, wind)
set(get(AX(1),'Ylabel'),'String','T(^{\circ}C)')
set(get(AX(2),'Ylabel'),'String','P(dbar)') 
xlim(AX(1),[timelim(1) timelim(2)]);
xlim(AX(2),[timelim(1) timelim(2)]);
ylim(AX(1),[0 3]);
ylim(AX(2),[60 70]);
set(AX(1),'YTick',[0 1 2 3])
set(AX(2),'YTick',[60 62 64 66 68 70])
set(AX(1),'ygrid','on')
set(AX(2),'ygrid','on')
set(AX, 'xTick',xtickmark) % delete the x-tick labels
set(AX, 'xTickLabel','') % delete the x-tick labels
set(gca, 'xgrid', 'on')
title('RBR at 65m')
adjust_space

subplot(4,1,4)
imagesc(time_adcp, z, Imean) 
set(gca, 'ydir', 'reverse')
xlim([timelim])
set(gca, 'xtick', xtickmark)
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
ylim([60 90])
caxis([60 180])
datetick('x', 7)
xlim([timelim])
xlabel('Date (sept-oct 2011)')
ylabel('d(m)')
title('Mean backscatter intensity')
adjust_space

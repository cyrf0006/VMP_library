

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3+both; % no. subplot row
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

subplot(nrow,1,1)
plot(time_tide, tide, 'k')
ylabel('\eta (m)')
xlim([timelim])
set(gca, 'xtick', xtickmark)
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
%set(gca, 'xminortick', 'on')
%set(gca, 'xminorgrid', 'on')
%title('Predicted tide heights at Rimouski')
adjust_space

subplot(nrow,1,2)
[AX,H1,H2] = plotyy(time_wind, wind, time_wind, winddir);
ylabel('U (m/s)')
set(get(AX(1),'Ylabel'),'String','U(m s^{-1})')
set(get(AX(2),'Ylabel'),'String','\theta') 
%xlim([timelim])
xlim(AX(1),[timelim(1) timelim(2)]);
xlim(AX(2),[timelim(1) timelim(2)]);
ylim(AX(2),[0 360]);
set(AX(2), 'yTick',[0 90 180 270])
set(AX(2), 'yTicklabel',['N'; 'E'; 'S'; 'W';])
%set(gca, 'xtick', xtickmark)
set(AX, 'xTick',xtickmark) 
%set(gca, 'xticklabel', [])
set(AX, 'xTickLabel','') % delete the x-tick labelsset(gca, 'xgrid', 'on')
set(AX(1),'ygrid','on')
set(AX(2),'ygrid','on')
%set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
%set(gca, 'xminortick', 'on')
%set(gca, 'xminorgrid', 'on')
%title('Wind at Ile Bicquette (6-h filtered)')
adjust_space

% T at 65m
% $$$ subplot(4,1,3)
% $$$ [AX,H1,H2] = plotyy(time_rbr, T_rbr, time_rbr, P_rbr);
% $$$ %plot(time_wind, wind)
% $$$ set(get(AX(1),'Ylabel'),'String','T(^{\circ}C)')
% $$$ set(get(AX(2),'Ylabel'),'String','P(dbar)') 
% $$$ xlim(AX(1),[timelim(1) timelim(2)]);
% $$$ xlim(AX(2),[timelim(1) timelim(2)]);
% $$$ ylim(AX(1),[0 3]);
% $$$ ylim(AX(2),[60 70]);
% $$$ set(AX(1),'YTick',[0 1 2 3])
% $$$ set(AX(2),'YTick',[60 62 64 66 68 70])
% $$$ set(AX(1),'ygrid','on')
% $$$ set(AX(2),'ygrid','on')
% $$$ set(AX, 'xTick',xtickmark) % delete the x-tick labels
% $$$ set(AX, 'xTickLabel','') % delete the x-tick labels
% $$$ set(gca, 'xgrid', 'on')
% $$$ title('RBR at 65m')
% $$$ keyboard
% $$$ adjust_space


if n080 == 1
    subplot(nrow,1,3)
    imagesc(time_vec_N080, hab, Tgrid_n080)
    xlim([timelim(1) timelim(2)]);
    %title('Interpolated temperature')
    set(gca, 'xtick', xtickmark)
    set(gca, 'xticklabel', [])
    set(gca, 'xgrid', 'on')
    set(gca, 'ygrid', 'on')
    set(gca, 'ydir', 'normal')
    ylabel('hab (m)')
    caxis([0 4])
    if both == 0
        datetick('x', 7)
        xlim([timelim(1) timelim(2)]);
        ylabel('sept/oct')
        %set(gca, 'xminortick', 'on')
        %set(gca, 'xminorgrid', 'on')
    end
    adjust_space
elseif riki == 1
    subplot(nrow,1,3)
    imagesc(time_vec_RIKI, zvec, Tgrid_riki)
    xlim([timelim(1) timelim(2)]);
    %title('Interpolated temperature')
    set(gca, 'xtick', xtickmark)
    set(gca, 'xticklabel', [])
    set(gca, 'xgrid', 'on')
    set(gca, 'ygrid', 'on')
    set(gca, 'ydir', 'reverse')
    ylabel('z (m)')
    caxis([0 4])
    if both == 0
        datetick('x', 7)
        xlim([timelim(1) timelim(2)]);
        ylabel('sept/oct')
        %set(gca, 'xminortick', 'on')
        %set(gca, 'xminorgrid', 'on')
    end
    adjust_space
end


if both ==1
    subplot(nrow,1,4)
    imagesc(time_vec_N080, hab, Tgrid_n080)
    xlim([timelim(1) timelim(2)]);
    %title('Interpolated temperature')
    set(gca, 'xtick', xtickmark)
    %set(gca, 'xticklabel', []);
    set(gca, 'xgrid', 'on')
    set(gca, 'ygrid', 'on')
    set(gca, 'ydir', 'normal')
    %set(gca, 'xminortick', 'on')
    %set(gca, 'xminorgrid', 'on')
    datetick('x', 7)
    xlim([timelim(1) timelim(2)]);
    ylabel('sept/oct')
    ylabel('hab (m)')
    caxis([0 4])
    adjust_space
end

% $$$ subplot(nrow,1,4)
% $$$ imagesc(time_adcp, z+59, Imean) 
% $$$ set(gca, 'ydir', 'reverse')
% $$$ datetick('x', 15)
% $$$ xlim([timelim])
% $$$ set(gca, 'xtick', xtickmark)
% $$$ set(gca, 'xticklabel', [])
% $$$ set(gca, 'xgrid', 'on')
% $$$ set(gca, 'ygrid', 'on')
% $$$ ylim([60 90])
% $$$ caxis([60 180])
% $$$ datetick('x', 7)
% $$$ xlim([timelim])
% $$$ xlabel('Date (sept-oct 2011)')
% $$$ ylabel('d(m)')
% $$$ title('Mean backscatter intensity')
% $$$ adjust_space


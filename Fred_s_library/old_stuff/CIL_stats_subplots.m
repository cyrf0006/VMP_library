clear
month = 4:11;
TIK = [datenum(999, month,1, 0,0,0)]; %ticklabel for plot
TIK15 = [datenum(999, month,15, 0,0,0)]; %ticklabel for plot
T_LIM = [datenum(999,4,15, 0,0,0) datenum(999,11,15, 0,0,0)]; %XLIM
                                                              %for
                                                              %plot(apr15
                                                              %-15nov)
T_LIM2 = [datenum(999,4,1, 0,0,0) datenum(999,11,30, 0,0,0)]; %XLIM for plot
no_days = T_LIM(2)-T_LIM(1);
no_days2 = T_LIM2(2)-T_LIM2(1);
%sub = 0; % 1 for subplot, 0 for independant plots
plotlegend = ['b'; 'r'; 'k'; 'y'; 'm'; 'g'; 'c'; 'b']; %for T-S plot

figw = 12; %cm width
figh = 4; %cm height

% Parameters for figure costumization
offset1 = 0.01; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 0.0005; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom
cbar_width = 0.03;
clabel_width = 0.08;
ti_cbar_frac = 15/16; %reduction of distance bet. colorbar and its
                      %title
ti_cbar_off1 = 3; % 300% from left of cbar
ti_cbar_off2 = 0.5; %middle of Cbar

%Ylab_offset = -18; %offset between yaxis and ylabel
Ylab_tightoffset = 0.09; %TightInset that must be added for title

% Different colormap
load BR_colormap;
load BR_colormap2;


for i = 1:length(month) 
    
    
    if month(i) <10
            tfilename  = sprintf('T_climato_0%d.dat', month(i));
            sfilename = sprintf('S_climato_0%d.dat', month(i));
            n2filename  = sprintf('N2_climato_0%d.dat', month(i));
            tstdfilename = sprintf('T_climatoSTD_0%d.dat', month(i));
            sstdfilename = sprintf('S_climatoSTD_0%d.dat', month(i));
            n2stdfilename = sprintf('N2_climatoSTD_0%d.dat', month(i));
            
            % model results
            tkobsfilename = sprintf('Tmodel_kobs_climato_0%d.dat', month(i));
            tkcstfilename = sprintf('Tmodel_kcst_climato_0%d.dat', month(i));
    else
            tfilename  = sprintf('T_climato_%d.dat', month(i));
            sfilename  = sprintf('S_climato_%d.dat', month(i));
            n2filename  = sprintf('N2_climato_%d.dat', month(i));
            tstdfilename = sprintf('T_climatoSTD_%d.dat', month(i));
            sstdfilename = sprintf('S_climatoSTD_%d.dat',month(i));
            n2stdfilename = sprintf('N2_climatoSTD_%d.dat', month(i));
            
            tkobsfilename = sprintf('Tmodel_kobs_climato_%d.dat', month(i));
            tkcstfilename = sprintf('Tmodel_kcst_climato_%d.dat', month(i));
    end
    
    tprofile = load(tfilename);
    sprofile = load(sfilename);
    n2profile = load(n2filename);
    tstdprofile = load(tstdfilename);
    sstdprofile = load(sstdfilename);
    n2stdprofile = load(n2stdfilename);
    
    % model results
    t_kobsprofile = load(tkobsfilename);
    t_kcstprofile = load(tkcstfilename);
    
    if i == 1
        depth = tprofile(:,1);
        depthN2 = n2profile(:,1);
    end
    
    n(i) = datenum(999, month(i), 15, 0,0,0); % climatology, 15th of the month, year 999
    n = n';
    Tmat(i,:) = tprofile(:,2);
    Smat(i,:) = sprofile(:,2);
    N2mat(i,:) = n2profile(:,2);
    TSTDmat(i,:) = tstdprofile(:,2);
    SSTDmat(i,:) = sstdprofile(:,2);
    N2STDmat(i,:) = n2stdprofile(:,2);
    
    % model results
    T_kobsmat(i,:) = t_kobsprofile(:,2);
    T_kcstmat(i,:) = t_kcstprofile(:,2);
       
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- plot preamble ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfig width/heigth
fw = 0.75;
fh = 0.18;
x0 = 0.1;
y0 = [.76 .54 .32 .1];

figure(1)
clf
%whitebg('black')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figw figh*3.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------   Tmin   -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('position', [x0 y0(1) fw fh])

load 'CIL_core.dat'
I = min(month):max(month);

plot(n, CIL_core(I,1), 'k', 'linewidth', 1)
xlim(T_LIM);
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
ylim([-2 2])

hold on
% shade area
y1 = CIL_core(I,1)+CIL_core(I,2);
y2 = CIL_core(I,1)-CIL_core(I,2);
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_core(I,1), 'k', 'linewidth', 1)

% fit curve
[P,S]=polyfit(I, CIL_core(I,1)', 1);
fit = P(2)+I.*P(1);
Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
R = corrcoef(fit, CIL_core(I,1)); % In fact here corr. makes no
                                  % sense because it uses too much pts
                                  
%text(n(6),0.1 ,{'\uparrow', '0.23^{\circ}C/month'},'VerticalAlignment','top', 'FontSize',8)
text(n(1)+15,.95*4-2 ,{'best fit: 0.23^{\circ}C mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
hold off

Ylab_pos = get(ylab, 'position') ;
ylab1 = Ylab_pos(1)-15;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(1)+.05*(Ylim(end)-Ylim(1)), 'a', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
% --------------------------------- %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Thickness ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
subplot('position', [x0 y0(2) fw fh])


plot(n, CIL_core(I,5), 'k', 'linewidth', 1)

ylab = ylabel('d (m)', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
ylim([0 120]);
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

% shade area
hold on
y1 = CIL_core(I,5)+CIL_core(I,6);
y2 = CIL_core(I,5)-CIL_core(I,6);
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_core(I,5), 'k', 'linewidth', 1)
hold off

% fit curve
[P,S]=polyfit(I, CIL_core(I,5)', 1);
fit = P(2)+I.*P(1);
hold on
Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);
hold off
set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
R = corrcoef(fit, CIL_core(I,1)); % In fact here corr. makes no sense because it uses too much pts

text(n(1)+15,.95*120 ,{'best fit: -9 m mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
hold off


Ylab_pos = get(ylab, 'position') ;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(1)+.05*(Ylim(end)-Ylim(1)), 'b', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
% --------------------------------- %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Heat content ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
subplot('position', [x0 y0(3) fw fh])

% load data
heat_int = load('HEAT_int_T.dat'); % in kJ/m3
heat_int_std = load('HEAT_int_std_T.dat');
% - If meanheat content (divided by thickness - %)
heat_int2 = heat_int./1000; %in MJ/m3
heat_int_std2 = heat_int_std./1000; 

plot(n, heat_int2(2,I), 'k', 'linewidth', 1)
xlim(T_LIM);
% - If mean (divided by thickness) - %
ylab = ylabel('{H (MJ m^{-3})}', 'fontsize', 10, 'VerticalAlignment', 'top');
% - If integrated - %
% $$$ ylab = ylabel('{H (GJ/m^3)}', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

% shade area
hold on
y1 = heat_int2(2,I)+heat_int_std2(2,I);
y2 = heat_int2(2,I)-heat_int_std2(2,I);
if size(n)~=size(y1)
    y1 = y1';
    y2 = y2';
end
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');

plot(n, heat_int2(2,I), 'k', 'linewidth', 1)
hold off

ylim([6 11])
%ylim([0 1])
% fit curve
[P,S]=polyfit(I, heat_int2(2,I), 1);
fit = P(2)+I.*P(1);
hold on

Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);

set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)

text(n(1)+15,.95*5+6 ,{'best fit: 0.6  MJ m^{-3} mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
hold off

Ylab_pos = get(ylab, 'position') ;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(1)+.05*(Ylim(end)-Ylim(1)), 'c', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
% --------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- core depth  ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
subplot('position', [x0 y0(4) fw fh])

plot(n, CIL_core(I,3), 'k', 'linewidth', 1)

ylab = ylabel('{Z (m)}', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

% shade area
hold on
y1 = CIL_core(I,3)+CIL_core(I,4);
y2 = CIL_core(I,3)-CIL_core(I,4);
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_core(I,3), 'k', 'linewidth', 1)
hold off
ylim([40 120]);

%replace datetick for custumization
for i=2:length(month)-1
    mm=datestr(n(i), 3);
    day_norm = (n(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

Ylab_pos = get(ylab, 'position') ;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(end)-.05*(Ylim(end)-Ylim(1)), 'd', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
% ------------    END   --------------------- %



% Save figure
set(gcf, 'renderer', 'painters')
print('-depsc2', 'CIL_plots.eps') % no colorbar
print('-dpng', '-r300',  'CIL_plots.png') % no colorbar




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% -------- For Quebec-Ocean  ----------- %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- plot preamble ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfig width/heigth
fw = 0.75;
fh = 0.25;
x0 = 0.1;
y0 = [.7 .4 .1];


figure(2)
clf
%whitebg('black')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figw figh*2.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------   Tmin   -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('position', [x0 y0(1) fw fh])

load 'CIL_core.dat'
I = min(month):max(month);

plot(n, CIL_core(I,1), 'k', 'linewidth', 1)
xlim(T_LIM);
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
ylim([-2 2])

hold on
% shade area
y1 = CIL_core(I,1)+CIL_core(I,2);
y2 = CIL_core(I,1)-CIL_core(I,2);
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_core(I,1), 'k', 'linewidth', 1)

% fit curve
[P,S]=polyfit(I, CIL_core(I,1)', 1);
fit = P(2)+I.*P(1);
Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
R = corrcoef(fit, CIL_core(I,1)); % In fact here corr. makes no
                                  % sense because it uses too much pts
                                  
%text(n(6),0.1 ,{'\uparrow', '0.23^{\circ}C/month'},'VerticalAlignment','top', 'FontSize',8)
%text(n(1)+15,.95*4-2 ,{'pente: 0.23^{\circ}C mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
text(n(1)+15,.95*4-2 ,{'0.23^{\circ}C mo^{-1}'},'VerticalAlignment','top', ...
     'HorizontalAlignment', 'left', 'FontSize',13, 'color', 'r', ...
     'fontweight', 'bold')
hold off

Ylab_pos = get(ylab, 'position') ;
ylab1 = Ylab_pos(1)-15;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 
% --------------------------------- %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Thickness ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
subplot('position', [x0 y0(2) fw fh])


plot(n, CIL_core(I,5), 'k', 'linewidth', 1)

ylab = ylabel('d (m)', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
ylim([0 120]);
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

% shade area
hold on
y1 = CIL_core(I,5)+CIL_core(I,6);
y2 = CIL_core(I,5)-CIL_core(I,6);
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_core(I,5), 'k', 'linewidth', 1)
hold off

% fit curve
[P,S]=polyfit(I, CIL_core(I,5)', 1);
fit = P(2)+I.*P(1);
hold on
Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);
hold off
set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
R = corrcoef(fit, CIL_core(I,1)); % In fact here corr. makes no sense because it uses too much pts

%text(n(1)+15,.95*120 ,{'pente: -9 m
%mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment',
%'left', 'FontSize',8)
text(n(1)+15,.95*120 ,{'-9 m mo^{-1}'},'VerticalAlignment','top', ...
     'HorizontalAlignment', 'left', 'FontSize',13, 'fontweight', ...
     'bold', 'color', 'r')
hold off


Ylab_pos = get(ylab, 'position') ;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 
% --------------------------------- %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Heat content ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
subplot('position', [x0 y0(3) fw fh])

% load data
heat_int = load('HEAT_int_T.dat'); % in kJ/m3
heat_int_std = load('HEAT_int_std_T.dat');
% - If meanheat content (divided by thickness - %)
heat_int2 = heat_int./1000; %in MJ/m3
heat_int_std2 = heat_int_std./1000; 

plot(n, heat_int2(2,I), 'k', 'linewidth', 1)
xlim(T_LIM);
% - If mean (divided by thickness) - %
ylab = ylabel('{H (MJ m^{-3})}', 'fontsize', 10, 'VerticalAlignment', 'top');
% - If integrated - %
% $$$ ylab = ylabel('{H (GJ/m^3)}', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

% shade area
hold on
y1 = heat_int2(2,I)+heat_int_std2(2,I);
y2 = heat_int2(2,I)-heat_int_std2(2,I);
if size(n)~=size(y1)
    y1 = y1';
    y2 = y2';
end
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');

plot(n, heat_int2(2,I), 'k', 'linewidth', 1)
hold off

ylim([6 11])
%ylim([0 1])
% fit curve
[P,S]=polyfit(I, heat_int2(2,I), 1);
fit = P(2)+I.*P(1);
hold on

Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);

set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)

%text(n(1)+15,.95*5+6 ,{'pente: 0.6  MJ m^{-3} mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
text(n(1)+15,.95*5+6 ,{'0.6  MJ m^{-3} mo^{-1}'}, ...
     'VerticalAlignment','top', 'HorizontalAlignment', 'left', ...
     'FontSize',13, 'fontweight', 'bold','color', 'r')
hold off

Ylab_pos = get(ylab, 'position') ;
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 

%replace datetick for custumization
for i=2:length(month)-1
    mm=datestr(n(i), 3);
    day_norm = (n(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end
% ------------    END   --------------------- %



% Save figure
set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'CIL_plots_QO.eps') % no colorbar
% $$$ print('-dpng', '-r300',  'CIL_plots_QO.png') % no colorbar
print('-depsc2', 'CIL_plots_agu.eps') % no colorbar

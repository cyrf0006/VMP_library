% Second version of a function used to subplot erosion variables of
% the CIl for the papaer submitted to JGR.

clear
month = 4:11;
TIK = [datenum(999, month,1, 0,0,0)]; %ticklabel for plot
TIK15 = [datenum(999, month,15, 0,0,0)]; %ticklabel for plot
T_LIM = [datenum(999,4,15, 0,0,0) datenum(999,11,15, 0,0,0)]; 

% some constants
rho_0 = 1.025e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
Tcil = 1; %degC, threshold for the CIL core
freeze_pt = -1.8;


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
    
    
    % -- Load and store bootstraped climatology -- %
    if month(i) <10
            tfilename  = sprintf('T_bootclim_0%d.dat', month(i));
            sfilename = sprintf('S_bootclim_0%d.dat', month(i));
            
            % model results
            %tkobsfilename = sprintf('Tmodel_kobs_climato_0%d.dat', month(i));
            %tkcstfilename = sprintf('Tmodel_kcst_climato_0%d.dat', month(i));
    else
            tfilename  = sprintf('T_bootclim_%d.dat', month(i));
            sfilename  = sprintf('S_bootclim_%d.dat', month(i));
             
            %tkobsfilename = sprintf('Tmodel_kobs_climato_%d.dat', month(i));
            %tkcstfilename = sprintf('Tmodel_kcst_climato_%d.dat', month(i));
    end
    
    tprofile = load(tfilename);
    sprofile = load(sfilename);
    
% $$$     % model results
% $$$     t_kobsprofile = load(tkobsfilename);
% $$$     t_kcstprofile = load(tkcstfilename);
    
    if i == 1
        depth = tprofile(:,1);
        dz = depth(2)-depth(1);
    end
    
    n(i) = datenum(999, month(i), 15, 0,0,0); % climatology, 15th of the month, year 999
    Tmat(i,:) = tprofile(:,2);
    Smat(i,:) = sprofile(:,2);
    TERRmat(i,:) = tprofile(:,3); %error on profile
    SERRmat(i,:) = sprofile(:,3);
    
% $$$     % model results
% $$$     T_kobsmat(i,:) = t_kobsprofile(:,2);
% $$$     T_kcstmat(i,:) = t_kcstprofile(:,2);
       
    % -- load Erosion variables -- %
    CIL_erosion = load('CIL_erosionCI.dat');
    
    
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

I = min(month):max(month);
%I = Imonth;

plot(n, CIL_erosion(1,I), 'k', 'linewidth', 1)
xlim(T_LIM);
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')

ylim([-2 2])

hold on
% shade area
y2 = [CIL_erosion(2,I)]';
y1 = [CIL_erosion(3,I)]';
n = n';

patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_erosion(1,I), 'k', 'linewidth', 1)

% fit curve
[P,S]=polyfit(I, CIL_erosion(1,I), 1);
fit = P(2)+I.*P(1);
disp(sprintf('Tmin slope is %d', P(1)));
Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
R = corrcoef(fit, CIL_erosion(1,I)); % In fact here corr. makes no
                                  % sense because it uses too much pts
                                  
%text(n(6),0.1 ,{'\uparrow', '0.23^{\circ}C/month'},'VerticalAlignment','top', 'FontSize',8)
text(n(1)+15,.95*4-2 ,{'best fit: 0.24^{\circ}C mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
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


plot(n, CIL_erosion(4,I), 'k', 'linewidth', 1)

ylab = ylabel('d (m)', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
ylim([0 120]);
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')

% shade area
hold on
y1 = [CIL_erosion(5,I)]';
y2 = [CIL_erosion(6,I)]';
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_erosion(4,I), 'k', 'linewidth', 1)
hold off

% fit curve
[P,S]=polyfit(I, CIL_erosion(4,I), 1);
fit = P(2)+I.*P(1);
disp(sprintf('thickness slope is %d', P(1)));
hold on
Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);
hold off
set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
R = corrcoef(fit, CIL_erosion(4,I)); % In fact here corr. makes no sense because it uses too much pts

text(n(1)+15,.95*120 ,{'best fit: -11m mo^{-1}'},'VerticalAlignment','top', 'HorizontalAlignment', 'left', 'FontSize',8)
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

plot(n, CIL_erosion(7,I), 'k', 'linewidth', 1)
xlim(T_LIM);
% - If mean (divided by thickness) - %
ylab = ylabel('{H (MJ m^{-3})}', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')

% shade area
hold on
y1 = [CIL_erosion(8,I)]';
y2 = [CIL_erosion(9,I)]';

if size(n)~=size(y1)
    y1 = y1';
    y2 = y2';
end

nann = find(isnan(y2)==1);
if ~isempty(nann)==1 % remove NaNs
    y1(nann)=[];
    y2(nann)=[];
    nn = n;
    nn(nann)=[];
    patch([nn; flipdim(nn,1); nn(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');    
    plot(n, CIL_erosion(7,I), 'k', 'linewidth', 1)
    hold off
    
    x =  CIL_erosion(7,I);
    II = I;
    x(nann) = [];
    II(nann) = [];
    
    [P,S]=polyfit(II, x, 1);
    fit = P(2)+II.*P(1);
    disp(sprintf('heat slope is %d', P(1)));
    hold on

    Rplot = plot(nn, fit);%, 'k', 'linewidth', 0.25);

    set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)
    
else
   
    patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(n, CIL_erosion(7,I), 'k', 'linewidth', 1)
    hold off 
    [P,S]=polyfit(I, CIL_erosion(7,I), 1);
    fit = P(2)+I.*P(1);
    disp(sprintf('heat slope is %d', P(1)));
    hold on

    Rplot = plot(n, fit);%, 'k', 'linewidth', 0.25);

    set(Rplot, 'color', [0 0 0], 'linewidth', 0.25)  
end

ylim([6 12])

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

plot(n, CIL_erosion(10,I), 'k', 'linewidth', 1)

ylab = ylabel('{Z (m)}', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')

% shade area
hold on
y1 = [CIL_erosion(11,I)]';
y2 = [CIL_erosion(12,I)]';
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_erosion(10,I), 'k', 'linewidth', 1)
hold off
ylim([40 120]);

%replace datetick for custumization
for i=2:length(month)-1
    mm=datestr(n(i), 3);
    day_norm = (n(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

Ylab_pos = get(ylab, 'position');
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ---  SECOND FIGURE --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
month = 4:11;
TIK = [datenum(999, month,1, 0,0,0)]; %ticklabel for plot
TIK15 = [datenum(999, month,15, 0,0,0)]; %ticklabel for plot
T_LIM = [datenum(999,4,15, 0,0,0) datenum(999,11,15, 0,0,0)]; 

% some constants
rho_0 = 1.025e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
Tcil = 1; %degC, threshold for the CIL core
freeze_pt = -1.8;


T_LIM2 = [datenum(999,4,1, 0,0,0) datenum(999,11,30, 0,0,0)]; %XLIM for plot
no_days = T_LIM(2)-T_LIM(1);
no_days2 = T_LIM2(2)-T_LIM2(1);
%sub = 0; % 1 for subplot, 0 for independant plots
plotlegend = ['b'; 'r'; 'k'; 'y'; 'm'; 'g'; 'c'; 'b']; %for T-S plot

figw = 13; %cm width
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

Ylab_offset = -18; %offset between yaxis and ylabel
%Ylab_tightoffset = 0.0; %TightInset that must be added for title



for i = 1:length(month) 
    
    
    % -- Load and store bootstraped climatology -- %
    if month(i) <10
            tfilename  = sprintf('T_bootclim_0%d.dat', month(i));
            sfilename = sprintf('S_bootclim_0%d.dat', month(i));
            
            % model results
            %tkobsfilename = sprintf('Tmodel_kobs_climato_0%d.dat', month(i));
            %tkcstfilename = sprintf('Tmodel_kcst_climato_0%d.dat', month(i));
    else
            tfilename  = sprintf('T_bootclim_%d.dat', month(i));
            sfilename  = sprintf('S_bootclim_%d.dat', month(i));
             
            %tkobsfilename = sprintf('Tmodel_kobs_climato_%d.dat', month(i));
            %tkcstfilename = sprintf('Tmodel_kcst_climato_%d.dat', month(i));
    end
    
    tprofile = load(tfilename);
    sprofile = load(sfilename);
    
% $$$     % model results
% $$$     t_kobsprofile = load(tkobsfilename);
% $$$     t_kcstprofile = load(tkcstfilename);
    
    if i == 1
        depth = tprofile(:,1);
        dz = depth(2)-depth(1);
    end
    
    n(i) = datenum(999, month(i), 15, 0,0,0); % climatology, 15th of the month, year 999
    Tmat(i,:) = tprofile(:,2);
    Smat(i,:) = sprofile(:,2);
    TERRmat(i,:) = tprofile(:,3); %error on profile
    SERRmat(i,:) = sprofile(:,3);
    
% $$$     % model results
% $$$     T_kobsmat(i,:) = t_kobsprofile(:,2);
% $$$     T_kcstmat(i,:) = t_kcstprofile(:,2);
       
    % -- load Erosion variables -- %
    CIL_erosion = load('CIL_erosionCI.dat');
    
    
    
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- plot preamble ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfig width/heigth
fw = 0.75;
fh = 0.18;
x0 = 0.1;
y0 = [.76 .54 .32 .1];




figure(2)
clf
%whitebg('black')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figw figh*3.5])

% -- param for text on graphics --%
n_fit = [datenum(999, 8, 1) datenum(999, 8, 31)];
y1_tmin = -0.4;% y value for  1-may-999
y1_thick = 70;
y1_heatc = 8;
slope_tmin = 0.24;
slope_thick = -11;
slope_heatc = 0.6;
y2_tmin = y1_tmin+slope_tmin;% y value for  31-may-999
y2_thick = y1_thick+slope_thick;
y2_heatc = y1_heatc+slope_heatc;
ytext_tmin = -0.5;
ytext_thick = 75;
ytext_heatc = 7.8;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------   Tmin   -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot('position', [x0 y0(1) fw fh])

I = min(month):max(month);

plot(n, CIL_erosion(1,I), 'k', 'linewidth', 1)
xlim(T_LIM);
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

ylim([-1 1.5])

hold on
% shade area
y2 = [CIL_erosion(2,I)]';
y1 = [CIL_erosion(3,I)]';
n = n';

patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_erosion(1,I), 'k', 'linewidth', 1)

plot(n_fit, [y1_tmin, y2_tmin], 'color', 'k',  'linewidth', 0.25)
%text(n_fit(1),ytext_tmin ,'${\rm \dot{T}_{min}=0.24\pm 0.04 ^{\circ}C mo^{-1}$}', 'interpreter','latex', 'VerticalAlignment','top', 'FontSize',8)
text(n_fit(1),ytext_tmin ,'T_{min}=0.24\pm 0.04 ^{\circ}C mo^{-1}','VerticalAlignment','top', 'FontSize',8)


hold off

Ylab_pos = get(ylab, 'position') ;
ylab1 = Ylab_pos(1)-12;
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


plot(n, CIL_erosion(4,I), 'k', 'linewidth', 1)

ylab = ylabel('d (m)', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
ylim([0 120]);
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

% shade area
hold on
y1 = [CIL_erosion(5,I)]';
y2 = [CIL_erosion(6,I)]';
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_erosion(4,I), 'k', 'linewidth', 1)

plot(n_fit, [y1_thick, y2_thick], 'color', 'k',  'linewidth', 0.25)
text(n_fit(1),ytext_thick ,{'d=-11\pm 2 m mo^{-1}'},'VerticalAlignment','bottom', 'FontSize',8)
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

plot(n, CIL_erosion(7,I), 'k', 'linewidth', 1)
xlim(T_LIM);
% - If mean (divided by thickness) - %
ylab = ylabel('{H (MJ m^{-3})}', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

% shade area
hold on
y1 = [CIL_erosion(8,I)]';
y2 = [CIL_erosion(9,I)]';

if size(n)~=size(y1)
    y1 = y1';
    y2 = y2';
end

nann = find(isnan(y2)==1);
if ~isempty(nann)==1 % remove NaNs
    y1(nann)=[];
    y2(nann)=[];
    nn = n;
    nn(nann)=[];
    patch([nn; flipdim(nn,1); nn(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');    
    plot(n, CIL_erosion(7,I), 'k', 'linewidth', 1)
    plot(n_fit, [y1_heatc y2_heatc], 'color', 'k',  'linewidth', 0.25)
    hold off
    
    
else
   
    patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(n, CIL_erosion(7,I), 'k', 'linewidth', 1)
    plot(n_fit, [y1_heatc y2_heatc], 'color', 'k',  'linewidth', 0.25)
    hold off 
end

ylim([6 12])
text(n_fit(1),ytext_heatc ,{'H=0.6 \pm 0.1 MJ m^{-3} mo^{-1}'},'VerticalAlignment','top', 'FontSize',8)
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

plot(n, CIL_erosion(10,I), 'k', 'linewidth', 1)

ylab = ylabel('{Z (m)}', 'fontsize', 10, 'VerticalAlignment', 'top');
xlim(T_LIM);
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

% shade area
hold on
y1 = [CIL_erosion(11,I)]';
y2 = [CIL_erosion(12,I)]';
patch([n; flipdim(n,1); n(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');
plot(n, CIL_erosion(10,I), 'k', 'linewidth', 1)
hold off
ylim([50 110]);

%replace datetick for custumization
for i=2:length(month)-1
    mm=datestr(n(i), 3);
    day_norm = (n(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

Ylab_pos = get(ylab, 'position');
set(ylab, 'position', [ylab1 Ylab_pos(2) Ylab_pos(3)]) 

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(end)-.05*(Ylim(end)-Ylim(1)), 'd', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
% ------------    END   --------------------- %



% Save figure
set(gcf, 'renderer', 'painters')
print('-depsc', 'CIL_plots_nofit.eps') % no colorbar
print('-dpng', '-r300',  'CIL_plots_nofit.png') % no colorbar



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ---  CONFIDENCE INTERVAL ON SLOPES--- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Computing confidence intervals...')
nboot = 1000;
I = min(month):max(month);
n_pts = 2; % number of point for confidence int.

% --------   Tmin   --------- %
ymean = [CIL_erosion(1,I)];
ymin = [CIL_erosion(2,I)];
ymax = [CIL_erosion(3,I)];

clear Y
for i=1:length(I)
    Y(:,i) = [ymin(i):(ymax(i)-ymin(i))/(n_pts-1):ymax(i)]';
    %Y(:,i) = [ymin(i) ymean(i) ymax(i)]';
end

N_month = size(Y,2);
N_CI = size(Y,1);

for b = 1:nboot
    
    for i=1:N_month
        r = rand(1,1);
        r = ceil(r*N_CI/1);
        y(i) = Y(r,i);
    end
    
    [P,S]=polyfit(I, y, 1);    
    Y_boot_b(b) = P(1); 

end

Y_boot_dot = nanmean(Y_boot_b);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

Y_boot_sort = sort(Y_boot_b);

% $$$ % 95%
% $$$ A = [Y_boot_dot Y_boot_sort(CI_2p5) Y_boot_sort(CI_97p5)];
% $$$ disp(sprintf(' -> Tmin 95p100 CI is: [%d  %d  %d]', A(1), A(2), A(3))); 
% $$$ 
% $$$ % 100%
% $$$ A = [Y_boot_dot Y_boot_sort(1) Y_boot_sort(end)];                        
% $$$ disp(sprintf(' -> Tmin 100p100 CI is: [%d  %d  %d]', A(1), A(2), A(3))); 

% standard error
Y_error = sqrt(sum((diff([Y_boot_b Y_boot_dot],1, 2)).^2, 2)./(nboot-1));
disp(sprintf(' -> Tmin error is: [%d  +- %d ]', Y_boot_dot, Y_error));

% --------- Thickness ------- %
ymean = [CIL_erosion(4,I)];
ymin = [CIL_erosion(5,I)];
ymax = [CIL_erosion(6,I)];

clear Y
for i=1:length(I)
    %Y(:,i) = [ymin(i):(ymax(i)-ymin(i))/(n_pts-1):ymax(i)]';
    Y(:,i) = [ymin(i) ymean(i) ymax(i)]';
end

N_month = size(Y,2);
N_CI = size(Y,1);

for b = 1:nboot
    
    for i=1:N_month
        r = rand(1,1);
        r = ceil(r*N_CI/1);
        y(i) = Y(r,i);
    end
    
    [P,S]=polyfit(I, y, 1);    
    Y_boot_b(b) = P(1); 

end

Y_boot_dot = nanmean(Y_boot_b);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

Y_boot_sort = sort(Y_boot_b);

% $$$ % 95%
% $$$ A = [Y_boot_dot Y_boot_sort(CI_2p5) Y_boot_sort(CI_97p5)];
% $$$ disp(sprintf(' -> Thickness 95p100 CI is: [%d  %d  %d]', A(1), A(2), A(3)));
% $$$ 
% $$$ % 100%
% $$$ A = [Y_boot_dot Y_boot_sort(1) Y_boot_sort(end)];                               
% $$$ disp(sprintf(' -> Thickness 100p100 CI is: [%d  %d  %d]', A(1), A(2), A(3)));

% standard error
Y_error = sqrt(sum((diff([Y_boot_b Y_boot_dot],1, 2)).^2, 2)./(nboot-1));
disp(sprintf(' -> Thickness error is: [%d  +- %d ]', Y_boot_dot, Y_error));


% -------- Heat content ------ %
ymean = [CIL_erosion(7,I)];
ymin = [CIL_erosion(8,I)];
ymax = [CIL_erosion(9,I)];

clear Y
for i=1:length(I)
    %Y(:,i) = [ymin(i):(ymax(i)-ymin(i))/(n_pts-1):ymax(i)]';
    Y(:,i) = [ymin(i) ymean(i) ymax(i)]';
end


N_month = size(Y,2);
N_CI = size(Y,1);

for b = 1:nboot
    
    for i=1:N_month
        r = rand(1,1);
        r = ceil(r*N_CI/1);
        y(i) = Y(r,i);
    end
    
    nann = find(isnan(y)==1);
    if ~isempty(nann)==1 % remove NaNs
        y(nann)=[];
        nn = I;
        nn(nann)=[];
    else
        nn=I;
    end

    [P,S]=polyfit(nn, y, 1);    
    Y_boot_b(b) = P(1); 

end

Y_boot_dot = nanmean(Y_boot_b);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

Y_boot_sort = sort(Y_boot_b);

% $$$ % 95%
% $$$ A = [Y_boot_dot Y_boot_sort(CI_2p5) Y_boot_sort(CI_97p5)];
% $$$ disp(sprintf(' -> Heat cont. 95p100 CI is: [%d  %d  %d]', A(1), A(2), A(3)))
% $$$ 
% $$$ % 100%
% $$$ A = [Y_boot_dot Y_boot_sort(1) Y_boot_sort(end)];                               
% $$$ disp(sprintf(' -> Heat cont. 100p100 CI is: [%d  %d  %d]', A(1), A(2), A(3)))

% standard error
Y_error = sqrt(sum((diff([Y_boot_b Y_boot_dot],1, 2)).^2, 2)./(nboot-1));
disp(sprintf(' -> Heat cont. error is: [%d  +- %d ]', Y_boot_dot, Y_error));

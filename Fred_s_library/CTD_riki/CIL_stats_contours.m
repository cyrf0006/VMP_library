
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

figw = 13; %cm width
figh = 4; %cm height

% Parameters for figure costumization
offset1 = 0.15; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 1.15; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom
cbar_width = 0.03;
clabel_width = 0.08;
ti_cbar_frac = 15/16; %reduction of distance bet. colorbar and its
                      %title
ti_cbar_off1 = 3; % 300% from left of cbar
ti_cbar_off2 = 0.5; %middle of Cbar
cbar_x0 = 0.1;
%Ylab_offset = -18; %offset between yaxis and ylabel
Ylab_tightoffset = 0.09; %TightInset that must be added for title

% customize subplot
out_dy = 1;
out_y01 = 0.7;
out_y02 = 0.4;
out_y03 = 0;
out_x0 = 0;
out_dx = 1;

pos_dy = 0.23;
pos_y01 = 0.7;
pos_y02 = 0.42;
pos_y03 = 0.15;
pos_x0 = 0.1;
pos_dx = 0.8;

ylab_offset = 10;
cbar_width = 0.02;
cbar_length = 0.6;
%cbar_offset = 0.01; % colorbar offset from figure
cbar_offset = 0.03;

% Different colormap
load BR_colormap;
load BR_colormap2;
load BRcolormap3;


%T_kcst_daily = load('T_diffus_daily_5p5e-5.dat');
T_kcst_daily = load('T_diffus_daily_Kveryveryall.dat'); % ATTENTION,PAS KCST
T_kobs_daily = load('T_diffus_daily_Kobs.dat');
N = load('N_diffus_daily.dat');

for i = 1:length(month) 
    
    
    if month(i) <10
            tfilename  = sprintf('T_bootclim_0%d.dat', month(i));
            %sfilename = sprintf('S_bootclim_0%d.dat', month(i));
    else
            tfilename  = sprintf('T_bootclim_%d.dat', month(i));
            %sfilename  = sprintf('S_bootclim_%d.dat', month(i));
    end
    
    % observations
    tprofile = load(tfilename);
    T = tprofile(:,2);
    depth = tprofile(:,1);
    
    % model results (pressure vec must be same as obser.)
    disp(sprintf('month %d', month(i)))
    II = find(str2num(datestr(N,5))==month(i));    
    
    Tmobs = nanmean(T_kobs_daily(:,II), 2);
    Tmcst = nanmean(T_kcst_daily(:,II), 2);

    if i == 1
        depth = tprofile(:,1);
    end
    
    n(i) = datenum(999, month(i), 15, 0,0,0); % climatology, 15th of the month, year 999
    n = n';
    Tmat(i,:) = tprofile(:,2);
    
    % model results
    T_kobsmat(i,:) = Tmobs;
    T_kcstmat(i,:) = Tmcst;
       
end


figure(1)
clf
%whitebg('black')
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 figw figh*3.5])


% ---------------------------- %
% --------- contour T -------- %
% ---------------------------- % 
subplot(3,1,1)
% Subplot calling
Out = [out_x0 out_y01 out_dx out_dy];
Pos = [pos_x0 pos_y01 pos_dx pos_dy];

%subplot(3, 1, 1, 'outerposition', Out, 'position', Pos);
%subplot(3, 1, 1, 'outerposition', Out);

contourf(n, depth, Tmat', [-1:0.1:5], 'linestyle', 'none')

% Xtra contour around CIL
hold on
contour(n, depth, Tmat',[1 1 1], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
hold off

set(gca, 'ydir', 'reverse')
xlim(T_LIM);
ylim([0 300])
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'xticklabel', [])
set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')  

%ylab = ylabel('{depth (m)}', 'fontsize', 10, 'VerticalAlignment', 'top');
ylab = ylabel(' ', 'fontsize', 10, 'VerticalAlignment', 'top');
ylabpos = get(ylab, 'position');
set(ylab, 'Pos', [ylabpos(1)-ylab_offset ylabpos(2) ylabpos(3)]);

% colorbar and its title
caxis([-2 5])
%colormap(BR_colormap)
colormap(mycolormap)

% dateticks
for i=2:length(month)-1
    mm=datestr(n(i), 3);
    day_norm = (n(i)-T_LIM(1))/no_days;
    text(day_norm,offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end


set(gca, 'fontsize', 10);

set(gca, 'Out', Out)
set(gca, 'Pos', Pos)

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(end)-.05*(Ylim(end)-Ylim(1)), 'a', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight','bold');

% ---------------------------- %
% -------- Tmodel - Kobs ------ %
% ---------------------------- % 
subplot(3,1,2)
Out(2) = out_y02;
Pos(2) = pos_y02;
%subplot(3, 1, 2, 'outerposition', Out, 'position', Pos);
%subplot(3, 1, 2, 'outerposition', Out);

contourf(n, depth, T_kobsmat', [-1:0.1:5], 'linestyle', 'none')

% Xtra contour around CIL
hold on
contour(n, depth, T_kobsmat',[1 1 1], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
hold off

set(gca, 'ydir', 'reverse')
xlim(T_LIM);
ylim([0 300])
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'XGrid', 'on')
set(gca, 'xticklabel', [])
set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')  

ylab = ylabel('P (dbar)', 'fontsize', 10, 'VerticalAlignment','top');
%no title but a space
%ylab = ylabel(' ', 'fontsize', 10, 'VerticalAlignment', 'top');
ylabpos = get(ylab, 'position');
set(ylab, 'Pos', [ylabpos(1)-ylab_offset ylabpos(2) ylabpos(3)]);

%colormap(BR_colormap);
colormap(mycolormap)
caxis([-2 5])
% $$$ c = colorbar('FontSize', 10, 'position', [Out(1)+Out(3)-(offset1+cbar_width+clabel_width) Pos(2)+offset2 cbar_width Pos(4)-2*offset2]);
% $$$ %c = colorbar('FontSize', 10, 'position', cbar_pos1);
% $$$ ti = ylabel(c,'{T(^{\circ}C)}', 'FontSize', 10);

set(gca, 'fontsize', 10);
set(gca, 'Out', Out)
set(gca, 'Pos', Pos)


% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(end)-.05*(Ylim(end)-Ylim(1)), 'b', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold');


% ---------------------------- % 
% -------- Tmodel - Kcst ------ %
% ---------------------------- % 
subplot(3,1,3)
Out(2) = out_y03;
Pos(2) = pos_y03;
%subplot(3, 1, 3, 'outerposition', Out, 'position', Pos);
%subplot(3, 1, 3, 'outerposition', Out);

contourf(n, depth, T_kcstmat', [-1:0.1:5], 'linestyle', 'none')

% Xtra contour around CIL
hold on
contour(n, depth, T_kcstmat',[1 1 1], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
hold off

set(gca, 'ydir', 'reverse')
xlim(T_LIM);
ylim([0 300])
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'XGrid', 'on')
set(gca, 'xticklabel', [])
set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')  

ylab = ylabel(' ', 'fontsize', 10, 'VerticalAlignment','top');
%no title but a space
%ylab = ylabel(' ', 'fontsize', 10, 'VerticalAlignment', 'top');
ylabpos = get(ylab, 'position');
set(ylab, 'Pos', [ylabpos(1)-ylab_offset ylabpos(2) ylabpos(3)]);

%colormap(BR_colormap);
colormap(mycolormap)
caxis([-2 5])
set(gca, 'fontsize', 10);

%colorbar
% $$$ c = colorbar('FontSize', 10, 'position', [Out(1)+Out(3)-(offset1+ ...
% $$$                                                   cbar_width+clabel_width) Pos(2)+offset2 cbar_width Pos(4)-2*offset2]);
c = colorbar('location', 'northoutside', 'fontsize', 10);
ti = ylabel(c,'{T(^{\circ}C)}', 'FontSize', 8);

set(c, 'position', [Pos(1)+cbar_x0 Pos(2)-cbar_offset-cbar_width ...
                    cbar_length cbar_width*.75]);
set(c, 'tickdir', 'out')

ti_pos = get(ti, 'position');
set(ti, 'rotation', 0)
set(ti, 'position', [1.5 ti_pos(2)-1 ti_pos(3)]); 
        
% $$$ %replace datetick for custumization
% $$$ for i=2:length(month)-1
% $$$     mm=datestr(n(i), 3);
% $$$     day_norm = (n(i)-T_LIM(1))/no_days;
% $$$     text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
% $$$          'center', 'verticalalignment', 'top', 'fontsize', 10)
% $$$ end

% --------------------------------- %

set(gca, 'Out', Out)
set(gca, 'Pos', Pos)

% Write letter identification
Ylim = ylim;
text(n(end)-4, Ylim(end)-.05*(Ylim(end)-Ylim(1)), 'c', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold');

set(gcf, 'renderer', 'painters'); % vectorial figure                              
print('-depsc2', 'CIL_contours.eps') % save without Xlabel


function obs_vs_model(month, model_Kobs, model_Kcst, model_N)

% function obs_vs_model(month, model_Kobs, model_Kcst, model_N)
% 
%    
% Used to plot a monthly climatology of CTD profiles with the
% uncertainty, together with the results of the diffusivity model
% for each monthly average.
%
% usage ex: 
% >> obs_vs_model([4 5 6 7 8 9 10 11], 'T_diffus_daily_Kobs.dat',
% 'T_diffus_daily_5p5e-5.dat', 'N_diffus_daily.dat')
% >> obs_vs_model([4 5 6 7 8 9 10 11], 'T_diffus_daily_Kobs.dat', 'T_diffus_daily_Kveryveryall.dat', 'N_diffus_daily.dat')
%
% If only one modeloutput is specified, the results will be a
% subplot for each month where the CTD observations are in a shade
% of gray while the modeled profile if superimposed. If another
% model output is specified, the second profile is also
% superimposed with a dashed line. Note that times of both output
% must be the same... (This is easy to adjust but I dont want to
% include too many parameters for the moment!)
%
% NOTE: to be run in  ~/WINDEX/data_processing/1D_mixing/
%
% Frederic Cyr - march 2011
%
% -------------------------------------------------------------- %

% % some constants
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.025e3;%kg/m^3
GAMMA = 0.2; %mixing efficiency
freeze_pt = -1.8; 
cp = 3.99; %Kj/Kg/K
CIL_def = 1;
dz=1;

% load files
%T=load('Tclim_matrix')';
TKobs=load(model_Kobs);
TKcst=load(model_Kcst);
Nmod=load(model_N)';
P = [1:300]';

offset3 = 0.04; % offset beteween figure and custom Xlabel 
                
%month = 4:11;
TIK = [datenum(999, month,1, 0,0,0)]; %ticklabel for plot
TIK15 = [datenum(999, month,15, 0,0,0)]; %ticklabel for plot
T_LIM = [datenum(999,4,15, 0,0,0) datenum(999,11,15, 0,0,0)]; %XLIM
                                                              %for
                                                              %plot(apr15                                                             %-15nov)
T_LIM2 = [datenum(999,4,1, 0,0,0) datenum(999,11,30, 0,0,0)]; %XLIM for plot
no_days = T_LIM(2)-T_LIM(1);
no_days2 = T_LIM2(2)-T_LIM2(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- Main loop (on months) ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on climatology
clear N
for i = 1:length(month)
    disp(sprintf('month %d', month(i)))
    II = find(str2num(datestr(Nmod,5))==month(i));
    
    % clim time vector
    N(i) = datenum(999, month(i), 15);
    
    if month(i) <10
        Tname  = sprintf('T_bootclim_0%d.dat', month(i));
        Sname  = sprintf('S_bootclim_0%d.dat', month(i));
               
    else
        Tname  = sprintf('T_bootclim_%d.dat', month(i));
        Sname  = sprintf('S_bootclim_%d.dat', month(i));
    end
    
    T = load(Tname);
    S = load(Sname);
    TmKo = nanmean(TKobs(:,II), 2);
    TmKc = nanmean(TKcst(:,II), 2);
    To_ave = T(:,2);
    To_2p5 = T(:,3);
    To_97p2 = T(:,4);
    So_ave = S(:,2);
    %So_err = S(:,3);
    P = T(:,1);   

    x1 = To_97p2;
    x2 = To_2p5;
    
    % --- Compute erosion variables --- %

    % T_min
    tcore_o(i) = min(To_ave);
    tcore_omin(i) = min(x2);
    tcore_omax(i) = min(x1);
    tcore_mko(i) = min(TmKo);
    tcore_mkc(i) = min(TmKc);
    
    % Thickness and Heat content
    CIL = find(To_ave < CIL_def);
    DENS = sw_dens(So_ave(CIL), To_ave(CIL), P(CIL));
    if ~isempty(CIL)==1
        thick_o(i) = length(CIL);
        heatc_o(i) =  cp*nansum(DENS.*(To_ave(CIL)-freeze_pt))*dz/length(CIL)/1000;%MJ    
    else
        thick_o(i) = 0;
        heatc_o(i) = NaN;
    end 
     
    CIL = find(x2 < CIL_def);
    DENS = sw_dens(So_ave(CIL), x2(CIL), P(CIL));
    if ~isempty(CIL)==1
        thick_omin(i) = length(CIL);
        heatc_omin(i) =  cp*nansum(DENS.*(x2(CIL)-freeze_pt))*dz/length(CIL)/1000;%MJ    
    else
        thick_omin(i) = 0;
        heatc_omin(i) = NaN;
    end
    
    CIL = find(x1 < CIL_def);
    DENS = sw_dens(So_ave(CIL), x1(CIL), P(CIL));
    if ~isempty(CIL)==1
        thick_omax(i) = length(CIL);
        heatc_omax(i) =  cp*nansum(DENS.*(x1(CIL)-freeze_pt))*dz/length(CIL)/1000;%MJ
    else
        thick_omax(i) = 0;
        heatc_omax(i) = NaN;
    end    

    CIL = find(TmKo < CIL_def);
    DENS = sw_dens(So_ave(CIL), TmKo(CIL), P(CIL));
    if ~isempty(CIL)==1
        thick_mko(i) = length(CIL);
        heatc_mko(i) =  cp*nansum(DENS.*(TmKo(CIL)-freeze_pt))*dz/length(CIL)/1000;%MJ
    else
        thick_mko(i) = 0;
        heatc_mko(i) = NaN;
    end    

    CIL = find(TmKc < CIL_def);
    DENS = sw_dens(So_ave(CIL), TmKc(CIL), P(CIL));
    if ~isempty(CIL)==1
        thick_mkc(i) = length(CIL);
        heatc_mkc(i) =  cp*nansum(DENS.*(TmKc(CIL)-freeze_pt))*dz/length(CIL)/1000;%MJ
    else
        thick_mkc(i) = 0;
        heatc_mkc(i) = NaN;
    end    
        
end %for i

% Rename vars...
Tmin_o = tcore_o;
Tmin_ko = tcore_mko;
Tmin_kc = tcore_mkc;
thick_ko = thick_mko;
thick_kc = thick_mkc;
heat_o = heatc_o;
heat_ko = heatc_mko;
heat_kc = heatc_mkc;


%%%%%%%%%%%%%%%%%%%%%%%%%
% -- find linear fit -- %
%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Tmin slope observed -- %
I = find(~isnan(Tmin_o)==1);
% -- for fit on grahp -- %
[PP,S]=polyfit(N(I), Tmin_o(I), 1);
Tmin_o_fit = PP(2)+N(I).*PP(1);
% -- real slope -- %
[PP,S]=polyfit(I, Tmin_o(I), 1);
slope_Tmin_o = PP(1);

% -- Tmin slope modeled_kobs -- %
I = find(~isnan(Tmin_ko)==1);
[PP,S]=polyfit(N(I), Tmin_ko(I), 1);
Tmin_ko_fit = PP(2)+N(I).*PP(1);
[PP,S]=polyfit(I, Tmin_ko(I), 1);
slope_Tmin_ko = PP(1);

% -- Tmin slope modeled_kcst -- %
I = find(~isnan(Tmin_kc)==1);
[PP,S]=polyfit(N(I), Tmin_kc(I), 1);
Tmin_kc_fit = PP(2)+N(I).*PP(1);
[PP,S]=polyfit(I, Tmin_kc(I), 1);
slope_Tmin_kc = PP(1);

% -- Thickness observed -- %
I = find(~isnan(thick_o)==1);
[PP,S]=polyfit(N(I), thick_o(I), 1);
thick_o_fit = PP(2)+N(I).*PP(1);
[PP,S]=polyfit(I, thick_o(I), 1);
slope_thick_o = PP(1);

% -- Thickness modeled_kobs -- %
I = find(~isnan(thick_ko)==1);
[PP,S]=polyfit(N(I), thick_ko(I), 1);
thick_ko_fit = PP(2)+N(I).*PP(1);
[PP,S]=polyfit(I, thick_ko(I), 1);
slope_thick_ko = PP(1);

% -- Thickness modeled_Kcst -- %
I = find(~isnan(thick_kc)==1);
[PP,S]=polyfit(N(I), thick_kc(I), 1);
thick_kc_fit = PP(2)+N(I).*PP(1);
[PP,S]=polyfit(I, thick_kc(I), 1);
slope_thick_kc = PP(1);

% -- heat content observed -- %
I = find(~isnan(heat_o)==1);
[PP,S]=polyfit(N(I), heat_o(I), 1);
heat_o_fit = PP(2)+N.*PP(1);
[PP,S]=polyfit(I, heat_o(I), 1);
slope_heat_o = PP(1);

% -- heat content modeled_kobs -- %
I = find(~isnan(heat_ko)==1);
[PP,S]=polyfit(N(I), heat_ko(I), 1);
heat_ko_fit = PP(2)+N.*PP(1);
[PP,S]=polyfit(I, heat_ko(I), 1);
slope_heat_ko = PP(1);

% -- heat content modeled_kcst -- %
I = find(~isnan(heat_kc)==1);
[PP,S]=polyfit(N(I), heat_kc(I), 1);
heat_kc_fit = PP(2)+N.*PP(1);
[PP,S]=polyfit(I, heat_kc(I), 1);
slope_heat_kc = PP(1);


%%%%%%%%%%%%%%%%%%%%%%%
% --  plot result  -- %
%%%%%%%%%%%%%%%%%%%%%%%
exy = 0.03;


% --- Version 1: all fits --- %
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 12 14])


% Tmin %
s1=subplot(3,1,1);
% Adjust space to minimize white space
pos1 = get(s1, 'position');
set(s1, 'position', [pos1(1) pos1(2)-exy pos1(3) pos1(4)+exy])

plot(N, Tmin_o, '*k')
hold on
plot(N, Tmin_ko, 'ok')
plot(N, Tmin_kc, 'dk')
legend('observ', '{K_i}', '{K_a}', ...
       'location', 'northwest')

Rplot = plot(N, Tmin_o_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, Tmin_ko_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, Tmin_kc_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
hold off

xlim(T_LIM);
ylim([-1 1.5])
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'box', 'on')
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')

ypos = get(ylab, 'position');
ypos1 = ypos(1)-15;
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_Tmin_o slope_Tmin_ko slope_Tmin_kc]


% Thickness %
s2=subplot(3,1,2);
pos2 = get(s2, 'position');
set(s2, 'position', [pos2(1) pos2(2)-exy/2 pos2(3) pos2(4)+exy])

plot(N, thick_o, '*k')
hold on
plot(N, thick_ko, 'ok')
plot(N, thick_kc, 'dk')

Rplot = plot(N, thick_o_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, thick_ko_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, thick_kc_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
hold off

set(gca, 'xticklabel', [])
xlim(T_LIM);
%ylim([-1 1])
ylab = ylabel('d(m)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'box', 'on')
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_thick_o slope_thick_ko slope_thick_kc]


% heat %
s3=subplot(3,1,3);
pos3 = get(s3, 'position');
set(s3, 'position', [pos3(1) pos3(2) pos3(3) pos3(4)+exy])


% mean heat content (divided by d)
plot(N, heat_o, '*k')
hold on
plot(N, heat_ko, 'ok')
plot(N, heat_kc, 'dk')

Rplot = plot(N, heat_o_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, heat_ko_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, heat_kc_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
hold off
set(gca, 'xticklabel', [])
xlim(T_LIM);
ylab = ylabel('H(MJ m^{-3})', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'box', 'on')
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

for i=2:length(month)-1
    mm=datestr(N(i), 3);
    day_norm = (N(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

[slope_heat_o slope_heat_ko slope_heat_kc]

set(gcf, 'renderer', 'painters')
print('-deps2', 'obs_vs_mod.eps')
% $$$ print('-dpng', '-r300', 'obs_vs_mod.png')


%--------------------------------------------- %
% --- Version 2: obs envelope + model fits --- %
% -------------------------------------------- %
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 12 14])


% Tmin %
s1=subplot(3,1,1);
% Adjust space to minimize white space
pos1 = get(s1, 'position');
set(s1, 'position', [pos1(1) pos1(2)-exy pos1(3) pos1(4)+exy])

% shade area
y1 = tcore_omin;
y2 = tcore_omax;

patch([N'; flipdim(N',1); N(1)'], [y1'; flipdim(y2',1); y1(1)'], [.8 .8 .8], 'edgecolor', 'none');

hold on
% $$$ plot(N, Tmin_ko, '.k')
% $$$ plot(N, Tmin_kc, '*k')

% $$$ legend('observations', 'K(z)', 'K=5\times10^{-5} m^2 s^{-1}', ...
% $$$        'location', 'northwest')

% $$$ Rplot = plot(N, Tmin_o_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, Tmin_ko_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, Tmin_kc_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
hold off

xlim(T_LIM);
ylim([-2 2])
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'box', 'on')
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')

ypos = get(ylab, 'position');
ypos1 = ypos(1)-15;
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

% Write letter identification
Ylim = ylim;
n=N;
text(n(1)+10, Ylim(1)+.02*(Ylim(end)-Ylim(1)), 'a', ...
         'verticalalignment', 'bottom', 'horizontalalignment', ...
     'right', 'fontsize', 10, 'fontweight','bold');

% Thickness %
s2=subplot(3,1,2);
pos2 = get(s2, 'position');
set(s2, 'position', [pos2(1) pos2(2)-exy/2 pos2(3) pos2(4)+exy])

% shade area
y1 = thick_omin;
y2 = thick_omax;
patch([N'; flipdim(N',1); N(1)'], [y1'; flipdim(y2',1); y1(1)'], [.8 .8 .8], 'edgecolor', 'none');

hold on
% $$$ plot(N, thick_ko, '.k')
% $$$ plot(N, thick_kc, '*k')

% $$$ legend('observations', 'K(z)', 'K=5\times10^{-5} m^2 s^{-1}', ...
% $$$        'location', 'northwest')

% $$$ Rplot = plot(N, thick_o_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, thick_ko_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, thick_kc_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
hold off

set(gca, 'xticklabel', [])
xlim(T_LIM);
ylim([0 120])
ylab = ylabel('d(m)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'box', 'on')
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

% Write letter identification
Ylim = ylim;
n=N;
text(n(1)+10, Ylim(1)+.02*(Ylim(end)-Ylim(1)), 'b', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight','bold');

% heat %
s3=subplot(3,1,3);
pos3 = get(s3, 'position');
set(s3, 'position', [pos3(1) pos3(2) pos3(3) pos3(4)+exy])

% shade area
y1 = heatc_omin;
y2 = heatc_omax;
I = find(~isnan(y2)==1);
patch([N(I)'; flipdim(N(I)',1); N(1)'], [y1(I)'; flipdim(y2(I)',1); y1(1)'], [.8 .8 .8], 'edgecolor', 'none');

hold on
% $$$ plot(N, heat_ko, '.k')
% $$$ plot(N, heat_kc, '*k')

% $$$ legend('observations', 'K(z)', 'K=5\times10^{-5} m^2 s^{-1}', ...
% $$$        'location', 'northwest')

% $$$ Rplot = plot(N, heat_o_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, heat_ko_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, heat_kc_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
hold off

set(gca, 'xticklabel', [])
xlim(T_LIM);
ylab = ylabel('H(MJ m^{-3})', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'box', 'on')
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
set(gca, 'yminortick', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

for i=2:length(month)-1
    mm=datestr(N(i), 3);
    day_norm = (N(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end


% Write letter identification
Ylim = ylim;
n=N;
text(n(1)+10, Ylim(1)+.02*(Ylim(end)-Ylim(1)), 'c', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight','bold');

set(gcf, 'renderer', 'painters')
print('-deps2', 'obs_vs_mod_shade.eps')
% $$$ print('-dpng', '-r300', 'obs_vs_mod.png')
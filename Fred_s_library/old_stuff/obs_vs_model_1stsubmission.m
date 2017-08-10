% % some constants
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.025e3;%kg/m^3
GAMMA = 0.2; %mixing efficiency
freeze_pt = -1.8; 
cp = 3.99; %Kj/Kg/K
CIL_def = 1;
dz=1;

T=load('Tclim_matrix')';
Tobs=load('T_kobsclim_matrix')';
Tcst=load('T_kcstclim_matrix')';
N=load('tclim_vector')';
P = [1:300]';

offset3 = 0.04; % offset beteween figure and custom Xlabel 
month = 4:11;
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

for i=1:length(N)

    % ------------------------------------------------- %
    % ----------------- Observations ------------------ %
    
    I = find(T(:,i) <= CIL_def);
    
    % Tmin
    Tmin_o(i) = min(T(:,i)); 
        
    if ~isempty(I)==1
        % thickness
        thick_o(i) = length(I);
    
        % heat content
        heat_o(i) = cp*nansum(rho_0.*(T(I,i)-freeze_pt))*dz/length(I);%in KJ/m3
                                                                      
        %heat_o(i) = cp*nansum(rho_0.*(T(I,i)-freeze_pt))*dz;%in kj/m2
       
    else
       thick_o(i)=0; 
       heat_o(i)=NaN;
    end
    % -------------------------------------------------- %
    

    % ------------------------------------------------- %
    % ------------------ K observed ------------------- %
    
    I = find(Tobs(:,i) <= CIL_def);
    
    % Tmin
    Tmin_ko(i) = min(Tobs(:,i)); 
    
    if ~isempty(I)==1    
        % thickness
        thick_ko(i) = length(I);
    
        % heat content
        heat_ko(i) = cp*nansum(rho_0.*(Tobs(I,i)-freeze_pt))*dz/length(I); %in KJ/m3
        %heat_ko(i) = cp*nansum(rho_0.*(Tobs(I,i)-freeze_pt))*dz; %in kj/m2
    else
       thick_ko(i)=0; 
       heat_ko(i)=NaN;
    end
    % -------------------------------------------------- %
   
      
    % -------------------------------------------------- %
    % ------------------ K constant ------------------- %
    
    I = find(Tcst(:,i) <= CIL_def);
    
    % Tmin
    Tmin_kc(i) = min(Tcst(:,i)); 
    
    if ~isempty(I)==1    
        % thickness
        thick_kc(i) = length(I);

        % heat content
        heat_kc(i) = cp*nansum(rho_0.*(Tcst(I,i)-freeze_pt))*dz/length(I); %in KJ/m3
        %heat_kc(i) = cp*nansum(rho_0.*(Tcst(I,i)-freeze_pt))*dz; %in kj/m2
    else
        thick_kc(i)=0; 
        heat_kc(i)=NaN;
    end
    % -------------------------------------------------- %   
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% -- find linear fit -- %
%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Tmin slope observed -- %
I = find(~isnan(Tmin_o)==1);
[PP,S]=polyfit(N(I), Tmin_o(I), 1);
Tmin_o_fit = PP(2)+N(I).*PP(1);
slope_Tmin_o = PP(1)*(N(2)-N(1)) ;

% -- Tmin slope modeled_kobs -- %
I = find(~isnan(Tmin_ko)==1);
[PP,S]=polyfit(N(I), Tmin_ko(I), 1);
Tmin_ko_fit = PP(2)+N(I).*PP(1);
slope_Tmin_ko = PP(1)*(N(2)-N(1)) ;

% -- Tmin slope modeled_kcst -- %
I = find(~isnan(Tmin_kc)==1);
[PP,S]=polyfit(N(I), Tmin_kc(I), 1);
Tmin_kc_fit = PP(2)+N(I).*PP(1);
slope_Tmin_kc = PP(1)*(N(2)-N(1)) ;

% -- Thickness observed -- %
I = find(~isnan(thick_o)==1);
[PP,S]=polyfit(N(I), thick_o(I), 1);
thick_o_fit = PP(2)+N(I).*PP(1);
slope_thick_o = PP(1)*(N(2)-N(1)) ;

% -- Thickness modeled_kobs -- %
I = find(~isnan(thick_ko)==1);
[PP,S]=polyfit(N(I), thick_ko(I), 1);
thick_ko_fit = PP(2)+N(I).*PP(1);
slope_thick_ko = PP(1)*(N(2)-N(1)) ;

% -- Thickness modeled_Kcst -- %
I = find(~isnan(thick_kc)==1);
[PP,S]=polyfit(N(I), thick_kc(I), 1);
thick_kc_fit = PP(2)+N(I).*PP(1);
slope_thick_kc = PP(1)*(N(2)-N(1)) ;

% -- heat content observed -- %
I = find(~isnan(heat_o)==1);
[PP,S]=polyfit(N(I), heat_o(I), 1);
heat_o_fit = PP(2)+N.*PP(1);
slope_heat_o = PP(1)*(N(2)-N(1)) ;

% -- heat content modeled_kobs -- %
I = find(~isnan(heat_ko)==1);
[PP,S]=polyfit(N(I), heat_ko(I), 1);
heat_ko_fit = PP(2)+N.*PP(1);
slope_heat_ko = PP(1)*(N(2)-N(1)) ;

% -- heat content modeled_kcst -- %
I = find(~isnan(heat_kc)==1);
[PP,S]=polyfit(N(I), heat_kc(I), 1);
heat_kc_fit = PP(2)+N.*PP(1);
slope_heat_kc = PP(1)*(N(2)-N(1)) ;


%%%%%%%%%%%%%%%%%%%%%%%
% --  plot result  -- %
%%%%%%%%%%%%%%%%%%%%%%%
exy = 0.03;

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 14 16])


% Tmin %
s1=subplot(3,1,1);
% Adjust space to minimize white space
pos1 = get(s1, 'position');
set(s1, 'position', [pos1(1) pos1(2)-exy pos1(3) pos1(4)+exy])

plot(N, Tmin_o, '*k')
hold on
plot(N, Tmin_ko, 'ok')
plot(N, Tmin_kc, 'dk')
legend('observations', 'K(z)', 'K=5\times10^{-5} m^2 s^{-1}', ...
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
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

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
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_thick_o slope_thick_ko slope_thick_kc]


% heat %
s3=subplot(3,1,3);
pos3 = get(s3, 'position');
set(s3, 'position', [pos3(1) pos3(2) pos3(3) pos3(4)+exy])

% if intergrated only
% $$$ plot(N, heat_o./1000000, '*k')
% $$$ hold on
% $$$ plot(N, heat_ko./1000000, 'ok')
% $$$ plot(N, heat_kc./1000000, 'dk')
% $$$ Rplot = plot(N, heat_o_fit./1000000);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
% $$$ Rplot = plot(N, heat_ko_fit./1000000);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
% $$$ Rplot = plot(N, heat_kc_fit./1000000);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
% $$$ hold off
% $$$ set(gca, 'xticklabel', [])
% $$$ xlim(T_LIM);
% $$$ ylab = ylabel('H(GJ m^{-2})', 'fontsize', 10, 'VerticalAlignment', 'top');
% $$$ set(gca, 'xtick', TIK)
% $$$ set(gca, 'fontsize', 10)
% $$$ set(gca, 'XGrid', 'on')
% $$$ set(gca, 'YGrid', 'on')


% mean heat content (divided by d)
plot(N, heat_o./1000, '*k')
hold on
plot(N, heat_ko./1000, 'ok')
plot(N, heat_kc./1000, 'dk')

Rplot = plot(N, heat_o_fit./1000);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 0.25)
Rplot = plot(N, heat_ko_fit./1000);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
Rplot = plot(N, heat_kc_fit./1000);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 0.25)
hold off
set(gca, 'xticklabel', [])
xlim(T_LIM);
ylab = ylabel('H(MJ m^{-3})', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

for i=2:length(month)-1
    mm=datestr(N(i), 3);
    day_norm = (N(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

[slope_heat_o slope_heat_ko slope_heat_kc]


print('-deps2', 'obs_vs_mod.eps')
print('-dpng', '-r300', 'obs_vs_mod.png')

% $$$ 
% $$$ % -- heat content -- %
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ plot(N, HC_m, '.-b')
% $$$ hold on
% $$$ 
% $$$ plot(N, HC_o, '.-r')
% $$$ datetick('x', 3)
% $$$ 
% $$$ Rplot = plot(N, HC_m_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'b','linestyle', '--', 'linewidth', 0.25)
% $$$ 
% $$$ hold on
% $$$ Rplot = plot(N, HC_o_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 0.25)
% $$$ 
% $$$ legend('modeled', 'observed', 'location', 'northwest')
% $$$ ylabel('kj/m^3/month')
% $$$ set(gca, 'xtick', N)
% $$$ datetick('x', 3, 'keepticks')
% $$$ title('Warming rate (Heat content)')
% $$$ 
% $$$ hold off
% $$$ 
% $$$ [slope_HC_m slope_HC_o]
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % -- Tmin of the core -- %
% $$$ figure(2)
% $$$ clf
% $$$ plot(N, Tmin_m, '.-b')
% $$$ datetick('x', 3)
% $$$ hold on
% $$$ 
% $$$ plot(N, Tmin_o, '.-r')
% $$$ datetick('x', 3)
% $$$ 
% $$$ Rplot = plot(N, Tmin_m_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'b','linestyle', '--', 'linewidth', 0.25)
% $$$ 
% $$$ hold on
% $$$ Rplot = plot(N, Tmin_o_fit);%, 'k', 'linewidth', 0.25);
% $$$ set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 0.25)
% $$$ 
% $$$ legend('modeled', 'observed', 'location', 'northwest')
% $$$ ylabel({'^{\circ}C/month'})
% $$$ set(gca, 'xtick', N)
% $$$ datetick('x', 3, 'keepticks')
% $$$ title('Warming rate (core temp.)')
% $$$ 
% $$$ hold off
% $$$ 
% $$$ [slope_Tmin_m slope_Tmin_o]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --  same but colorplot -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exy = 0.03;

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 14 16])


% Tmin %
s1=subplot(3,1,1);
% Adjust space to minimize white space
pos1 = get(s1, 'position');
set(s1, 'position', [pos1(1) pos1(2)-exy pos1(3) pos1(4)+exy])

plot(N, Tmin_o, '.b')
hold on
plot(N, Tmin_ko, '.k')
plot(N, Tmin_kc, '.r')
legend('observations', 'K(z)', 'K=5\times10^{-5} m^2 s^{-1}', ...
       'location', 'northwest')

Rplot = plot(N, Tmin_o_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'b','linestyle', '-', 'linewidth', 1)
Rplot = plot(N, Tmin_ko_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 1)
Rplot = plot(N, Tmin_kc_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 1)
hold off

xlim(T_LIM);
ylim([-1 1.5])
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
ypos1 = ypos(1)-15;
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_Tmin_o slope_Tmin_ko slope_Tmin_kc]


% Thickness %
s2=subplot(3,1,2);
pos2 = get(s2, 'position');
set(s2, 'position', [pos2(1) pos2(2)-exy/2 pos2(3) pos2(4)+exy])

plot(N, thick_o, '.b')
hold on
plot(N, thick_ko, '.k')
plot(N, thick_kc, '.r')

Rplot = plot(N, thick_o_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'b','linestyle', '-', 'linewidth', 1)
Rplot = plot(N, thick_ko_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 1)
Rplot = plot(N, thick_kc_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 1)
hold off

set(gca, 'xticklabel', [])
xlim(T_LIM);
%ylim([-1 1])
ylab = ylabel('d(m)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_thick_o slope_thick_ko slope_thick_kc]


% heat %
s3=subplot(3,1,3);
pos3 = get(s3, 'position');
set(s3, 'position', [pos3(1) pos3(2) pos3(3) pos3(4)+exy])

% if intergrated only
% $$$ plot(N, heat_o./1000000, '*k')
% $$$ hold on
% $$$ plot(N, heat_ko./1000000, 'ok')
% $$$ plot(N, heat_kc./1000000, 'dk')
% $$$ Rplot = plot(N, heat_o_fit./1000000);%, 'k', 'linewidth', 1);
% $$$ set(Rplot, 'color', 'k','linestyle', '-', 'linewidth', 1)
% $$$ Rplot = plot(N, heat_ko_fit./1000000);%, 'k', 'linewidth', 1);
% $$$ set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 1)
% $$$ Rplot = plot(N, heat_kc_fit./1000000);%, 'k', 'linewidth', 1);
% $$$ set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 1)
% $$$ hold off
% $$$ set(gca, 'xticklabel', [])
% $$$ xlim(T_LIM);
% $$$ ylab = ylabel('H(GJ m^{-2})', 'fontsize', 10, 'VerticalAlignment', 'top');
% $$$ set(gca, 'xtick', TIK)
% $$$ set(gca, 'fontsize', 10)
% $$$ set(gca, 'XGrid', 'on')
% $$$ set(gca, 'YGrid', 'on')


% mean heat content (divided by d)
plot(N, heat_o./1000, '.b')
hold on
plot(N, heat_ko./1000, '.k')
plot(N, heat_kc./1000, '.r')

Rplot = plot(N, heat_o_fit./1000);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'b','linestyle', '-', 'linewidth', 1)
Rplot = plot(N, heat_ko_fit./1000);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'k','linestyle', '--', 'linewidth', 1)
Rplot = plot(N, heat_kc_fit./1000);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 1)
hold off
set(gca, 'xticklabel', [])
xlim(T_LIM);
ylab = ylabel('H(MJ m^{-3})', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

for i=2:length(month)-1
    mm=datestr(N(i), 3);
    day_norm = (N(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

[slope_heat_o slope_heat_ko slope_heat_kc]


print('-deps2', 'obs_vs_modc.eps')
print('-dpng', '-r300', 'obs_vs_modc.png')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --  same colorplot but for AGU (no model with k(z))-- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exy = 0.03;

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 14 10])


% Tmin %
s1=subplot(3,1,1);
% Adjust space to minimize white space
pos1 = get(s1, 'position');
set(s1, 'position', [pos1(1) pos1(2)-exy pos1(3) pos1(4)+exy])

plot(N, Tmin_o, '.b')
hold on
plot(N, Tmin_kc, '.r')
legend('observations', 'K=5\times10^{-5} m^2 s^{-1}', ...
        'fontsize', 10, 'location', 'northwest')

Rplot = plot(N, Tmin_o_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'b','linestyle', '-', 'linewidth', 1)
Rplot = plot(N, Tmin_kc_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 1)
hold off

xlim(T_LIM);
ylim([-1 2])
ylab = ylabel('T_{min}(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)

ypos = get(ylab, 'position');
ypos1 = ypos(1)-15;
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_Tmin_o slope_Tmin_ko slope_Tmin_kc]


% Thickness %
s2=subplot(3,1,2);
pos2 = get(s2, 'position');
set(s2, 'position', [pos2(1) pos2(2)-exy/2 pos2(3) pos2(4)+exy])

plot(N, thick_o, '.b')
hold on
plot(N, thick_kc, '.r')

Rplot = plot(N, thick_o_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'b','linestyle', '-', 'linewidth', 1)
Rplot = plot(N, thick_kc_fit);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 1)
hold off

set(gca, 'xticklabel', [])
xlim(T_LIM);
%ylim([-1 1])
ylab = ylabel('d(m)', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

[slope_thick_o slope_thick_ko slope_thick_kc]


% heat %
s3=subplot(3,1,3);
pos3 = get(s3, 'position');
set(s3, 'position', [pos3(1) pos3(2) pos3(3) pos3(4)+exy])

% mean heat content (divided by d)
plot(N, heat_o./1000, '.b')
hold on
plot(N, heat_kc./1000, '.r')

Rplot = plot(N, heat_o_fit./1000);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'b','linestyle', '-', 'linewidth', 1)
Rplot = plot(N, heat_kc_fit./1000);%, 'k', 'linewidth', 1);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 1)
hold off
set(gca, 'xticklabel', [])
xlim(T_LIM);
ylim([6 12])
ylab = ylabel('H(MJ m^{-3})', 'fontsize', 10, 'VerticalAlignment', 'top');
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')

ypos = get(ylab, 'position');
set(ylab, 'position', [ypos1 ypos(2) ypos(3)])

for i=2:length(month)-1
    mm=datestr(N(i), 3);
    day_norm = (N(i)-T_LIM(1))/no_days;
    text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
         'center', 'verticalalignment', 'top', 'fontsize', 10)
end

[slope_heat_o slope_heat_ko slope_heat_kc]


print('-depsc2', 'obs_vs_modc_agu.eps')
print('-dpng', '-r300', 'obs_vs_modc_agu.png')

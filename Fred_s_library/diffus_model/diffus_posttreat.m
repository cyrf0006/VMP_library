clear
% % some constants
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.035e3;%kg/m^3
GAMMA = 0.2; %mixing efficiency
freeze_pt = -1.8; 
cp = 3.99; %Kj/Kg/K
CIL_def = 1;
dz=1;


% $$$ T=load('T_diffus.dat');
% $$$ S=load('S_diffus.dat');
% $$$ N=load('N_diffus.dat');

T=load('T_diffus_daily.dat');
S=load('S_diffus_daily.dat');
N=load('N_diffus_daily.dat');

P = [1:300]';


% monthly average quantities

% Monthly average for which months?
month = unique(str2num(datestr(N,5)));
disp('monthly average from model...')
for i = 1:length(month)
    disp(sprintf('  month %d', month(i)))
    II = find(str2num(datestr(N,5))==month(i));
    climT(:,i) = nanmean(T(:,II),2);
    climS(:,i) = nanmean(S(:,II),2);
    climN(i) = datenum(999, month(i), 15);
end

T = climT;
S = climS;
N = climN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Heat cont and Tmin slope -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(N)

    % -- find heat content modeled -- %
    I = find(T(:,i) <= CIL_def);
        
    if ~isempty(I)==1
        A = min(I); %shallowest limit of CIL
        B = max(I); %deepest limit of CIL
        DENS = sw_dens(S(A:B,i), T(A:B,i), P(A:B)); %density for bin into CIL core

        HC_m(i) = cp*nansum(DENS.*(T(A:B,i)-freeze_pt))*dz/length(A:B)/1000; %in KJ/m3
    else
        HC_m(i)=NaN;
    end
    
    % -- find heat content observed -- %
    month = str2num(datestr(N(i),5));    
    
    if month <10
        Toutname  = sprintf('T_bootclim_0%d.dat', month);
        Soutname  = sprintf('S_bootclim_0%d.dat', month);
    else
        Toutname  = sprintf('T_bootclim_%d.dat', month);
        Soutname  = sprintf('S_bootclim_%d.dat', month);
    end
    
    TT = load(Toutname);
    SS = load(Soutname);

    I = find(TT(:,2) <= CIL_def);
        
    if ~isempty(I)==1
        A = min(I); %shallowest limit of CIL
        B = max(I); %deepest limit of CIL
        DENS = sw_dens(SS(A:B,2), TT(A:B,2), P(A:B)); %density for bin into CIL core

        HC_o(i) = cp*nansum(DENS.*(TT(A:B,2)-freeze_pt))*dz/length(A:B)/1000; %in KJ/m3
    else
        HC_o(i)=NaN;
    end    
    
   % -- find Tmin slope modeled -- %
   Tmin_m(i) = min(T(:,i)); 
    
   % -- find Tmin slope observed -- %
   Tmin_o(i) = min(TT(:,2)); 
    
end


    
%%%%%%%%%%%%%%%%%%%%%%%%%
% -- find linear fit -- %
%%%%%%%%%%%%%%%%%%%%%%%%%

% -- heat content modeled -- %
[PP,S]=polyfit(N, HC_m, 1);
HC_m_fit = PP(2)+N.*PP(1);
slope_HC_m = PP(1)*(N(2)-N(1));

% -- heat content observed -- %
[PP,S]=polyfit(N, HC_o, 1);
HC_o_fit = PP(2)+N.*PP(1);
slope_HC_o = PP(1)*(N(2)-N(1));

% -- Tmin slope modeled -- %
[PP,S]=polyfit(N, Tmin_m, 1);
Tmin_m_fit = PP(2)+N.*PP(1);
slope_Tmin_m = PP(1)*(N(2)-N(1));

% -- Tmin slope observed -- %
[PP,S]=polyfit(N, Tmin_o, 1);
Tmin_o_fit = PP(2)+N.*PP(1);
slope_Tmin_o = PP(1)*(N(2)-N(1));


%%%%%%%%%%%%%%%%%%%%%%%
% --  plot result  -- %
%%%%%%%%%%%%%%%%%%%%%%%


% -- heat content -- %

figure(1)
clf
plot(N, HC_m, '.-b')
hold on

plot(N, HC_o, '.-r')
datetick('x', 3)

Rplot = plot(N, HC_m_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'b','linestyle', '--', 'linewidth', 0.25)

hold on
Rplot = plot(N, HC_o_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 0.25)

legend('modeled', 'observed', 'location', 'northwest')
ylabel('kj/m^3/month')
set(gca, 'xtick', N)
datetick('x', 3, 'keepticks')
title('Warming rate (Heat content)')

hold off

[slope_HC_m slope_HC_o]



% -- Tmin of the core -- %
figure(2)
clf
plot(N, Tmin_m, '.-b')
datetick('x', 3)
hold on

plot(N, Tmin_o, '.-r')
datetick('x', 3)

Rplot = plot(N, Tmin_m_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'b','linestyle', '--', 'linewidth', 0.25)

hold on
Rplot = plot(N, Tmin_o_fit);%, 'k', 'linewidth', 0.25);
set(Rplot, 'color', 'r','linestyle', '--', 'linewidth', 0.25)

legend('modeled', 'observed', 'location', 'northwest')
ylabel({'^{\circ}C/month'})
set(gca, 'xtick', N)
datetick('x', 3, 'keepticks')
title('Warming rate (core temp.)')

hold off

[slope_Tmin_m slope_Tmin_o]


% saving Tcore
%dlmwrite('T_core_mod.dat', [N' Tmin_m'],'delimiter',' ','precision',6)
%dlmwrite('T_core_5e-5.dat', [N' Tmin_m'],'delimiter',' ','precision',6)
%dlmwrite('T_core_obs.dat', [N' Tmin_o'],'delimiter',' ','precision',6)


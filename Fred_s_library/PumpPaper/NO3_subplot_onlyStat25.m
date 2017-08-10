clear

% Few infos on plots
zbin = 5;
zMax = 300;
zVec = [zbin/2:zbin:zMax]';



% datasets
VMP_riki = load('~/WINDEX/data_processing/boot_VMP_riki.mat');
VMP_hlc = load('~/WINDEX/data_processing/Mission_tadoussac_2009/boot_VMP_Stat25.mat');

Stat25 = load(['/home/cyrf0006/PhD/Nitrates/Stat25_20km/' ...
               'smooth_scatter_stat25.mat']);
Stat23 = load(['/home/cyrf0006/PhD/Nitrates/riki/' ...
               'smooth_scatter_stat23.mat']);

% bin and rename variables and compute fluxes
epsRiki = nan(size(zVec)); epsRiki1 = nan(size(zVec)); epsRiki2 = nan(size(zVec));
epsHlc = nan(size(zVec)); epsHlc1 = nan(size(zVec)); epsHlc2 = nan(size(zVec));
KRiki = nan(size(zVec)); KRiki1 = nan(size(zVec)); KRiki2 = nan(size(zVec));
KHlc = nan(size(zVec)); KHlc1 = nan(size(zVec)); KHlc2 = nan(size(zVec));
N2Riki = nan(size(zVec)); N2Riki1 = nan(size(zVec)); N2Riki2 = nan(size(zVec));
N2Hlc = nan(size(zVec)); N2Hlc1 = nan(size(zVec)); N2Hlc2 = nan(size(zVec));

NO3Riki = nan(size(zVec)); NO3Riki1 = nan(size(zVec)); NO3Riki2 = nan(size(zVec));
NO3Hlc = nan(size(zVec)); NO3Hlc1 = nan(size(zVec)); NO3Hlc2 = nan(size(zVec));

for i = 1:length(zVec)
    I = find(VMP_riki.P_bin>=zVec(i)-zbin/2 & VMP_riki.P_bin<zVec(i)+zbin/2);
    epsRiki(i) = nanmean(VMP_riki.eps_ave(I));
    epsRiki1(i) = nanmean(VMP_riki.eps_2p5(I));
    epsRiki2(i) = nanmean(VMP_riki.eps_97p5(I));
    KRiki(i) = nanmean(VMP_riki.K_ave(I));
    KRiki1(i) = nanmean(VMP_riki.K_2p5(I));
    KRiki2(i) = nanmean(VMP_riki.K_97p5(I));   
    N2Riki(i) = nanmean(VMP_riki.N2_ave(I));
    N2Riki1(i) = nanmean(VMP_riki.N2_2p5(I));
    N2Riki2(i) = nanmean(VMP_riki.N2_97p5(I));  
    I = find(Stat23.Pfit>=zVec(i)-zbin/2 & Stat23.Pfit<zVec(i)+zbin/2);
    NO3Riki(i) = nanmean(Stat23.Sfit_ave(I));
    NO3Riki1(i) = nanmean(Stat23.Sfit2p5(I)); 
    NO3Riki2(i) = nanmean(Stat23.Sfit97p5(I));
        
    I = find(VMP_hlc.P_bin>=zVec(i)-zbin/2 & VMP_hlc.P_bin<zVec(i)+zbin/2);
    epsHlc(i) = nanmean(VMP_hlc.eps_ave(I));
    epsHlc1(i) = nanmean(VMP_hlc.eps_2p5(I));
    epsHlc2(i) = nanmean(VMP_hlc.eps_97p5(I));
    KHlc(i) = nanmean(VMP_hlc.K_ave(I));
    KHlc1(i) = nanmean(VMP_hlc.K_2p5(I));
    KHlc2(i) = nanmean(VMP_hlc.K_97p5(I));   
    N2Hlc(i) = nanmean(VMP_hlc.N2_ave(I));
    N2Hlc1(i) = nanmean(VMP_hlc.N2_2p5(I));
    N2Hlc2(i) = nanmean(VMP_hlc.N2_97p5(I));  
    I = find(Stat25.Pfit>=zVec(i)-zbin/2 & Stat25.Pfit<zVec(i)+zbin/2);
    NO3Hlc(i) = nanmean(Stat25.Sfit_ave(I));
    NO3Hlc1(i) = nanmean(Stat25.Sfit2p5(I)); 
    NO3Hlc2(i) = nanmean(Stat25.Sfit97p5(I));
end

dNO3dzHlc = gradient(NO3Hlc, zVec);
dNO3dzHlc1 = gradient(NO3Hlc1, zVec);
dNO3dzHlc2 = gradient(NO3Hlc2, zVec);
dNO3dzRiki = gradient(NO3Riki, zVec);
dNO3dzRiki1 = gradient(NO3Riki1, zVec);
dNO3dzRiki2 = gradient(NO3Riki2, zVec);

FNO3Riki = KRiki.*dNO3dzRiki*86400;
FNO3Riki1 = KRiki1.*dNO3dzRiki1*86400;
FNO3Riki2 = KRiki2.*dNO3dzRiki2*86400;
FNO3Hlc = KHlc.*dNO3dzHlc*86400;
FNO3Hlc1 = KHlc1.*dNO3dzHlc1*86400;
FNO3Hlc2 = KHlc2.*dNO3dzHlc2*86400;



Stat25 = load(['/home/cyrf0006/PhD/Nitrates/Stat25_20km/' ...
               'smooth_scatter_stat25.mat']);
Stat23 = load(['/home/cyrf0006/PhD/Nitrates/riki/' ...
               'smooth_scatter_stat23.mat']);




figure(1)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 12])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 3; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.08; % very left of figure
rigs = 0.04; % very right of figure
tops = 0.02; % top of figure
bots = 0.15; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

s8 = subplot(231);
x1 = KHlc1;
x2 = KHlc2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(KHlc(I), zVec(I), 'color', [1 1 1]*.6)
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([1e-6 1e-1])
set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
ylabel('Depth (m)')
xlabel('K (m^{2} s^{-1})')
adjust_space


subplot(232)
plot(Stat25.NO3RawVec, Stat25.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
hold on

% tanh fit on bootstrapped bins
x1 = Stat25.Sfit2p5;
x2 = Stat25.Sfit97p5;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [Stat25.Pfit(I); flipud(Stat25.Pfit(I)); Stat25.Pfit(I(1))], [1 1 1]*.6, 'edgecolor', 'none');

plot(Stat25.NO3RawVec, Stat25.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
I = find(~isnan(Stat25.NO3_boot_dot)==1);
herrorbar(Stat25.NO3_ave(I), Stat25.P(I), Stat25.NO3_ave(I)-Stat25.NO3_2p5(I), Stat25.NO3_97p5(I)-Stat25.NO3_ave(I), '.k');

% ---- MS-2012 Nitrates bottles ---- %
disp('Add bottles from MS2012')
data1 = load('~/PhD/Nitrates/CTD_station_fixe_MS_2012/NO3_SFA.dat');
data2 = load('~/PhD/Nitrates/CTD_station_fixe_MS_2012/NO3_SFB.dat');
% Average all bottles at same station (uncommoent below to check)
I = find(diff(data1(:,1))==0 & diff(data1(:,2))==0); 
while ~isempty(I)
    data1(I(1),:) = nanmean(data1(I:I+1,:), 1);
    data1(I(1)+1, :) = [];
    %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
    I = find(diff(data1(:,1))==0 & diff(data1(:,2))==0);     
end
I = find(diff(data2(:,1))==0 & diff(data2(:,2))==0); 
while ~isempty(I)
    data2(I(1),:) = nanmean(data2(I:I+1,:), 1);
    data2(I(1)+1, :) = [];
    %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
    I = find(diff(data2(:,1))==0 & diff(data2(:,2))==0);     
end
data = [data1; data2];

% vectors
zVec2 = data(:,2);
NO3Vec2 = data(:,4);

for i = 1:length(zVec)
    plot(NO3Vec2, zVec2, 'color',[.3 .3 .3], 'linestyle', '.')
end
% --------------------------------- %
    
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
hold off
ylim([0 zMax])
xlim([6 25])
hold off
xlabel('[NO_3] (mmol m^{-3})')
adjust_space


% -- S10 -- %
s10 = subplot(2,3,3);
x1 = FNO3Hlc1;
x2 = FNO3Hlc2;
I = find(x1>0 & x2>0);
semilogx(FNO3Hlc(I), zVec(I), 'color', [1 1 1]*.6)
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([1e-4 5e2])
set(gca, 'xtick', [1e-4 1e-2 1e0 1e2])
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
xlabel('F_{NO_3} (mmol m^{-2} d^{-1})')
adjust_space





% $$$ 
% $$$ 
% $$$ % -- S3 -- %
% $$$ s3 = subplot(234);
% $$$ x1 = KRiki1;
% $$$ x2 = KRiki2;
% $$$ I = find(~isnan(x1)==1 & ~isnan(x2)==1);
% $$$ semilogx(KRiki(I), zVec(I), 'color', [1 1 1]*.6)
% $$$ hold on
% $$$ patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
% $$$ hold off
% $$$ ylim([0 zMax])
% $$$ xlim([1e-6 1e-1])
% $$$ set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
% $$$ set(gca, 'ydir', 'reverse')
% $$$ set(gca, 'xgrid', 'on')
% $$$ set(gca, 'xminorgrid', 'off')
% $$$ set(gca, 'tickdir', 'out')
% $$$ xlabel('K (m^{2} s^{-1})')
% $$$ set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
% $$$ adjust_space
% $$$ 
% $$$ 
% $$$ % -- S4 -- %
% $$$ s4 = subplot(2,3,5);
% $$$ 
% $$$ plot(Stat23.NO3RawVec, Stat23.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
% $$$ hold on
% $$$ 
% $$$ % tanh fit on bootstrapped bins
% $$$ x1 = Stat23.Sfit2p5;
% $$$ x2 = Stat23.Sfit97p5;
% $$$ I = find(~isnan(x1)==1 & ~isnan(x2)==1);
% $$$ patch([x1(I); flipud(x2(I)); x1(I(1))], [Stat23.Pfit(I); flipud(Stat23.Pfit(I)); Stat23.Pfit(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
% $$$ 
% $$$ plot(Stat23.NO3RawVec, Stat23.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
% $$$ I = find(~isnan(Stat23.NO3_boot_dot)==1);
% $$$ herrorbar(Stat23.NO3_ave(I), Stat23.P(I), Stat23.NO3_ave(I)- ...
% $$$           Stat23.NO3_2p5(I), Stat23.NO3_97p5(I)-Stat23.NO3_ave(I),'.k');
% $$$ 
% $$$ ylim([0 zMax])
% $$$ xlim([6 25])
% $$$ set(gca, 'yticklabel', [])
% $$$ set(gca, 'ydir', 'reverse')
% $$$ set(gca, 'xgrid', 'on')
% $$$ set(gca, 'xminorgrid', 'off')
% $$$ set(gca, 'tickdir', 'out')
% $$$ set(gca, 'ytick', [0:50:300])
% $$$ xlabel('[NO_3] (mmol m^{-3})')
% $$$ hold off
% $$$ hold off
% $$$ adjust_space
% $$$ 
% $$$ 
% $$$ % -- S5 -- %
% $$$ s5 = subplot(236);
% $$$ x1 = FNO3Riki1;
% $$$ x2 = FNO3Riki2;
% $$$ 
% $$$ I = find(x1>0 & x2>0);
% $$$ semilogx(FNO3Riki(I), zVec(I), 'color', [1 1 1]*.6)
% $$$ hold on
% $$$ patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
% $$$ %I = find(x1<0 & x2<0);
% $$$ %patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1], 'edgecolor', 'none');
% $$$ hold off
% $$$ ylim([0 zMax])
% $$$ xlim([1e-4 5e2])
% $$$ set(gca, 'xtick', [1e-4 1e-2 1e0 1e2])
% $$$ xlabel('F_{NO_3} (mmol m^{-2} d^{-1})')
% $$$ %set(gca, 'xminortick', [])
% $$$ set(gca, 'yticklabel', [])
% $$$ set(gca, 'ydir', 'reverse')
% $$$ set(gca, 'xgrid', 'on')
% $$$ set(gca, 'xminorgrid', 'off')
% $$$ % $$$ set(gca, 'ygrid', 'on')
% $$$ % $$$ set(gca, 'yminorgrid', 'off')
% $$$ set(gca, 'tickdir', 'out')
% $$$ adjust_space
% $$$ 







print('-dpng', '-r300', 'NO3_subplot_only25.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'NO3_subplot_only25.eps')

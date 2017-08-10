clear

% Few infos on plots
zbin = 5;
zMax = 300;
zVec = [zbin/2:zbin:zMax]';



% datasets
VMP_riki = load('~/WINDEX/data_processing/boot_VMP_riki.mat');
VMP_hlc = load('~/WINDEX/data_processing/Mission_tadoussac_2009/boot_VMP_Stat25.mat');
NO3_stat23 = load('~/PhD/Nitrates/Stat25_20km/boot_NO3_stat25.dat');
NO3_stat25 = load('~/PhD/Nitrates/riki/boot_NO3_stat23.dat');

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
    I = find(NO3_stat23(:,1)>=zVec(i)-zbin/2 & NO3_stat23(:,1)<zVec(i)+zbin/2);
    NO3Riki(i) = nanmean(NO3_stat23(I,2));
    NO3Riki1(i) = nanmean(NO3_stat23(I,3));
    NO3Riki2(i) = nanmean(NO3_stat23(I,4));
    
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
    I = find(NO3_stat25(:,1)>=zVec(i)-zbin/2 & NO3_stat25(:,1)<zVec(i)+zbin/2);
    NO3Hlc(i) = nanmean(NO3_stat25(I,2));
    NO3Hlc1(i) = nanmean(NO3_stat25(I,3));
    NO3Hlc2(i) = nanmean(NO3_stat25(I,4));
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


figure(1)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 20])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 5; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.02; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %



% --  S1 -- %
s1 = subplot(251);
x1 = N2Riki1;
x2 = N2Riki2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(N2Riki(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
%text(4e-4, 172, 'Stat. 23')
text(4e-4, 290, 'Stat. 23')
hold off
ylabel('Depth (m)')
ylim([0 zMax])
xlim([4e-5 5e-3])
set(gca, 'ydir', 'reverse')
set(gca, 'xtick', [1e-4 1e-3])
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
%et(s1, ')
adjust_space

% -- S2 -- %
s2 = subplot(252);
x1 = epsRiki1;
x2 = epsRiki2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(epsRiki(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([1e-9 1e-3])
set(gca, 'xtick', [1e-8 1e-7 1e-6 1e-5 1e-4])
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space

% -- S3 -- %
s3 = subplot(253);
x1 = KRiki1;
x2 = KRiki2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(KRiki(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([1e-6 1e-1])
set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


% -- S4 -- %
s4 = subplot(254);
x1 = NO3Riki1;
x2 = NO3Riki2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
plot(NO3Riki(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([6 25])
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


% -- S5 -- %
s5 = subplot(255);
x1 = FNO3Riki1;
x2 = FNO3Riki2;

I = find(x1>0 & x2>0);
semilogx(FNO3Riki(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
%I = find(x1<0 & x2<0);
%patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1], 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([5e-3 5e2])
set(gca, 'xtick', [1e-2 1e-1 1e0 1e1 1e2])
%set(gca, 'xminortick', [])
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
% $$$ set(gca, 'ygrid', 'on')
% $$$ set(gca, 'yminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space





% --  S6-- %
s6 = subplot(256);
x1 = N2Hlc1;
x2 = N2Hlc2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(N2Hlc(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
%text(4e-4, 172, 'Stat. 25')
text(4e-4, 290, 'Stat. 25')
hold off
ylim([0 zMax])
xlim([4e-5 5e-3])
xlabel('N^2 (s^{-2})')
set(gca, 'xtick', [1e-4 1e-3])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


% -- S7 -- %
s7 = subplot(257);
x1 = epsHlc1;
x2 = epsHlc2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(epsHlc(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([1e-9 1e-3])
xlabel('\epsilon (W kg^{-1})')
set(gca, 'xtick', [1e-8 1e-7 1e-6 1e-5 1e-4])
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space

% -- S8 -- %
s8 = subplot(258);
x1 = KHlc1;
x2 = KHlc2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
semilogx(KHlc(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([1e-6 1e-1])
xlabel('K (m^{2} s^{-1})')
set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


% -- S9 -- %
s9 = subplot(259);
x1 = NO3Hlc1;
x2 = NO3Hlc2;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
plot(NO3Hlc(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([6 25])
xlabel('[NO_3] (mmol m^{-3})')
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


% -- S10 -- %
s10 = subplot(2,5,10);
x1 = FNO3Hlc1;
x2 = FNO3Hlc2;
I = find(x1>0 & x2>0);
semilogx(FNO3Hlc(I), zVec(I), 'k')
hold on
patch([x1(I); flipud(x2(I)); x1(I(1))], [zVec(I); flipud(zVec(I)); zVec(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
hold off
ylim([0 zMax])
xlim([5e-3 5e2])
xlabel('F_{NO_3} (mmol m^{-2} d^{-1})')
set(gca, 'xtick', [1e-2 1e-1 1e0 1e1 1e2])
set(gca, 'yticklabel', [])
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


print(gcf, '-dpng', 'NO3_subplot.png')
set(gcf, 'renderer', 'painters')
print(gcf, '-depsc2', 'NO3_subplot.eps')

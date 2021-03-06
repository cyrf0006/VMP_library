% Figure for tanh fit on NOx concentration at Stat. 23 and 25
clear   
   
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 11 13])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.02; % very right of figure
tops = 0.03; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

% these data from: boot_NO3_from_S.m
Stat25 = load(['/home/cyrf0006/PhD/Nitrates/Stat25_20km/' ...
               'smooth_scatter_stat25.mat']);
Stat23 = load(['/home/cyrf0006/PhD/Nitrates/riki/' ...
               'smooth_scatter_stat23.mat']);

VMP_riki = load('~/WINDEX_X1/data_processing/boot_VMP_riki.mat');
VMP_hlc = load('~/WINDEX_X1/data_processing/Mission_tadoussac_2009/boot_VMP_Stat25.mat');

zNO3 = VMP_riki.P_bin;
NO3riki = VMP_riki.NO3_ave;
NO3riki97p5 = VMP_riki.NO3_97p5;
NO3riki2p5 = VMP_riki.NO3_2p5;
NO3hcl97p5 = VMP_hlc.NO3_97p5;
NO3hlc2p5 = VMP_hlc.NO3_2p5;


%% S1 %%
subplot(121)
plot(Stat25.NO3RawVec, Stat25.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
hold on

% tanh fit on bootstrapped bins
% $$$ x1 = Stat25.Sfit2p5;
% $$$ x2 = Stat25.Sfit97p5;
% $$$ zfit = Stat25.Pfit;
% From S bootstrapped in boot_NO3...
x1 = NO3hlc2p5;
x2 = NO3hcl97p5;
zfit = zNO3;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [zfit(I); flipud(zfit(I)); zfit(I(1))], [1 1 1]*.6, 'edgecolor', 'none');

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
zVec = data(:,2);
PO4Vec = data(:,3);
NO3Vec = data(:,4);
SiVec = data(:,5);   

% save for further use in scatterNO3.m
save NO3FromMS2012.mat zVec PO4Vec NO3Vec SiVec

for i = 1:length(zVec)
    plot(NO3Vec, zVec, '.m')
end
% error bar for these profiles only
Pbin = Stat25.P;
dp = Pbin(2)-Pbin(1);
CORave = nan(length(Pbin),1);
CORmin = nan(length(Pbin),1);
CORmax = nan(length(Pbin),1);
for i = 1:length(Pbin)
    I = find(zVec>=Pbin(i)-dp/2 & zVec<Pbin(i)+dp/2);
    if ~isempty(I)
        CORave(i) = nanmean(NO3Vec(I));
        CORmin(i) = nanmin(NO3Vec(I));
        CORmax(i) = nanmax(NO3Vec(I));
    end
end
I = find(~isnan(CORave)==1);
%herrorbar(CORave(I), Pbin(I), CORave(I)-CORmin(I), CORmax(I)-CORave(I), '.m');
%plot(CORave, Pbin, 'm', 'linewidth', 2)

% --------------------------------- %
title('Station 25', 'fontWeight', 'bold') 
set(gca, 'ydir', 'reverse')
set(gca, 'xgrid', 'on')
hold off
xlab = xlabel('C_{NO_3} (mmol m^{-3})');
xlab_pos = get(xlab, 'pos');
xlab_pos(2) = xlab_pos(2)-1;
set(xlab, 'pos', xlab_pos)
ylabel('Depth(m)')
ylim([0 350])
xlim([0 28])
hold off
adjust_space


%% S2 %%
subplot(122)
plot(Stat23.NO3RawVec, Stat23.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
hold on

% tanh fit on bootstrapped bins
% $$$ x1 = Stat23.Sfit2p5;
% $$$ x2 = Stat23.Sfit97p5;
% $$$ zfit = Stat23.Pfit;
% From S bootstrapped in boot_NO3...
x1 = NO3riki2p5;
x2 = NO3riki97p5;
zfit = zNO3;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [zfit(I); flipud(zfit(I)); zfit(I(1))], [1 1 1]*.6, 'edgecolor', 'none');

plot(Stat23.NO3RawVec, Stat23.zRawVec, 'color',[.3 .3 .3], 'linestyle', '.')
I = find(~isnan(Stat23.NO3_boot_dot)==1);
herrorbar(Stat23.NO3_ave(I), Stat23.P(I), Stat23.NO3_ave(I)- ...
          Stat23.NO3_2p5(I), Stat23.NO3_97p5(I)-Stat23.NO3_ave(I),'.k');
title('Station 23', 'fontWeight', 'bold')
set(gca, 'ydir', 'reverse')
set(gca, 'yticklabel', '')
xlab = xlabel('C_{NO_3} (mmol m^{-3})');
xlab_pos = get(xlab, 'pos');
xlab_pos(2) = xlab_pos(2)-1;
set(xlab, 'pos', xlab_pos)
set(gca, 'xgrid', 'on')
hold off
ylim([0 350])
xlim([0 28])
hold off
adjust_space

set(gcf, 'renderer', 'painters')
print('-depsc2', 'CNOx.eps')

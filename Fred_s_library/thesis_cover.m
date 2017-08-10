clear
% Which wave
iWave = 3;

fname(:,1) = '/home/cyrf0006/LaTeX/MS/BMix/matlab_files/plot_ISW/windex20120808_212335.mat';
fname(:,2) = '/home/cyrf0006/LaTeX/MS/BMix/matlab_files/plot_ISW/windex20121025_135931.mat';
fname(:,3) = '/home/cyrf0006/LaTeX/MS/BMix/matlab_files/plot_ISW/windex20121025_144730.mat';

iMin(1) = 7000;
iMax(1) = 9722;
kMin(1) = 1;
kMax(1) = 2600;

iMin(2) = 1600;
iMax(2) = 3600;
kMin(2) = 1;
kMax(2) = 6500;

% $$$ iMin(3) = 1;
% $$$ iMax(3) = 2000;
% $$$ kMin(3) = 1;
% $$$ kMax(3) = 5000;
iMin(3) = 1;
iMax(3) = 2000;
kMin(3) = 1;
kMax(3) = 5000;

% z index offset to find the bottom without finding 
% mid-water colum maximum (e.g.: VMP trace).
kOffset(1) = 2000;
kOffset(2) = 5000;
kOffset(3) = 4000;


% caxis limit
cmin = 0;
cmax = 4;

% Fontsize
fs = 10;

% FIGURE1 (2 subplot)
lengthFig = 20;
heightFig = 15;
figure(1);
set(gcf,'PaperUnits','centimeters','PaperSize',[lengthFig heightFig]);
set(gcf,'PaperPosition',[0 0 lengthFig heightFig]);
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.06 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.02; % very right of figure
tops = 0.02; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
       
% Load data 
load(fname(:,iWave));
    
% Reduce the data
int   = int';
int   = int(iMin(iWave):iMax(iWave),kMin(iWave):kMax(iWave));
%int   = int';
mtime = mtime(iMin(iWave):iMax(iWave));
z     = z(kMin(iWave):kMax(iWave));

% Find bottom and replace with NaNs;
for i = 1:length(mtime)
    kbot = find(int(i,kOffset(iWave):end) == max(int(i,kOffset(iWave):end)));
    kbot = kbot(1) + kOffset(iWave);   
    int(i,kbot:end) = NaN;
end

xticks = floor(mtime(1)*1440/5)*5/1440:5/1440:ceil(mtime(end)*1440/5)*5/1440;

colormap(flipud(gray));
    
imagesc(mtime,z,int');
caxis([cmin cmax])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca,'yminortick','off','xminortick','off');
xlim([mtime(1) mtime(end)])
adjust_space

ylim([2.5 40])
xlim([datenum(2012, 10, 25, 14, 50, 0), datenum(2012, 10, 25, 15, 5, 0)])
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
caxis([1 3])
set(gca, 'outer', [0 0 1 1])
set(gca, 'pos', [0 0 1 1])
%set(gca,'FontSize',12,'FontWeight','bold','linewidth',5)
%set(gca,'box', 'on','linewidth',5)
set(gca,'box', 'off')
set(gcf, 'renderer', 'painters')
print('-deps2', 'thesis_cover.eps');
print('-dpdf', 'thesis_cover.pdf');


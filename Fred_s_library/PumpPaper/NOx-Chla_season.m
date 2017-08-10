function NOx-Chla_season(datFiles, zMin, zMax, varargin)
    
% function NOx_relationship(datFiles, zMin, zMax, varargin)
%
% Will build a scatterplot of any of the following variables given
% as input: 'P', 'T', 'S', NO2', 'NOx', 'PO4', 'SIO4', 'O2' in any
% order. zMin and zMax will be ignored whe trying to obtain a
% linear best fit.
%
% usage ex:
% NOx_season('datFiles.list', 0, 25, 'NOx', 'S') 
% to be run in /home/cyrf0006/PhD/Nitrates/allStat 
%
% Or
%
% NOx_season('LSLE_stats.list', 0, 25, 'NOx', 'S') 
% to be run in /home/cyrf0006/PhD/Nitrates/StatsLSLE
%
% F. Cyr - September 2013


% Few params that could be passed as function inputs
yearmin = 2000;
%zMin = 25;
%zMax = 300;

% Varargin test
if size(varargin,2)~=2
    disp('Default relationship, NOx vs S')
    varNOx = 1;
    varS = 1;
    varNO2 = 0; varT = 0; varO2 = 0; varPO4 = 0; varSIO4 = 0;
else
    varNOx = strcmp(varargin{1}, 'NOx') + strcmp(varargin{2}, 'NOx');
    varNO2 = strcmp(varargin{1}, 'NO2') + strcmp(varargin{2}, 'NO2');
    varO2 = strcmp(varargin{1}, 'O2') + strcmp(varargin{2}, 'NO2');
    varT = strcmp(varargin{1}, 'T') + strcmp(varargin{2}, 'T');
    varS = strcmp(varargin{1}, 'S') + strcmp(varargin{2}, 'S');
    varPO4 = strcmp(varargin{1}, 'PO4') + strcmp(varargin{2}, 'PO4');
    varSIO4 = strcmp(varargin{1}, 'SIO4') + strcmp(varargin{2}, 'SIO4');
end


% Columns order in .dat files (from ~/shellscripts/ODF2ASCII_bottle_NO3-O2rel)
colOrder{1} = 'P';
colOrder{2} = 'T';
colOrder{3} = 'S';
colOrder{4} = 'NO2';
colOrder{5} = 'NOx';
colOrder{6} = 'PO4';
colOrder{7} = 'SIO4';
colOrder{8} = 'O2';

% variable order for plotting (1st variable will be x, 2nd y)
varOrder{1} = 'S';
varOrder{2} = 'T';
varOrder{3} = 'NO2';
varOrder{4} = 'NOx';
varOrder{5} = 'PO4';
varOrder{6} = 'SIO4';
varOrder{7} = 'O2';
varOrder{8} = 'P';

% which variables
if find(strcmp(varOrder, varargin{1})) < find(strcmp(varOrder, varargin{2}))
    XVAR_COL = find(strcmp(colOrder, varargin{1})==1);
    YVAR_COL = find(strcmp(colOrder, varargin{2})==1);
else
    XVAR_COL = find(strcmp(colOrder, varargin{2})==1);
    YVAR_COL = find(strcmp(colOrder, varargin{1})==1);
end

% load file names
fid = fopen(datFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

n = datenum(ctd(:,1:11)); % in files from ODF2ASCII_bottle.sh


% year restriction
I = find(str2num(datestr(n, 10))<yearmin);
n(I) = [];
ctd(I,:) = [];


%% Store and clean data
noFiles = size(ctd,1); 
zVec = [];
xVec = [];
yVec = [];
timeVec = [];
for j = 1:length(n)       
    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=length(colOrder)
            disp('Problem with file size')
            keyboard
        end
            
        pres = file(:,1);
        xVar = file(:,XVAR_COL);
        yVar = file(:,YVAR_COL);

        % store variables
        zVec = [zVec; pres];
        xVec = [xVec; xVar];
        yVec = [yVec; yVar];                
        timeVec = [timeVec; pres*0+n(j)];
    end
    
end  % for j

% Quality control
I = find(xVec==-99);
xVec(I) = [];
yVec(I) = [];
zVec(I) = [];
timeVec(I) = [];

I = find(yVec==-99);
xVec(I) = [];
yVec(I) = [];
zVec(I) = [];
timeVec(I) = [];

I = find(zVec==-99);
xVec(I) = [];
yVec(I) = [];
zVec(I) = [];
timeVec(I) = [];

% Surface pts
xVecSurf = xVec;
yVecSurf = yVec;
zVecSurf = zVec;
timeVecSurf = timeVec;
I = find(zVec<zMin);
xVecSurf(I) = [];
yVecSurf(I) = [];
zVecSurf(I) = [];
timeVecSurf(I) = [];

I = find(zVecSurf>zMax);
xVecSurf(I) = [];
yVecSurf(I) = [];
zVecSurf(I) = [];
timeVecSurf(I) = [];

% Deep pts
xVecDeep = xVec;
yVecDeep = yVec;
zVecDeep = zVec;
timeVecDeep = timeVec;
I = find(xVecDeep<32);
xVecDeep(I) = [];
yVecDeep(I) = [];
zVecDeep(I) = [];
timeVecDeep(I) = [];

I = find(xVecDeep>34.3);
xVecDeep(I) = [];
yVecDeep(I) = [];
zVecDeep(I) = [];
timeVecDeep(I) = [];

noFiles = size(ctd,1); 

xxVec = [];
yyVec = [];

% Build timesereis
monthVec = 4:11;
aveSurf = nan(size(monthVec));
aveDeep = nan(size(monthVec));
stdSurf = nan(size(monthVec));
stdDeep = nan(size(monthVec));
monthVecSurf = str2num(datestr(timeVecSurf, 5));
monthVecDeep = str2num(datestr(timeVecDeep, 5));
for i = 1:length(monthVec)
    I = find(monthVecSurf == monthVec(i));
    aveSurf(i) = nanmean(xVecSurf(I));
    stdSurf(i) = nanstd(xVecSurf(I));
    I = find(monthVecDeep == monthVec(i));
    aveDeep(i) = nanmean(xVecDeep(I));
    stdDeep(i) = nanstd(xVecDeep(I));
end

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 10])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
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
% $$$ disp(sprintf('%d pts plotted from %d files', length(zVec), length(n)))
% $$$ disp(sprintf([' -> falls to %d pts when ignoring bottles shallower than %dm'], length(zVec2), zMin))
% $$$ 
% NO2 plot

errorbar(monthVec, aveSurf, stdSurf, '.', 'color', [1 1 1]*.5, 'linewidth', 2, ...
         'markerSize', 15) 
hold on
errorbar(monthVec, aveDeep, stdDeep, '.k', 'linewidth', 2, 'markerSize', 15)
legend('0-25m','32-34.3 g kg^{-1}', 'location', 'southEast')
ylabel('C_{NO_3}', 'fontSize', 12)
xlabel('Month of year', 'fontSize', 12)
xlim([3.5 11.5])
ylim([24 35])
set(gca, 'ygrid', 'on')
set(gca, 'fontSize', 12)


print(gcf, '-dpng', 'NOxAnnualCycle.png')
set(gcf, 'renderer', 'painters')
print(gcf, '-depsc2', 'NOxAnnualCycle.eps')





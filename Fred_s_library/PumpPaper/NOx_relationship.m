function NOx_relationship(datFiles, zMin, zMax, varargin)
    
% function NOx_relationship(datFiles, zMin, zMax, varargin)
%
% Will build a scatterplot of any of the following variables given
% as input: 'P', 'T', 'S', NO2', 'NOx', 'PO4', 'SIO4', 'O2' in any
% order. zMin and zMax will be ignored whe trying to obtain a
% linear best fit.
%
% usage ex:
% NOx_relationship('datFiles.list', 25, 300, 'NOx', 'O2')
% NOx_relationship('datFiles.list', 50, 250, 'NOx', 'S') 
% NOx_relationship('datFiles.list', 50, 300, 'NOx', 'S') 
% NOx_relationship('datFiles.list', 60, 250, 'NOx', 'S') < --- Zze one
% to be run in /home/cyrf0006/PhD/Nitrates/allStat 
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

% month restriction
I = find(str2num(datestr(n, 5)) < 4 | str2num(datestr(n, 5)) == 12);
count = length(I);
n2 = n(I);
ctd2 = ctd(I,:);
n(I) = [];
ctd(I,:) = [];

%% Store and clean data
noFiles = size(ctd,1); 
zVec = [];
xVec = [];
yVec = [];
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
    end
    
end  % for j

% Quality control
I = find(xVec==-99);
xVec(I) = [];
yVec(I) = [];
zVec(I) = [];

I = find(yVec==-99);
xVec(I) = [];
yVec(I) = [];
zVec(I) = [];

I = find(zVec==-99);
xVec(I) = [];
yVec(I) = [];
zVec(I) = [];

% Remove shallow pts
I = find(zVec<zMin);
I = find(xVec<32);
xVec2 = xVec;
yVec2 = yVec;
zVec2 = zVec;
xVec2(I) = [];
yVec2(I) = [];
zVec2(I) = [];

% Removed shallow pts
I = find(zVec2>zMax);
I = find(xVec2>34.3);
xVec3 = xVec2;
yVec3 = yVec2;
zVec3 = zVec2;
xVec3(I) = [];
yVec3(I) = [];
zVec3(I) = [];

noFiles = size(ctd,1); 

xxVec = [];
yyVec = [];



save plot_NOxRel_info.mat xVec yVec zVec xVec2 yVec2 zVec2 xVec3 ...
    yVec3 zVec3

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

disp(sprintf('%d pts plotted from %d files', length(zVec), length(n)))
disp(sprintf([' -> falls to %d pts when ignoring bottles shallower than %dm'], length(zVec2), zMin))

% NO2 plot
plot(xVec, yVec, '.r')
if strcmp(colOrder(XVAR_COL), 'NOx')==1
    colOrder{XVAR_COL} = '[NO_3] (mmol m^{-3})';
elseif strcmp(colOrder(YVAR_COL), 'NOx')==1
    colOrder{YVAR_COL} = '[NO_3] (mmol m^{-3})';
end
xlabel(colOrder(XVAR_COL))
ylabel(colOrder(YVAR_COL))
hold on
plot(xVec2, yVec2, '.m')
plot(xVec3, yVec3, '.k')
hold off

if strcmp(colOrder(YVAR_COL), 'P')==1
    set(gca, 'ydir', 'reverse')
end

[xVecRaw, I] = sort(xVec3);
yVecRaw = yVec3(I);

p = polyfit(xVecRaw, yVecRaw, 1)

keyboard

%xVecItp = min(xVecRaw):max(xVecRaw);
yVecItp = p(2)+p(1).*xVecRaw;
hold on
plot(xVecRaw, yVecItp, 'color', [1 1 1]*.5)

[r,p] = corrcoef(yVecRaw,yVecItp);
r
keyboard

xlim([min(xVec) max(xVec)])
ylim([min(yVec) max(yVec)])



keyboard
print(gcf, '-dpng', 'NOxRel.png')
set(gcf, 'renderer', 'painters')
print(gcf, '-depsc2', 'NOxRel.eps')

keyboard





% $$$ %% Store and clean ctd2
% $$$ for j = 1:length(n2)       
% $$$     file = load(ctd2(j,:));
% $$$     
% $$$     % Quality control + save into matrix T and S
% $$$     if ~isempty(file)==1 % check if file empty
% $$$         
% $$$         if size(file, 2)~=length(colOrder)
% $$$             disp('Problem with file size')
% $$$             keyboard
% $$$         end
% $$$             
% $$$         xxVar = file(:,XVAR_COL);
% $$$         yyVar = file(:,YVAR_COL);
% $$$ 
% $$$         % store variables
% $$$         xxVec = [xxVec; xxVar];
% $$$         yyVec = [yyVec; yyVar];                
% $$$     end
% $$$     
% $$$ end  % for j
% $$$ 
% $$$ % Quality control
% $$$ I = find(xxVec==-99);
% $$$ xxVec(I) = [];
% $$$ yyVec(I) = [];
% $$$ 
% $$$ I = find(yyVec==-99);
% $$$ xxVec(I) = [];
% $$$ yyVec(I) = [];
% $$$ 
% $$$ I = find(zVec==-99);
% $$$ xxVec(I) = [];
% $$$ yyVec(I) = [];
% $$$ 
% $$$ hold on
% $$$ plot(xxVec, yyVec, '.b')

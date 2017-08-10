function NOx_tanh_fit(datFiles, zMin, zMax, varargin)
    
% NOx_tanh_fit('datFiles.list', 0, 350, 'NOx', 'P')
%
% to be run in /home/cyrf0006/PhD/Nitrates/allStat 
%
%

% Few params that could be passed as function inputs
yearmin = 1996;

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
varOrder{1} = 'NO2';
varOrder{2} = 'NOx';
varOrder{3} = 'PO4';
varOrder{4} = 'SIO4';
varOrder{5} = 'O2';
varOrder{6} = 'S';
varOrder{7} = 'T';
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
xVec2 = xVec;
yVec2 = yVec;
zVec2 = zVec;
xVec2(I) = [];
yVec2(I) = [];
zVec2(I) = [];

% Remove shallow pts
I = find(zVec2>zMax);
xVec3 = xVec2;
yVec3 = yVec2;
zVec3 = zVec2;
xVec3(I) = [];
yVec3(I) = [];
zVec3(I) = [];

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 15])
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
xlabel(colOrder(XVAR_COL))
ylabel(colOrder(YVAR_COL))
hold on
plot(xVec2, yVec2, '.m')
plot(xVec3, yVec3, '.k')
hold off
set(gca, 'ydir', 'reverse')




global S Z


[yVecRaw, I] = sort(yVec3);
xVecRaw = xVec3(I);

S = xVecRaw;
Z = yVecRaw;
%zVec = zMin:zMax;
hold on

a0 = [-5 .01 10 .1 100];
asol = fminsearch ( 'fitFun1', a0)

S0 = a0(1) + a0(2).*Z + a0(3).*tanh((Z-a0(4))./a0(5));
S1 = asol(1) + asol(2).*Z + asol(3).*tanh((Z-asol(4))./asol(5));

%plot(S0, Z, 'm--')
plot(S1, Z, 'r--')

% 2nd fit
a0 = [25 -2 80];
%a0 = [1027.5 -.0924 15.4]; % for sept2011 N080
asol = fminsearch ( 'fitFun2', a0)

S2 = asol(1).*exp(asol(2)./(Z+asol(3)));
plot(S2, Z, 'k--')


% 3rd fit
a0 = [5 10 .1 100];
%a0 = [1027.5 -.0924 15.4]; % for sept2011 N080
asol = fminsearch ( 'fitFun10', a0)

S3 = asol(1) + asol(2).*tanh((Z-asol(3))./asol(4));

plot(S3, Z, 'm--')


keyboard

p = polyfit(xVecRaw, yVecRaw, 1);

%xVecItp = min(xVecRaw):max(xVecRaw);
yVecItp = p(2)+p(1).*xVecRaw;
hold on
plot(xVecRaw, yVecItp, 'k')

[r,p] = corrcoef(yVecRaw,yVecItp)
keyboard





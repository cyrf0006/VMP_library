function NOx_Chla_season(datFiles, zMin, zMax)
    
% Few params that could be passed as function inputs
yearmin = 2000;
%zMin = 25;
%zMax = 300;


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
        
        if size(file, 2) < 6
            continue
        end            
        
        pres = file(:,1);
        xVar = file(:,3);
        yVar = file(:,6);
        
        if isempty(yVar)
            continue
        end
        
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


noFiles = size(ctd,1); 

xxVec = [];
yyVec = [];


% Build timesereis
monthVec = 4:11;
aveSurf = nan(size(monthVec));
stdSurf = nan(size(monthVec));
monthVecSurf = str2num(datestr(timeVecSurf, 5));
aveSurf2 = nan(size(monthVec));
stdSurf2 = nan(size(monthVec));
for i = 1:length(monthVec)
    I = find(monthVecSurf == monthVec(i));
    aveSurf(i) = nanmean(xVecSurf(I));
    stdSurf(i) = nanstd(xVecSurf(I));
    I = find(monthVecSurf == monthVec(i));
    aveSurf2(i) = nanmean(yVecSurf(I));
    stdSurf2(i) = nanstd(yVecSurf(I));
end

keyboard

plot(monthVec, aveSurf)
hold on
plot(monthVec, aveSurf2, 'r')
plot(monthVec, aveSurf./aveSurf2, 'm')

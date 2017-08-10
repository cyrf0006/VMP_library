function wind_climato(windlist)


% ------------------- extract profiles -------------------- %
fid = fopen(windlist);
C = textscan(fid, '%s', 'delimiter', '\n');
windFiles = char(C{1});

noFiles = size(windFiles, 1); %number of eps_files 

% raw matrix to fill
noYear = noFiles/12;
windSpeed = nan(noYear, 365*24);
windDir = nan(noYear, 365*24);
nMat = nan(noYear, 365*24);


count = 1;
year = 1;
mtime = [];
wind17Vec = [];
wind15Vec = [];

for i = 1:noFiles
    
    fname = windFiles(i, :);
    I = find(fname==' ');   
    fname(I) = [];    
    disp(fname)
    data = load(fname);
    
 
    n = datenum(data(:,1), data(:,2), data(:,3), data(:,4), 0, 0);
    wind17 = data(:,17);
    wind15 = data(:,15);
    
    mtime = [mtime n'];
    wind17Vec = [wind17Vec wind17'];
    wind15Vec = [wind15Vec wind15'];
       
end

dir = wind15Vec*10;
vel = wind17Vec/3.6;
save windWGR.mat mtime dir vel
% ---------------------------------------------------------- %




% ------------------- Climato -------------------- %
fid = fopen('windClimato.mat');
if fid == -1 % doesnt exist yet!
    disp('Computing climato....')
    dirClimato = [];
    velClimato = [];
    timeClimato = [];
    for i = 1:12
        disp(sprintf('month %d', i))
        I = find(str2num(datestr(mtime,5)) == i);
        timeVec = mtime(I);
        
        J = find(str2num(datestr(mtime(I),7)) < 7); 
        dirClimato = [dirClimato nanmean(dir(I(J)))];
        velClimato = [velClimato nanmean(vel(I(J)))];    
        timeClimato = [timeClimato datenum(0000, i, 3.5)];
        
        J = find(str2num(datestr(mtime(I),7)) >= 7 & str2num(datestr(mtime(I),7)) < 14); 
        dirClimato = [dirClimato nanmean(dir(I(J)))];
        velClimato = [velClimato nanmean(vel(I(J)))];    
        timeClimato = [timeClimato datenum(0000, i, 10.5)];
        
        J = find(str2num(datestr(mtime(I),7)) >= 14 & str2num(datestr(mtime(I),7)) < 21); 
        dirClimato = [dirClimato nanmean(dir(I(J)))];
        velClimato = [velClimato nanmean(vel(I(J)))];    
        timeClimato = [timeClimato datenum(0000, i, 17.5)];    
        
        J = find(str2num(datestr(mtime(I),7)) >= 21); 
        dirClimato = [dirClimato nanmean(dir(I(J)))];
        velClimato = [velClimato nanmean(vel(I(J)))];    
        timeClimato = [timeClimato datenum(0000, i, 24.5)];    
    end
    save windClimato.mat dirClimato velClimato timeClimato

else
    load windClimato
end
% --------------------------------------------------- %

% ------------------- Climato 2012 -------------------- %
fid = fopen('windClimato2012.mat');
if fid == -1 % doesnt exist yet!
    disp('Computing climato 2012....')
    
    I = find(str2num(datestr(mtime,10)) == 2012);
    mtime = mtime(I);
    
    dirClimato2012 = [];
    velClimato2012 = [];
    timeClimato2012 = [];
    for i = 1:12
        disp(sprintf('month %d', i))
        I = find(str2num(datestr(mtime,5)) == i);
        timeVec = mtime(I);
        
        J = find(str2num(datestr(mtime(I),7)) < 7); 
        dirClimato2012 = [dirClimato2012 nanmean(dir(I(J)))];
        velClimato2012 = [velClimato2012 nanmean(vel(I(J)))];    
        timeClimato2012 = [timeClimato2012 datenum(0000, i, 3.5)];
        
        J = find(str2num(datestr(mtime(I),7)) >= 7 & str2num(datestr(mtime(I),7)) < 14); 
        dirClimato2012 = [dirClimato2012 nanmean(dir(I(J)))];
        velClimato2012 = [velClimato2012 nanmean(vel(I(J)))];    
        timeClimato2012 = [timeClimato2012 datenum(0000, i, 10.5)];
        
        J = find(str2num(datestr(mtime(I),7)) >= 14 & str2num(datestr(mtime(I),7)) < 21); 
        dirClimato2012 = [dirClimato2012 nanmean(dir(I(J)))];
        velClimato2012 = [velClimato2012 nanmean(vel(I(J)))];    
        timeClimato2012 = [timeClimato2012 datenum(0000, i, 17.5)];    
        
        J = find(str2num(datestr(mtime(I),7)) >= 21); 
        dirClimato2012 = [dirClimato2012 nanmean(dir(I(J)))];
        velClimato2012 = [velClimato2012 nanmean(vel(I(J)))];    
        timeClimato2012 = [timeClimato2012 datenum(0000, i, 24.5)];    
    end
    save windClimato2012.mat dirClimato2012 velClimato2012 timeClimato2012

else
    load windClimato2012
end
% --------------------------------------------------- %

% ------------------- plotting -------------------- %
figure(1)
clf
velClimato = runmean(velClimato, 4);
dirClimato = runmean(dirClimato, 4);
I = find(isnan(velClimato2012)==1);
velClimato2012(I) = [];
dirClimato2012(I) = [];
timeClimato2012(I) = [];
velClimato2012 = runmean(velClimato2012, 4);
dirClimato2012 = runmean(dirClimato2012, 4);


xticks = datenum(0000, 1:12, 15);
xlims = [min(timeClimato) max(timeClimato)];

[AX,H1,H2] = plotyy(timeClimato, velClimato, timeClimato, dirClimato);

set(H1, 'lineWidth', 2)
set(H2, 'lineWidth', 2)

set(AX(1), 'xtick', xticks, 'xlim', xlims);
set(AX(2), 'xtick', xticks, 'xlim', xlims, 'xticklabel', []);
datetick('x', 3)
set(AX(1), 'xtick', xticks, 'xlim', xlims);
set(AX(2), 'xtick', xticks, 'xlim', xlims, 'xticklabel', []);

ylabel(AX(1), 'm s^{-1}');
%ylabel(AX(2), 'wind direction');

set(AX(2), 'ytick', [0 45 90 135 180 225 270 315], 'ylim', [170 280])
set(AX(2), 'yticklabel', ['E '; 'NE'; 'N '; 'NW'; 'W '; 'SW'; 'S '; 'SE'])

hold(AX(1), 'on')
plot(AX(1), timeClimato2012, velClimato2012, 'k', 'linewidth', 2)
hold(AX(1), 'off')
hold(AX(2), 'on')
plot(AX(2), timeClimato2012, dirClimato2012, '--k', 'linewidth', 2)
hold(AX(2), 'off')

% --------------------------------------------------- %

keyboard
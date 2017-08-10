function sl_surface_ctd(ctdFiles, gpsFile, coord)

% Will plot TS drawings on a map, like thermosalinograph 'En Route!'
% usage ex: 
%  >> surface_ctd('SBE19plus_01906786_2013_06_10_0077.mat', 'SLEIWEX2013GPS_06-09.TXT', [48.20 48.25 -69.9167 -69.8333])
%  >> surface_ctd('transect_june10_71-73', 'SLEIWEX2013GPS_06-09.TXT', [48.20 48.25 -69.9167 -69.8333])
%  >> surface_ctd('SBE19plus_01906786_2013_06_08_0056.mat', 'SLEIWEX2013GPS_06-09.TXT', [48.195 48.25 -69.9167 -69.8333])
% F. Cyr - june 2013

% Single profile or list?
if strcmp(ctdFiles(end-3:end),'.mat') == 1
    command = sprintf('ls %s > /tmp/tmp', ctdFiles);
    system(command);
    ctdFiles='/tmp/tmp';
end


% Open Profiles list    
fid = fopen(ctdFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});
noProfiles = size(files,1);

% deal with GPs files
fid = fopen(gpsFile);
C = textscan(fid, '%s', 'delimiter', '\n');
gps = char(C{1});
gps(1:8,:) = []; % weak way to do it!

% find comma position
I=findstr(gps(1,:),',');

%gps = str2num(gps);

lat = str2num(gps(:,1:I(1)-1));
lon = str2num(gps(:,I(1)+1:I(2)-1));
dategps = gps(:,I(2)+1:I(3)-1);
hourgps = gps(:,I(3)+1:end);

mtimeGPS = datenum([dategps hourgps], 'yyyymmddHHMMSS');

% get TS
mtimeCTD = [];
TCTD = [];
SCTD = [];
PCTD = [];
for i = 1:noProfiles
    
    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];

    disp(fname);
    cast = load(fname);
    mtimeCTD = [mtimeCTD cast.mtime];
    TCTD = [TCTD cast.data(:,3)'];
    SCTD = [SCTD cast.data(:,6)'];
    PCTD = [PCTD cast.data(:,2)'];

end

I = find(PCTD<.2);
PCTD(I) = NaN;
TCTD(I) = NaN;
SCTD(I) = NaN;


% time average
dt = 10; %sec.
timeVec = min(mtimeCTD):dt/86400:max(mtimeCTD);
TVec = [];
SVec = [];
for i = 1:length(timeVec)
    I = find(mtimeCTD>=timeVec(i)-dt/86400 & mtimeCTD<timeVec(i)+dt/86400);
    TVec(i) = nanmean(TCTD(I));    
    SVec(i) = nanmean(SCTD(I));
end

% interpol. GPS position
latVecGPS = interp1(mtimeGPS, lat, timeVec);
lonVecGPS = interp1(mtimeGPS, lon, timeVec);


% plot Figure
figure(1)
clf
paperwidth = 16;%cm
paperheight = 14;%cm
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator','long',coord(3:4),'lat',coord(1:2));
%m_gshhs_h('patch',[.5 .5 .5]); %coastlines
%m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy')
m_gshhs_h('patch',[.5 .5 .5]); %coastlines

moorLat = 48+12.642/60;
moorLon = -69-53.719/60;

m_line(moorLon,moorLat, 'marker','p','MarkerFaceColor',[0 1 0]*.5, 'markersize',6,'color','k');
hold on

% define color
Tmax = max(TVec);
Tmin = min(TVec);
Smax = max(SVec);
Smin = min(SVec);

[x, y]=m_ll2xy(lonVecGPS, latVecGPS);
rectbar = [0.92 0.15 0.01 0.75]; 
keyboard
cb=colour_profile(x,y,TVec,Tmin,Tmax,1, rectbar, x.*0); 
hold off


keyboard
% Simple script that uses a list of SBE casts in .mat to plot over
% Georeferenced pictures. See inputs below.
% Could be easily modified for a function.
% F. Cyr june 2013

figure(1)
clf
geoFname = '~/PhD/SLEIWEX2013/Pour_Fred/g_rect.mat';
imgFname = '~/PhD/SLEIWEX2013/Pour_Fred/IMG_0441.JPG';
%imgFname = '~/PhD/SLEIWEX2013/Pour_Fred/IMG_0531.JPG';

g_viz_field(imgFname,geoFname);


files = ['~/PhD/SLEIWEX2013/SBE19plus_psgHeader/SBE19plus_01906786_2013_06_10_0071.mat';...
'~/PhD/SLEIWEX2013/SBE19plus_psgHeader/SBE19plus_01906786_2013_06_10_0072.mat';...
'~/PhD/SLEIWEX2013/SBE19plus_psgHeader/SBE19plus_01906786_2013_06_10_0073.mat';...
];

% $$$ files = ['../SBE19plus_psgHeader/SBE19plus_01906786_2013_06_10_0077.mat'];

gpsFile = ['~/PhD/SLEIWEX2013/SBE19plus_psgHeader/SLEIWEX2013GPS_07-09.TXT'];
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


% $$$ % plot Figure
% $$$ figure(1)
% $$$ clf
% $$$ paperwidth = 16;%cm
% $$$ paperheight = 14;%cm
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
% $$$ m_proj('mercator','long',coord(3:4),'lat',coord(1:2));
% $$$ %m_gshhs_h('patch',[.5 .5 .5]); %coastlines
% $$$ %m_coast('patch',[.7 .7 .7],'edgecolor','none');
% $$$ m_grid('box','fancy')
% $$$ m_gshhs_h('patch',[.5 .5 .5]); %coastlines

moorLat = 48+12.642/60;
moorLon = -69-53.719/60;


m_line(moorLon,moorLat, 'marker','p','MarkerFaceColor',[1 1 1], 'markersize',12,'color','k');


% define color
Tmax = max(TVec);
Tmin = min(TVec);
Smax = max(SVec);
Smin = min(SVec);

[x, y]=m_ll2xy(lonVecGPS, latVecGPS);
rectbar = [0.92 0.15 0.01 0.75]; 

freezeColors
cb=colour_profile(x,y,TVec,Tmin,11,1, rectbar, x.*0); 
hold off


function sl_timeseries4video(adcpFile, gpsFile)
 
% Generates a .mat in which time, Temp. and vert_vel at mooring
% latlon position are given. 
% For Fred, run in /home/cyrf0006/PhD/SLEIWEX2013/dataprocessing/Mooring
% usage ex:
% sl_timeseries4video('./ADCP_8336/ADCP_8336a.mat','gps8336a.list'); 
%
% F. Cyr - june 2013    

% Few infos
zrange = 1:10;

    
    
% ADCP File
disp('process ADCP file...')
load(adcpFile);
Wadcp = ADCP.vert_vel;
Tadcp = ADCP.temperature;
mtimeADCP = ADCP.mtime;
clear ADCP
Wadcp = nanmean(Wadcp(zrange,:));



    
% GPS Files
disp('process GPS files...')
fid = fopen(gpsFile);
C = textscan(fid, '%s', 'delimiter', '\n');
gps = char(C{1});

noProfiles = size(gps,1);
lat = [];
lon = [];
mtimeGPS = [];
for i = 1:noProfiles
    
    fname= gps(i, :);
    I = find(fname==' ');   
    fname(I) = [];
   
    % deal with GPs files
    fid = fopen(fname);
    C = textscan(fid, '%s', 'delimiter', '\n');
    GPS = char(C{1});
    
    % remove GPs header
    notRemoved = 1;
    while notRemoved
        I = find(strcmp(GPS(1,1), '#'));
        if I == 1
            GPS(1,:) = [];
        else
            notRemoved = 0;
        end
    end
    
    % find comma position
    I=findstr(GPS(1,:),',');

    %gps = str2num(gps);

    lat = [lat; str2num(GPS(:,1:I(1)-1))];
    lon = [lon; str2num(GPS(:,I(1)+1:I(2)-1))];
    dategps = GPS(:,I(2)+1:I(3)-1);
    hourgps = GPS(:,I(3)+1:end);

    mtimeGPS = [mtimeGPS; datenum([dategps hourgps], 'yyyymmddHHMMSS')];
    
end



% PictureTime
disp('Referencing to camera snapshot...')
dt = 1/1440;
timeVec = round(mtimeADCP(1)/dt)*dt:dt:round(mtimeADCP(end)/dt)*dt;
latVec = nan(size(timeVec));
lonVec = nan(size(timeVec));
wVec = nan(size(timeVec));
TVec = nan(size(timeVec));

for i = 1:length(timeVec)
    
    % mooring position
    [Y, I] = min(abs(timeVec(i)-mtimeGPS));
    if Y <= dt
        lonVec(i) = lon(I);
        latVec(i) = lat(I);
    end

    % mooring position
    I = find(mtimeADCP>=timeVec(i)-dt/2 & mtimeADCP<timeVec(i)+dt/2);
    if ~isempty(I) == 1
        TVec(i) = nanmean(Tadcp(I));
        wVec(i) = nanmean(Wadcp(I));
    end
end


save mooring_framereferenced_LatLonTw.mat timeVec lonVec latVec TVec ...
    wVec


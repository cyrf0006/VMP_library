function [u, v, u_gps, v_gps] = adcp_velcor(east_vel, north_vel, mtime_adcp, gps_file, varargin)

% function [u v] = adcp_velcor(east_vel, north_vel, mtime_adcp, gps_file)
% 
% Designed to correct ADCP velocities for boat
% displacement/drift. Inputs are easy to understand. Outputs u and
% v are corrected ADCP velolcities. Inputs velocities can be
% whether matrix or vector, but one of the dimension must
% correspond to mtime_adcp
%
% usage ex: >> [u, v] = adcp_velcor(east_vel, north_vel, mtime_adcp, '../GPS/20110714')
%    
%
% F. Cyr, july 2011
%
% Modified: - june 2013: Take into account Fugawi .TXT files
% directly. USAGE FOR MELANY: 
% >> [u, v] = adcp_velcor(east_vel, north_vel, mtime_adcp, 'SLEIWEX2013GPS_06-11.TXT', 'fugawi')

% ----------------------------------------------------------- %   
% $$$ 
% $$$ disp('velcor')


% check if input velocities are needed to be flipped
if size(east_vel,2)~=length(mtime_adcp)
    
    if size(east_vel,2)==length(mtime_adcp)
        east_vel = east_vel';
        north_vel = north_vel';
    else
        disp('time not corresponding to velocities')
        return
    end    
end


% Load GPS infos
if strcmp(varargin, 'fugawi') == 1 % Take directly .TXT files, no
                                   % shell script manipulation...
    gps = fugawi_gps(gps_file);
    
    lat_gps = gps.lat;
    lon_gps = gps.lon;
    mtime_gps = gps.mtime;
 
else % Fred's original version
    track = load(gps_file);
    mtime_gps = datenum(track(:,7), track(:,6), track(:,5), track(:,8), track(:,9), track(:,10));
    lat_gps = track(:,1)+track(:,2)/60;
    lon_gps = track(:,3)+track(:,4)/60;
end

% ---- compute boat displacement (U, V) ---- %
dt_gps = mtime_gps(2) - mtime_gps(1); % in mtime units (.i.e. /86400sec)
mtime_gps = mtime_gps(1:end-1)+dt_gps; % shift mtime to center pts

dlat = diff(lat_gps);
dlon = diff(lon_gps);

range = m_lldist(lon_gps,lat_gps)*1000; % in m 
vel_norm = range./dt_gps/86400; % in m/s

% IMPORTANT NOTE:
% The angle should be calculated with the four-quadrant inverse tangent!
% Otherwise, the angle is caught between -90 and 90 degrees.
% Must use here the function atan2d and not atand
%vel_angle = atand(dlat./dlon); % NaNs if dlon==0
vel_angle = atan2d(dlat,dlon); % NaNs if dlon==0

% IMPORTANT NOTE: This condition is not required here when
% using the atan2d.
% remove NaNs 
%I = find(isnan(vel_angle)==1);
%vel_angle(I)=90; % tan not defined if angle = 90deg

u_gps = vel_norm.*cosd(vel_angle);
v_gps = vel_norm.*sind(vel_angle);

% ------------------------------------------- %

% interpolation of GPS velocity to adcp_mtime
%u_gps = interp1(mtime_gps, u_gps, mtime_adcp);
%v_gps = interp1(mtime_gps, v_gps, mtime_adcp);

% NOTE: Comparing the length of mtime_adcp and mtime_gps is not a robust
% way to determine whether one time vector has a higher sampling rate,
% unless both vectors start and end at the same time which may not be
% generally the case. A more robust way is to compare the sampling period 
% (dt_adcp compared with dt_gps).

dt_adcp = mtime_adcp(2)-mtime_adcp(1);
if dt_adcp < dt_gps
%if length(mtime_adcp)>length(mtime_gps)
    u_gps = interp1(mtime_gps, u_gps, mtime_adcp);
    v_gps = interp1(mtime_gps, v_gps, mtime_adcp);
else
    u_gps_dec = nan(length(mtime_adcp),1);
    v_gps_dec = nan(length(mtime_adcp),1);
    dt_adcp = mtime_adcp(2)-mtime_adcp(1);
    for i = 1:length(mtime_adcp)
        I = find(mtime_gps>=mtime_adcp(i)-dt_adcp/2 & mtime_gps< ...
                 mtime_adcp(i)+dt_adcp/2);
        u_gps_dec(i) = nanmean(u_gps(I))
        v_gps_dec(i) = nanmean(v_gps(I));
    end
    u_gps = u_gps_dec;
    v_gps = v_gps_dec;    
end

% correction of velocities
u = nan(size(east_vel));
v = nan(size(north_vel));

% flip u_gps and v_gps if needed
if size(east_vel,2) ~= size(u_gps,2)
    u_gps = u_gps';
    v_gps = v_gps';
end

for i = 1:size(east_vel, 1)
    u(i,:) = east_vel(i,:)+u_gps;
    v(i,:) = north_vel(i,:)+v_gps;
end
    
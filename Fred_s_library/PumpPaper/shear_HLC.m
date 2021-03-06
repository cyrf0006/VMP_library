clear

load('~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/ADCP/adcp_2009-10-01.mat');
gpsFile = '~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/GPS/20091001';


east_vel = ADCP.east_vel;
north_vel = ADCP.north_vel;
mtime = ADCP.mtime;
z = ADCP.config.ranges;

% $$$ I = find(mtime > datenum(2009, 10, 25, 14, 45, 0));
% $$$ east_vel = east_vel(:,I);
% $$$ north_vel = north_vel(:,I);
% $$$ mtime = mtime(I);

I = find(~isnan(nansum(east_vel,1))==1);
east_vel = east_vel(:,I);
north_vel = north_vel(:,I);
mtime = mtime(I);

[u v] = adcp_velcor(east_vel, north_vel, mtime, gpsFile);


% Filtering
filterTime = 1; % in minutes 
dt = (mtime(2)-mtime(1));
% Low pass 
fs = 1/(dt*86400);
freq_low = 1/(filterTime*60); %Hz
Wn_low = freq_low/(fs/2);
[b,a] = butter(4, Wn_low);


u_filt = nan(size(u));
v_filt = nan(size(v));

for i = 1:size(u, 1)
    % remove NANs
    uvec = u(i,:);
    I = find(isnan(uvec)==1);
    xitp = 1:length(uvec);
    x = xitp;
    uvec(I) = [];
    x(I) = [];
    uvec = interp1(x, uvec, xitp);    
    u_filt(i, :) = filtfilt(b, a, uvec);
    
    vvec = v(i,:);
    I = find(isnan(vvec)==1);
    xitp = 1:length(vvec);
    x = xitp;
    vvec(I) = [];
    x(I) = [];
    vvec = interp1(x, vvec, xitp);    
    v_filt(i, :) = filtfilt(b, a, vvec);
end


du = diff(u_filt,1,1);
dv = diff(v_filt,1,1);
dz = z(2) - z(1);

S2Mat = (du/dz).^2 + (dv/dz).^2;
zS2 = z(1:end-1)+dz/2;


save ADCP_shear.mat S2Mat zS2 mtime


I = find(mtime > datenum(2012, 10, 25, 14, 48, 0) & mtime < datenum(2012, 10, 25, 16, 0, 0) );
J = find(zS2<80);

% Must now find the bottom
Hobs = [79 79 65]+1; % +1m ADCP depth
Xobs = [datenum(2012, 10, 25, 14, 48, 0) datenum(2012, 10, 25, 14, 50, 0) datenum(2012, 10, 25, 16, 0, 0)];

H = interp1(Xobs, Hobs, mtime(I));
S2 = S2Mat(:,I);
%S2 = S2Mat;

for i = 1:length(I)
    II = find(zS2>H(i));
    S2(II, i) = NaN;
end

    


%imagesc(mtime(I), zS2(J), log10(S2))
%contourf(mtime(I), zS2(J), log10(S2(J,I)))
contourf(mtime(I), zS2(J), log10(S2(J,:)), 'linestyle', 'none')
%imagesc(mtime(I), zS2(J), log10(S2(J,:)))

caxis([-5 -2])
set(gca, 'ydir', 'reverse')
datetick('x', 15) 
xlim([datenum(2012, 10, 25, 14, 48, 0) datenum(2012, 10, 25, 16,0, 0)])
colorbar

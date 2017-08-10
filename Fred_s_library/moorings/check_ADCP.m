function check_ADCP(adcp_matfile, gpsfile)
% ex: 

% ---- parameters ---- %
theta = 33.5; %deg.

load(adcp_matfile);


east_vel  = ADCP.east_vel;
north_vel  = ADCP.north_vel;
mtime =  ADCP.mtime;
z = ADCP.config.ranges;
w = ADCP.vert_vel;

% correction for ADCP year (0011 instead 2011)
n = datestr(mtime,31);
mtime = datenum(2011, str2num(n(:,6:7)),str2num(n(:,9:10)),str2num(n(:,12:13)),str2num(n(:,15:16)),str2num(n(18:19)));

% correct for magnetic deviation mistake
[u, v] = rotate_vecd(east_vel, north_vel, 0.8);

% correct for boat drift
[ucor, vcor] = adcp_velcor(u, v, mtime, gpsfile);



% ---- Clean velovities (remove NaNs) --- %
I = find(sum(isnan(ucor),1)<size(ucor,1)/2); % remove entire column
ucor = ucor(:,I);
vcor = vcor(:,I);
wcor = w(:,I);
mtime = mtime(I);

% patch for missing values
for i = 1:size(ucor,1)
    I = find(~isnan(ucor(i, :))==1);
    if length(I) < size(ucor, 2)
        y1 = ucor(i,I);
        y2 = vcor(i,I);
        y3 = wcor(i,I);
        x = 1:length(I);
        xi = 1:length(mtime);
        ucor(i,:) = interp1(x, y1, xi);
        vcor(i,:) = interp1(x, y2, xi);
        wcor(i,:) = interp1(x, y3, xi); 
    end
end
% ---------------------------------------- %



% rotate along-cross shore
[valong vcross] = rotate_vecd(ucor, vcor, theta);

% ---- Apply a smoothing filter ---- %
a = 1;
window = 25; % no. of pts
b = ones(1,window);
b = b/length(b);
valong = filtfiltm(b,a,valong, 2);
vcross = filtfiltm(b,a,vcross, 2);
wfilt = filtfiltm(b,a,wcor, 2);
% ------------------------------------ %



% ---- INTENSITY ---- %
A = ADCP.intens;

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 20])
%subplot(3,1,1)
B = reshape(A(:,2,:), size(A,1), size(A,3));
imagesc(mtime, z, B) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
%ylim([5 100])
%caxis([60 180])
%caxis([60 70])
colorbar
title('backscatter3')


keyboard

subplot(3,1,2)
contourf(mtime, z, wfilt) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
colorbar
caxis([-0.2 0.2])
title('vertical velocities (m/s)')

keyboard

% ----------------------- %

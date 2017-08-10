function adcp_quickdisplay(adcpfile, gpsfile)

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
% ---- parameters ---- %
theta = 33.5; %deg.
outfile = 'along_vel3.avi';


! rm /tmp/*.png

load(adcpfile)


east_vel  = ADCP.east_vel;
north_vel  = ADCP.north_vel;
mtime =  ADCP.mtime;
z = ADCP.config.ranges;
w = ADCP.vert_vel;

% Backscatter intensity
I1 = ADCP.intens(:,1,:);
I2 = ADCP.intens(:,2,:);
I3 = ADCP.intens(:,3,:);
I4 = ADCP.intens(:,4,:);
 
I1 = squeeze(I1);
I2 = squeeze(I2);
I3 = squeeze(I3);
I4 = squeeze(I4);
Imean = (I1 + I2 + I3 + I4)/4;

% correction for ADCP year (0011 instead 2011)
n = datestr(mtime,31);
mtime = datenum(2011, str2num(n(:,6:7)),str2num(n(:,9:10)),str2num(n(:,12:13)),str2num(n(:,15:16)),str2num(n(18:19)));

% correct for magnetic deviation mistake
[u, v] = rotate_vecd(east_vel, north_vel, 0.8);

% correct for boat drift
[ucor, vcor] = adcp_velcor(u, v, mtime, gpsfile);



% ---- Clean velovities (remove NaNs) --- %
raw_mtime = mtime;
raw_Imean = Imean;
raw_ucor = ucor;
raw_vcor = vcor;
raw_w = w;

I = find(sum(isnan(ucor),1)<size(ucor,1)/2); % remove entire column
ucor = ucor(:,I);
vcor = vcor(:,I);
wcor = w(:,I);
Imean = Imean(:,I);
mtime = mtime(I);

% patch for missing values
for i = 1:size(ucor,1) % horizontal
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

for i = 1:size(ucor,2) % vertical
    I = find(~isnan(ucor(:, i))==1);
    if length(I) < size(ucor, 1)
        y1 = ucor(I,i);
        y2 = vcor(I,i);
        y3 = wcor(I,i);
        x = 1:length(I);
        xi = 1:length(z);
        ucor(:,i) = interp1(x, y1, xi);
        vcor(:,i) = interp1(x, y2, xi);
        wcor(:,i) = interp1(x, y3, xi); 
    end
end
% ---------------------------------------- %


% ------- rotate along-cross shore ------- %
[valong vcross] = rotate_vecd(ucor, vcor, theta);
[raw_valong raw_vcross] = rotate_vecd(raw_ucor, raw_vcor, theta);
% ---------------------------------------- %


% ------- Apply a smoothing filter ------- %
a = 1;
window = 50; % no. of pts
b = ones(1,window);
b = b/length(b);
valong_filt = filtfiltm(b,a,valong, 2);
vcross_filt = filtfiltm(b,a,vcross, 2);
wfilt = filtfiltm(b,a,wcor, 2);
% ---------------------------------------- %



figure(1) % Interpolated fields
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 20])

subplot(3,1,1)
contourf(mtime, z, valong_filt) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
colorbar
%caxis([-0.2 0.2])
title('along-shore velocities (m/s)')
adjust_space

subplot(3,1,2)
contourf(mtime, z, vcross_filt) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
colorbar
%caxis([-0.2 0.2])
title('cross-shore velocities (m/s)')
adjust_space

subplot(3,1,3)
contourf(mtime, z, wfilt) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
colorbar
%caxis([-0.2 0.2])
title('vertical velocities (m/s)')
adjust_space
% ----------------------- %



figure(2) % raw fields
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 20])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

subplot(3,1,1)
imagesc(raw_mtime, z, raw_Imean) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
caxis([60 180])
%caxis([60 70])
colorbar
title('mean backscatter')
adjust_space

subplot(3,1,2)
imagesc(raw_mtime, z, raw_valong) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
colorbar
%caxis([-0.2 0.2])
title('along-shore velocities (m/s)')
adjust_space

subplot(3,1,3)
imagesc(raw_mtime, z, raw_w) 
set(gca, 'ydir', 'reverse')
datetick('x', 15)
%xlim([datenum(2011, 07, 14, 14, 30, 0) datenum(2011, 07, 14, 19, 30, 0)] )
ylim([5 100])
colorbar
%caxis([-0.2 0.2])
title('vertical velocities (m/s)')
adjust_space
% ----------------------- %


keyboard

count = 1;
for i = 1:size(valong,2)
    
    % if nansum(valong(:,i))>1
        
        figure(1)
        clf
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 6 8])
        
        plot(valong(:,i), z, 'linewidth', 2);
        hold on
        plot([0 0],[min(z) max(z)], '--k')
        set(gca, 'ydir', 'reverse')
        axis([-.5 .5 0 110]);
        title(sprintf(datestr(mtime(i))))
        xlabel('along shore velocity (m/s)')
        ylabel('depth (m)')
        
        if count>=1000
            fig = sprintf('/tmp/fig%d', count);
        elseif count>=100
            fig = sprintf('/tmp/fig0%d', count);
        elseif count>=10
            fig = sprintf('/tmp/fig00%d', count);
        else 
            fig = sprintf('/tmp/fig000%d', count);
        end
        
        print('-dpng', '-r300', fig)
        count = count + 1;
        %end
    

end

!mencoder mf:///tmp/*.png -mf w=1000:h=1000:fps=25:type=png -ovc copy -oac copy -o output.avi
eval(['! mv output.avi ' outfile])

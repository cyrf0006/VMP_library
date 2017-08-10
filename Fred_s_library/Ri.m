function Ri(adcpfile, gpsfile, profiles)
    
% function Ri(adcpfile, gpsfile, profiles)    
% Compute a time serie of Gradient Richardson number. Inputs are :
%   - adcpfile: a .mat file containing adcp_strct (obtained through rdradcp.m)     
%   - gpsfile: standard gps.txt file (as in my other scripts)     
%   - profiles: a list of profiles obtained with
%      -> ls -1 profile*.mat | sed 's/\.mat//' > profiles
%  or  -> ls -1 ../VMP/profile*.mat | sed 's/\.mat//' > profiles
%                           (if in ./ADCP)
%
%
%   usage ex: Ri('D0714000.mat', '../GPS/20110714', 'profiles')
%
% F. Cyr (August 2011)
% ----------------------------------------------- %
    
% ******************** Few constants ********************** % 
g = 9.81; %m2/s
rho_0 = 1025; %kg/m3
% ********************************************************** % 
    
    
    
% ********************  Work on ADCP ********************** % 
disp('Compute velocity from ADCP...')
load(adcpfile)

east_vel  = ADCP.east_vel;
north_vel  = ADCP.north_vel;
mtime =  ADCP.mtime;
z = ADCP.config.ranges;
zbin = z(2)-z(1);
w = ADCP.vert_vel;
intens = reshape(ADCP.intens(:,3,:), size(ADCP.intens,1), size(ADCP.intens,3));

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
intens = intens(:,I);
mtime_adcp = mtime(I);

% patch for missing values
for i = 1:size(ucor,1)
    I = find(~isnan(ucor(i, :))==1);
    if length(I) < size(ucor, 2)
        y1 = ucor(i,I);
        y2 = vcor(i,I);
        y3 = wcor(i,I);
        x = 1:length(I);
        xi = 1:length(mtime_adcp);
        ucor(i,:) = interp1(x, y1, xi);
        vcor(i,:) = interp1(x, y2, xi);
        wcor(i,:) = interp1(x, y3, xi); 
    end
end
% ---------------------------------------- %

% ********************************************************** % 





% ********************  Work on CTD ********************** % 
disp('Compute density from CTD...')
fid = fopen(profiles);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

% empty matrix
rho_mat = nan(size(ucor,1), size(files, 1));
timectd = nan(size(ucor,1), size(files, 1));

for i = 1:size(files, 1)
    
    
    fname = files(i, :); 
    load(fname);
   
    %    timectd(i) = mtime(1);
    
    DENS = sw_dens(SBS, SBT, P);
    
    % bin and save density
    for j = 1:size(z)
        I = find(P>=z(j)-zbin/2 & P<=z(j)+zbin/2);
        rho_mat(j,i) = nanmean(DENS(I));
        timectd(j,i) = mean(mtime(I));
    end
    
    
end
% ********************************************************** % 


% ********************  Compute Richardson number ********************** % 
disp('Compute Ri...')
% shear prod.
dudz = diff(ucor,1)./zbin;
dvdz = diff(vcor,1)./zbin;
S2 = dudz.^2 + dvdz.^2;

% buoy. freq.
N2 = g./rho_0./zbin.*diff(rho_mat,1);

% Ri
Ri = nan(size(N2));
for i = 1:size(N2,1) % loop on z
    for j = 1:size(N2,2) %loop on time
        [Y, I] = min(abs(timectd(i,j)-mtime_adcp));
        Ri(i,j) = N2(i,j)./S2(i,I);
    end
end

%    Ri = N2./S2;

% pressure associated to Ri
z_ri = z(1:end-1)+zbin/2;
% ********************************************************************** % 

% ******************** Plotting ********************** % 
disp('Plotting...')
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

% --- ADCP baskscatter --- %
imagesc(mtime_adcp, z, intens)
set(gca, 'ydir', 'reverse')
ylim([0 100])
datetick('x',15)
xlim([datenum(2011,7,14,9,45,0), datenum(2011,7,14,23,0,0)]) 
xlabel('July 14 2011 (UTC)')
ylabel('depth(m)')

% --- add Ri --- %
colormap(gray)
freezeColors
epsMin = 0;
epsMax = 1;

hold on
for j = 1:size(Ri,2) 
    
    %data = eps_mat(:,j);
    data = Ri(:,j);
    x(1:length(data)) = timectd(1:end-1, j);
    y = z_ri;
      
    rectbar = [0.92 0.15 0.01 0.75]; 

    if j == size(Ri,2) % with cbar
        cb = colour_profile(x',y',data,epsMin,epsMax,1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',data,epsMin,epsMax,0,rectbar, x'.*0);
    end
 
end
hold off
print('-dpng', '-r300','isopyc_Ri.png')
% **************************************************** % 


% ******************** Plotting 2 ********************** % 
disp('Plotting2...')
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

% --- ADCP baskscatter --- %
V = 1018:0.2:1028;
contourf(timectd(1,:), z_ri, rho_mat(1:end-1, :), V, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
ylim([0 100])
datetick('x',15)
xlim([datenum(2011,7,14,9,45,0), datenum(2011,7,14,23,0,0)]) 
xlabel('July 14 2011 (UTC)')
ylabel('depth(m)')

% --- add Ri --- %
freezeColors
hold on
contour(timectd(1,:), z_ri, Ri, [.25 .25], 'edgecolor', 'k', 'LineStyle', '-','linewidth', 2 )
hold off
print('-dpng', '-r300','isopyc_Ri0p25.png')

% **************************************************** %

% ******************** Plotting 2 ********************** % 
disp('Plotting3 (shear)...')
figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

% --- ADCP S2 --- %
V = 1018:0.2:1028;
imagesc(timectd(1,:), z_ri, S2)
set(gca, 'ydir', 'reverse')
ylim([0 100])
datetick('x',15)
xlim([datenum(2011,7,14,9,45,0), datenum(2011,7,14,23,0,0)]) 
xlabel('July 14 2011 (UTC)')
ylabel('depth(m)')
caxis([0 0.005])

print('-dpng', '-r300','S2.png')

% **************************************************** %
keyboard
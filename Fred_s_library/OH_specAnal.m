%Infos on moorings in /home/cyrf0006/PhD/Hydroxide/courantometresSGDO
% $$$ mooringLat = [48.050000 48.050000 47.747200 47.747200 47.747800 47.743500,];
% $$$ mooringLon = [-61.538300 -61.538300  -59.596200 -59.596200 -59.599300 -59.622000];
% $$$ depth = [61 15 300 300 250 250]



% mooring mcm_bio70036_12_91_1200.odf
dt = 1200; %s
t0 = datenum(1970, 10, 23, 11, 31, 0);
tf = datenum(1970, 11, 15, 23, 11, 0);
data = load('mooring91.dat');

u = data(:,1);
angle = data(:,2);
freq = 1/dt;
timeVec = t0:dt/86400:tf;

nx = max(size(u)); % Hanning
na = 10;
w = hanning(floor(nx/na));

u = detrend(u);

[ps, f] = pwelch(u, w, 0, [], freq*86400);

figure(1)
clf
loglog(f, ps);
xlabel('cpd')
ylabel('psd') 
hold on
plot_harmonics



function Ea_tide(tidefile, Ea_file)
  
% function IW_tide(tidefile, IW_file)
%
% ex: Ea_tide('../tide_shear/tide_2009-2012.dat', 'Ea_raw.mat')
%
%

springtide = 3.5; %m
    
tide  = load(tidefile);
    
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);


% find high tide time 
count = 1;
for i = 2:length(mtime)-1

    if level(i)>level(i-1) & level(i)>level(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        count = count+1;
    end
end
% T contains the hour of each high tide



load(Ea_file)

dt = 1/48;
t1 = round(Ea_raw(2,1)*24)/24;
t2 = round(Ea_raw(2,end)*24)/24;
t_bin = t1:dt:t2;


% bin the Energy to time vector
for i = 1:length(t_bin);
    I = find(Ea_raw(2,:)>=t_bin(i)-dt/2 & Ea_raw(2,:)<t_bin(i)+dt/2);
    Ea_bin(i) = mean(Ea_raw(1, I)); 
end


% For each Ea, find position relative to cloasest high tide
count2=1;
count3=1;
for i = 1:length(t_bin)
    [Y, I] = min(abs(T-t_bin(i)));
    A(i) = (t_bin(i)-T(I))*24;
    B(i) = L(I); %level of the closest hightide
    if L(I) > springtide 
        spring(count2) = i;
        count2 = count2+1;
    else
        neap(count3) = i;
        count3 = count3+1;
    end
end


Ea_smooth = loess(A,Ea_bin,sort(A),0.3,1);
Ea_smooth_spring = loess(A(spring),Ea_bin(spring),sort(A(spring)),0.3,1);
Ea_smooth_neap = loess(A(neap),Ea_bin(neap),sort(A(neap)),0.3,1);


% $$$ 
% $$$ figure(1)
% $$$ plot(sort(A), Ea_smooth)
% $$$ title('Sept. 12 - Oct. 19')
% $$$ xlabel('time to hightide')
% $$$ ylabel('smooth Ea')
% $$$ print('-dpng', '-r300','Ea_M2_all.png')
% $$$ 
% $$$ 
% $$$ figure(2)
% $$$ plot(sort(A(spring)), Ea_smooth_spring)
% $$$ title('Spring tide')
% $$$ xlabel('time to hightide')
% $$$ ylabel('smooth Ea')
% $$$ print('-dpng', '-r300','Ea_M2_spring.png')
% $$$ 
% $$$ figure(3)
% $$$ plot(sort(A(neap)), Ea_smooth_neap)
% $$$ title('Neap tide')
% $$$ xlabel('time to hightide')
% $$$ ylabel('smooth Ea')
% $$$ print('-dpng', '-r300','Ea_M2_neap.png')
% $$$ 
% $$$ figure(4)
% $$$ plot(sort(A(neap)), Ea_smooth_neap, '--k')
% $$$ hold on
% $$$ plot(sort(A(spring)), Ea_smooth_spring, 'k')
% $$$ plot(sort(A), Ea_smooth, 'k', 'linewidth', 2)
% $$$ hold off
% $$$ xlabel('time to hightide')
% $$$ ylabel('smooth Ea')
% $$$ print('-dpng', '-r300','Ea_M2_3curves.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'Ea_M2_3curves.eps')


% save in .mat
time_Ea = sort(A);
time_Ea_spring = sort(A(spring));
time_Ea_neap = sort(A(neap));

save Ea_M2.mat time_Ea time_Ea_spring time_Ea_neap Ea_smooth Ea_smooth_neap ...
    Ea_smooth_spring




% ---------------- Full Ea cycle... ------------------- %
% Low pass velocity (13 hours)
% band pass filter (tides and noise)
dt = 3; %sec
fs = 1/dt;
freq_low = 1/(60*60*13); %Hz
Wn_low = freq_low/(fs/2);
[b,a] = butter(4, Wn_low);
Ea_filt = filtfilt(b, a, Ea_raw(1,:));
Ea_time = Ea_raw(2,:);

I = find(Ea_time<datenum(2011, 09, 21));
Ea_time(I) = [];
Ea_filt(I) = [];



I = find(mtime>=Ea_time(1) & mtime<=Ea_time(end));

figure(5)
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.05 ; % horiz. space between subplots
dy = 0.06; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 15 12])
subplot(211)
plot(mtime(I), level(I), 'k', 'linewidth', 2)
xlim([mtime(I(1)) mtime(I(end))])
ylabel('\eta (m)')
ylim([0 5])
set(gca, 'xticklabel', '')
set(gca, 'xtick', datenum(2011, 09, 22):datenum(2011, 10, 12));
adjust_space

subplot(212)
plot(Ea_time, Ea_filt, 'k', 'linewidth', 2)
hold on
plot([Ea_time(1) Ea_time(end)], [nanmean(Ea_filt) nanmean(Ea_filt)], '--k', 'linewidth', 2)
datetick('x',7)
xlim([mtime(I(1)) mtime(I(end))])
xlabel('Sept./Oct. 2011')
ylabel('E_a (J m^{-2})')
%set(gca, 'xtick', xticks)
adjust_space

print('-depsc2', 'Eafull_tide.eps')

% ---------------------------------------------------------------- %


% ---------------- Zoom and M2 cycle ------------------- %
Izoom = find(mtime>=datenum(2011, 10, 01) & mtime<datenum(2011, 10, 2));
IzoomEa = find(Ea_time>=datenum(2011, 10, 01) & Ea_time<datenum(2011, 10, 2));

figure(6)
clf
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.08 ; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 15 12])
subplot(311)
plot(mtime(Izoom), level(Izoom), 'k', 'linewidth', 2)
xlim([mtime(Izoom(1)) mtime(Izoom(end))])
ylabel('\eta (m)')
ylim([0 5])
set(gca, 'xticklabel', '')
xticks = get(gca, 'xtick');
adjust_space

subplot(312)
plot(Ea_time(IzoomEa), Ea_raw(1,IzoomEa), 'k', 'linewidth', 1)
datetick('x', 15)
xlabel('Oct. 1 2011')
ylabel('E_a (J m^{-2})')
%set(gca, 'xtick', xticks)
xlim([mtime(Izoom(1)) mtime(Izoom(end))])
adjust_space

subplot(313)
plot(sort(A(neap)), Ea_smooth_neap, 'color', [.4 .4 .4])
hold on
plot(sort(A(spring)), Ea_smooth_spring, 'k')
plot(sort(A), Ea_smooth, '--k', 'linewidth', 2)
hold off
xlabel('time to hightide')
ylabel('E_a (J m^{-2})')
xlim([-7 7])
adjust_space

print('-depsc2', 'Eazoom_tide.eps')
% ---------------------------------------------------------------- %


% ---------------- Zoom and M2 cycle (hope the good one) ------------------- %
Izoom = find(mtime>=datenum(2011, 10, 01, 03, 09, 0) & mtime<datenum(2011, 10, 2, 03, 59, 0));
IzoomEa = find(mtime>=datenum(2011, 10, 01, 03, 09, 0) & mtime<datenum(2011, 10, 2, 03, 59, 0));

% W adcp
load('M_N080_vel.mat')
load('M_N080_PTzt.mat')
load('./RBR_mat/019655_79m.mat')

bot = 50;
I = find(z<bot);
w = w(I,:);
Iw = find(time_adcp>=datenum(2011, 09, 21) & time_adcp<datenum(2011, 10, 12));
w = w(:,Iw);
time_adcp = time_adcp(Iw);
IzoomEa = find(time_adcp>=datenum(2011, 10, 01, 03, 09, 0) & time_adcp<datenum(2011, 10, 2, 03, 59, 0));

w_raw = nanmean(w, 1);

% Ea spectrum
%X = detrend(Ea_raw(1,:));
X = detrend(w_raw);
nx = max(size(X));
nx = nx/10;
na = 8;
w = hanning(floor(nx/na));
w = hanning(7200);
[ps, f] = pwelch(X, w, 0, [], fs*86400); 
Ispec = find(f < 1440); % 1min.

% $$$ % Temp spectrum
% $$$ Trbr = RBR.data;
% $$$ t1 = RBR.starttime;
% $$$ t2 = RBR.endtime;
% $$$ t1 = datenum(2011, 09, 19, 16, 0, 0);
% $$$ t2 = datenum(2011, 10, 19, 0, 0, 0);
% $$$ timeRBR = t1:RBR.sampleperiod/86400:t2;
% $$$ Irbr = find(timeRBR>=datenum(2011, 09, 21) & timeRBR<datenum(2011, 10, 12));
% $$$ Trbr = Trbr(Irbr);
% $$$ freqRBR = 1/10;
% $$$ 
% $$$ XX = detrend(Trbr);
% $$$ nx = max(size(XX));
% $$$ nx = nx/10;
% $$$ na = 4;
% $$$ w = hanning(floor(nx/na));
% $$$ [psR, fR] = pwelch(XX, w, 0, [], freqRBR*86400); 
% $$$ loglog(fR, psR, 'k', 'linewidth', 2, 'color', [.4 .4 .4])  

figure(7)
clf
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.1 ; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 20 15])
subplot(2,2,[1 2])
plot(Ea_time(IzoomEa), Ea_raw(1,IzoomEa), 'k', 'linewidth', 1)
hold on
LT1 = datenum(2011, 10, 1, 3, 09, 0);
HT1 = datenum(2011, 10, 1, 9, 25, 0);
LT2 = datenum(2011, 10, 1, 15, 59, 0);
HT2 = datenum(2011, 10, 1, 21, 41, 0);
LT3 = datenum(2011, 10, 2, 3, 59, 0);

patch([LT1 LT1 HT1 HT1 LT1], [0+eps 5e-3-eps 5e-3-eps 0+eps 0+eps], ...
      [1 1 1]*.85, 'edgecolor', 'k')
patch([LT2 LT2 HT2 HT2 LT2], [0+eps 5e-3-eps 5e-3-eps 0+eps 0+eps], ...
      [1 1 1]*.85, 'edgecolor', 'k')
plot(Ea_time(IzoomEa), Ea_raw(1,IzoomEa), 'k', 'linewidth', 1)

datetick('x', 15)
%title('Oct. 1 2011')
ylabel('E_a (J m^{-2})')
%set(gca, 'xtick', xticks)
xlim([mtime(Izoom(1)) mtime(Izoom(end))])
axis on
box on
set(gca, 'tickdir', 'in')

adjust_space
pos1 = get(gca, 'pos');
%pos1(1) = pos(1)-xtra_offset;
adjust_space
pos2 = get(gca, 'pos');
pos1(3) = pos2(1)+pos2(3)-pos1(1);
set(gca, 'pos', pos1)


subplot(223)
%plot(sort(A(neap)), Ea_smooth_neap, 'color', [.4 .4 .4])
%hold on
%plot(sort(A(spring)), Ea_smooth_spring, 'k')
plot(sort(A), Ea_smooth, 'k', 'linewidth', 2)
%hold off
xlabel('Time to hightide (h)')
ylabel('E_a (J m^{-2})')
set(gca, 'xtick', [-6:6]);
set(gca, 'xgrid', 'on')
xlim([-7 7])
adjust_space


subplot(224)
loglog(f(Ispec), ps(Ispec), 'k', 'linewidth', 2)  
hold on
patch([48 48 144 144 48], [2e-9 8e-7 8e-7 2e-9 2e-9], [.85 .85 .85],'edgecolor', 'k')
loglog(f(Ispec), ps(Ispec), 'k', 'linewidth', 2)  

ylabel('PSD (m^2 s^{-2} cpd^{-1})')
xlabel('\sigma (cpd)')
xlim([2e0 2e3])
ylim([1e-9 1e-6])
hold off
adjust_space

print('-depsc2', 'Ea_Wspectrum.eps')
% ----------------------------------------------------------------
% %

keyboard
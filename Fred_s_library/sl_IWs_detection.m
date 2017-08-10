% Basic script, not a function. Must have somewhere the path for a
% xtide output...
% F. Cyr - june 2013
% Example for mooring ADCP 1st file:


clear
load('/home/cyrf0006/PhD/SLEIWEX2013/ADCP_8336a.mat')
tidefile =  '/home/cyrf0006/PhD/SLEIWEX2013/tide_Rimouski2013.dat';

mtime = ADCP.mtime(1213:43990);
w = ADCP.vert_vel(1:65,1213:43990);
u = ADCP.east_vel(1:65,1213:43990);
v = ADCP.north_vel(1:65,1213:43990);
hab = ADCP.config.ranges(1:65);
z = hab;

tbin = 30/86400;


timeVec = mtime(1):tbin:mtime(end);

fid = fopen('./UVW_bin.mat');
if fid == -1 % doesnt exist yet!
             % Time binning
    u_bin = nan(length(hab), length(timeVec));
    v_bin = nan(length(hab), length(timeVec));
    w_bin = nan(length(hab), length(timeVec));
    for i = 1:length(timeVec)
        I = find(mtime >= timeVec(i)-tbin/2 & mtime <= timeVec(i)+tbin/2);
        u_bin(:,i) = nanmean(u(:,I), 2);
        v_bin(:,i) = nanmean(v(:,I), 2);
        w_bin(:,i) = nanmean(w(:,I), 2);
    end
    save UVW_bin.mat u_bin v_bin w_bin timeVec hab z
else    
    load('./UVW_bin.mat');
end



dz = z(2)-z(1);   
rho_0 = 1025;    
Ea = rho_0*dz*nansum((w_bin./100).^2, 1); % in J/m^2 


du = diff(u_bin, 1, 1);
dv = diff(v_bin, 1, 1);
habS2 = hab(1:end-1)+dz;
S2 = (du./dz).^2; + (dv./dz).^2;

% ----- Vertical averages ----- %
I = find(habS2<=20);
S2_20m = nanmean(S2(I,:),1);
I = find(habS2<=5);
S2_5m = nanmean(S2(I,:),1);
S2_U = (du./dz).^2;
S2_V = (dv./dz).^2;  
S2_U = nanmean(S2_U(I,:),1); 
S2_V = nanmean(S2_V(I,:),1);
% ----------------------------- %

% ----- Ea relative to Hightide ----- %
tide  = load(tidefile);
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);
levelFilt = runmean(level, 60);
levelFilt = level;

% find high tide time 
count = 1;
for i = 2:length(mtime)-1

    if levelFilt(i)>=levelFilt(i-1) & levelFilt(i)>=levelFilt(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        if level(i) < .9
            keyboard
        end
        
        count = count+1;
    end
end
% T contains the hour of each high tide
clear mtime

for i = 1:length(timeVec)
    [Y, I] = min(abs(timeVec(i)-T));
    time2(i) = (timeVec(i)-T(I))*24;
    B(i) = L(I); %level of the closest hightide
end

dt = 1;
t1 = -6;
t2 = 6;
t_bin = t1:dt:t2;


% bin the Energy to time vector
Ea_bin = nan(1, length(t_bin));
for i = 1:length(t_bin);
    I = find(time2>=t_bin(i)-dt/2 & time2<t_bin(i)+dt/2);
    Ea_bin(i) = nanmean(Ea(I)); 
end
% ----------------------------------- %


figure(7)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.09 ; % horiz. space between subplots
dy = 0.09; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.02; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %clf


% identify sub-timeseries (24-h)
IHigh = find(timeVec >= datenum(2013, 06, 05, 11, 0, 0) & timeVec <= ...
             datenum(2013, 06, 05, 21, 0, 0));
ILow = find(timeVec >= datenum(2013, 06, 07, 14, 0, 0) & timeVec <= ...
            datenum(2013, 06, 08, 0, 0, 0));


timeLow = timeVec(ILow);
EaLow = Ea(ILow);
S2_5mLow = S2_5m(ILow);
timeHigh = timeVec(IHigh);
EaHigh = Ea(IHigh);
S2_5mHigh = S2_5m(IHigh);

xLow = [timeLow(1) timeLow(end) timeLow(end) timeLow(1) timeLow(1)];
xHigh = [timeHigh(1) timeHigh(end) timeHigh(end) timeHigh(1) timeHigh(1)];
y = [1e-5 1e-5 1 1 1e-5];



% thresholds for Ea and S2
thresh = 1e-3;
thresh = mean(Ea) + 5*std(Ea);



m = nanmean(log(S2_5m));
s2 = nanvar(log(S2_5m));
M = exp(m+s2/2);
VAR = exp(2*m+s2).*(exp(s2)-1);
threshS2 = M+5*sqrt(VAR);
% $$$ threshS2 = 1e-2;
% $$$ threshS2 = nanmean(S2_5m)*10;

LETTER_Y= 4.5e-1; 
GRAY_LEV = .4;

% sub1
sub1 = subplot(3,2,[1 2]);
semilogy(timeVec, S2_5m, 'k')   
hold on
patch(xLow, y, [1 1 1]*.8)
patch(xHigh, y, [1 1 1]*.8)
semilogy(timeVec, S2_5m, 'k') 
semilogy(timeVec, Ea, 'color', [1 1 1]*GRAY_LEV)
plot([timeVec(1) timeVec(end)], [thresh thresh], 'color', [1 1 1]*GRAY_LEV, 'linestyle', '--')
plot([timeVec(1) timeVec(end)], [threshS2 threshS2], '--k')
hold off
ylim([1e-5 1])
datetick('x', 7)
xlim([timeVec(1) timeVec(end)])
xlabel('Sept/Oct 2011')
ylabel('E_a, S^2 (Kg s^{-2}, s^{-2})')
set(gca, 'ytick', [1e-4 1e-3 1e-2 1e-1])
set(gca, 'ygrid', 'on')
text(timeVec(end)-0.5, LETTER_Y, 'a', 'fontsize', 11, 'fontweight', 'bold')
adjust_space
pos1 = get(sub1, 'pos');

adjust_space
pos2 = get(sub1, 'pos');
offset = pos2(1)-pos1(1)-pos1(3);
pos1(3) = 2*pos1(3)+offset;
set(sub1, 'pos', pos1)

% sub2
sub2 = subplot(323);
semilogy(timeHigh, S2_5mHigh, 'k') 
hold on
semilogy(timeHigh, EaHigh,  'color', [1 1 1]*GRAY_LEV)
plot([timeHigh(1) timeHigh(end)], [thresh thresh], 'color', [1 1 1]*GRAY_LEV, 'linestyle', '--')
plot([timeHigh(1) timeHigh(end)], [threshS2 threshS2], '--k')
hold off
ylim([1e-5 1])
datetick
xlim([timeHigh(1) timeHigh(end)])
xlabel(datestr(timeHigh(1), 1))
ylabel('E_a, S^2 (kg s^{-2}, s^{-2})')
set(gca, 'ytick', [1e-4 1e-3 1e-2 1e-1])
set(gca, 'ygrid', 'on')
set(gca, 'xtick', datenum(2011, 09, 22, 12:4:20, 0, 0))
text(timeHigh(end)-0.03, LETTER_Y, 'b', 'fontsize', 11, 'fontweight', 'bold')
adjust_space

% sub3
sub3 = subplot(3,2,4);
semilogy(timeLow, S2_5mLow, 'k') 
hold on
semilogy(timeLow, EaLow,   'color', [1 1 1]*GRAY_LEV)
plot([timeLow(1) timeLow(end)], [thresh thresh], 'color', [1 1 1]*GRAY_LEV, 'linestyle', '--')
plot([timeLow(1) timeLow(end)], [threshS2 threshS2], '--k')
hold off
ylim([1e-5 1])
datetick
xlim([timeLow(1) timeLow(end)])
xlabel(datestr(timeLow(1), 1))
set(gca, 'ytick', [1e-4 1e-3 1e-2 1e-1])
set(gca, 'yticklabel', [])
set(gca, 'ygrid', 'on')
set(gca, 'xtick', datenum(2011, 10, 01, 15:4:23, 0, 0))
text(timeLow(end)-0.03, LETTER_Y, 'c', 'fontsize', 11, 'fontweight', 'bold')
adjust_space


% sub4
sub4 = subplot(3,2,5);

I = find(Ea>=thresh);
J = find(S2_5m(I)>=threshS2);

semilogy(Ea(I), S2_5m(I), '.k')
hold on
plot([min(Ea(I)) max(Ea(I))],[M M], 'k')
plot([min(Ea(I)) max(Ea(I))],[threshS2 threshS2], '--k')
hold off
xlim([min(Ea(I)) max(Ea(I))])
ylim([5e-4 1])
ylabel('S^2 (s^{-2})')
xlabel('E_a (kg s^{-2})')
set(gca, 'ytick', [1e-3 1e-2 1e-1])
set(gca, 'ygrid', 'on')
set(gca, 'yminorgrid', 'off')

disp(sprintf('%d out of %d IWs enhance shear (%d percent)', length(J), length(I), round(length(J)/length(I)*100)));
text(7.2e-3, LETTER_Y, 'd', 'fontsize', 11, 'fontweight', 'bold')
adjust_space


% sub5
sub5 = subplot(3,2,6);
% $$$ plot(t_bin, Ea_bin, 'k')
% $$$ xlabel('time to hightide (h)')
% $$$ ylabel('E_a (J m^{-2})')
% $$$ xlim([-6.5 6.5])

clear relativeNo;
I = find(Ea>=thresh);
dclass = 1;
classEdge = -7:dclass:7;
for j = 1:length(classEdge)-1
    J = find(time2(I)>=classEdge(j) & time2(I)<classEdge(j+1));
    relativeNo(j) = length(J)/length(I);
end

bar(classEdge(1:end-1)+dclass/2, relativeNo, 'facecolor', [.2 .2 .2]);
xlim([-7 7])
%ylim([0 .35])
ylabel('rel. freq.')
xlabel('time to hightide (h)')
II = find(S2_5m(I)>=threshS2);
clear relativeNo;

for j = 1:length(classEdge)-1
    J = find(time2(I(II))>=classEdge(j) & time2(I(II))<classEdge(j+1));
    relativeNo(j) = length(J)/length(I);
end
hold on
bar(classEdge(1:end-1)+dclass/2, relativeNo, 'facecolor', [.2 .2 .2]*2);
hold off
text(6.35, .32, 'e', 'fontsize', 11, 'fontweight', 'bold')
adjust_space


print(gcf, '-dpng', 'Ea_S2_subplot.png')
set(gcf, 'renderer', 'painters')
print(gcf, '-depsc2', 'Ea_S2_subplot.eps')





Tadcp = ADCP.temperature(2000:41000);
Td = detrend(Tadcp);

% Spectral analysis
I = find(z>=10 & z<=30);
ww = nanmean(w(I,:));
freq = 1/10;

Wd = detrend(ww);


nx = max(size(Wd)); % Hanning
na = 10;
h = hanning(floor(nx/na));
[ps, f] = pwelch(Wd, h, 0, [], freq*86400); 

nx = max(360*6); % Hanning
na = 10;
h = hanning(floor(nx/na));
[psT, fT] = pwelch(Td, h, 0, [], freq*86400); 

figure(10)
clf
loglog(f, ps, 'r', 'linewidth', 2)  
hold on
loglog(fT, psT, 'b', 'linewidth', 2)  

ylabel('PSD')
%xlabel('f (day^{-1})')
%xlim([1e-1 3e2])
%ylim([1e-8 1e0])


plot_harmonics3
hold off



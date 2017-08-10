clear


% exec in /home/cyrf0006/PhD/hypox/mouillage


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

% [P T S sigT O2]
data = load('mooring_hypox_data_2008-2010.dat');

fid = fopen('mooring_hypox_datefile_2008-2010.dat');
C = textscan(fid, '%s', 'delimiter', '\n');
date = char(C{1});


mtime = datenum(date);
time_reg = mtime(1):1/24:mtime(end);
orig_time = mtime;

O2 = data(:,5);
P = data(:,1);
T = data(:,2);
S = data(:,3);

% sort profiles
[mtime, I] = sort(orig_time);
O2 = O2(I);
P = P(I);
T = T(I);
S = S(I);

% $$$ % ml/l to mg/l
% $$$ O2 = O2.*1.42903;
% $$$ % mg/l to umol/l
% $$$ O2 = O2.*1000/32;



% remove bad values and interpolate
I = find(O2==-99);
O2(I) = [];
t = mtime;
t(I) = [];
O2 = interp1(t, O2, time_reg);

I = find(T==-99);
T(I) = [];
t = mtime;
t(I) = [];
T = interp1(t, T, time_reg);

I = find(S==-99);
S(I) = [];
t = mtime;
t(I) = [];
S = interp1(t, S, time_reg);

I = find(P==-99);
P(I) = [];
t = mtime;
t(I) = [];
P = interp1(t, P, time_reg);


mtime = time_reg;
timelim = [mtime(1) mtime(end)];
%xtickmark = timelim(1):timelim(2);




figure(1) % Interpolated fields
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 20])

% convert O2
% ml/l to mg/l
O2 = O2.*1.42903;
% mg/l to umol/l
O2 = O2.*1000/32;

subplot(3,1,1)
plot(mtime, P)
%ylim([5 100])
ylabel('P (dBar)')
xlim(timelim)
datetick('x', 12)
xlim(timelim)
%set(gca, 'xtick', xtickmark)
set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
%set(gca, 'ygrid', 'on')
title('Hypoxie Mooring')
adjust_space

subplot(3,1,2)
[AX, H1, H2] = plotyy(mtime, T, mtime, S)
set(get(AX(1),'Ylabel'),'String','T(^{\circ}C)')
set(get(AX(2),'Ylabel'),'String','S') 
xlim(AX(1),[timelim(1) timelim(2)]);
xlim(AX(2),[timelim(1) timelim(2)]);
%ylim(AX(1),[0 3]);
%ylim(AX(2),[60 70]);
%set(AX(1),'YTick',[0 1 2 3])
%set(AX(2),'YTick',[60 62 64 66 68 70])
%set(AX(1),'ygrid','on')
%set(AX(2),'ygrid','on')
%set(AX, 'xTick',xtickmark) % delete the x-tick labels
datetick('x', 12)
xlim(timelim)
set(AX, 'xTickLabel','') % delete the x-tick labels
set(gca, 'xgrid', 'on')

adjust_space

subplot(3,1,3)
plot(mtime, O2)
ylabel('dissol. oxyg. (\mu mol/l)')
%ylim([5 100])
xlim(timelim)
set(gca, 'xgrid', 'on')
datetick('x', 12)
xlim(timelim)
X = runmean(O2, 360);
X = runmean(O2, 180);
hold on
plot(mtime, X, 'r')


adjust_space
% ----------------------- %



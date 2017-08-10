clear


drange = [60 80];



% ----------- M_N080 ------------ %
load('../M_N080/M_N080_vel.mat')
load('../M_N080/M_N080_PTzt.mat')
z_adcp = 59; %m; found from look at Pressure



% vertical averages
I = find(z+z_adcp > drange(1) & z+z_adcp < drange(2));
E600 = nanmean(east_vel(I, :), 1);
N600 = nanmean(north_vel(I, :), 1);


% time average
dt = 1/48; %(time in days, 1/48 is 30 minutes...)
ave_time600 = round(time_adcp(1)):dt:round(time_adcp(end));

for i = 1:length(ave_time600)
    I = find(time_adcp > ave_time600(i)-dt/2 &  time_adcp < ave_time600(i)+dt/2);
    ave_E600(i) = nanmean(E600(I));
    ave_N600(i) = nanmean(N600(I));
end

% -------------------------------- %


% ----------- M_RIKI ------------ %
load('../M_RIKI/M_RIKI_vel.mat')
load('../M_RIKI/M_RIKI_PTzt.mat')
z_adcp = 149; %m; found from look at Pressure


% vertical averages
I = find(z_adcp-z > drange(1) & z_adcp-z < drange(2));
E300 = nanmean(east_vel(I, :), 1);
N300 = nanmean(north_vel(I, :), 1);


% $$$ % correct time_adcp 
% $$$ disp('correct time for 2011')
% $$$ tocmp(1:length(time_adcp), 1) = '0';
% $$$ tocmp(1:length(time_adcp), 2) = '0';
% $$$ tocmp(1:length(time_adcp), 3) = '1';
% $$$ tocmp(1:length(time_adcp), 4) = '1';
% $$$ tocmp = cellstr(tocmp);
% $$$ 
% $$$ A = datestr(time_adcp, 31);
% $$$ AA = cellstr(A(:,1:4));
% $$$ 
% $$$ bo = strcmp(tocmp, AA);
% $$$ 
% $$$ I = find(bo==1);
% $$$ A(I,1) = '2';
% $$$ time_adcp = datenum(A);


% time average
dt = 1/48; %(time in days, 1/48 is 30 minutes...)
ave_time300 = round(time_adcp(1)):dt:round(time_adcp(end));



for i = 1:length(ave_time300)
    I = find(time_adcp > ave_time300(i)-dt/2 &  time_adcp < ave_time300(i)+dt/2);
    ave_E300(i) = nanmean(E300(I));
    ave_N300(i) = nanmean(N300(I));
end


% -------------------------------- %

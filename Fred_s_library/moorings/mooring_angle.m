function mooring_angle(tidefile, rbrfile, hab, timelim)
    
% ex. for M_N080, RBR at ~30m, i.e., around 50m above bottom:
% 
% mooring_angle('wind_tides/tide_moorings.dat','M_N080/RBR_mat/014766_30m.mat' , 50, [datenum(2011, 09, 20 ,0, 0 ,0) datenum(2011, 10, 12, 12, 0, 0)])
% 
% ex: for M_RIKI at ~60m (executer dans M_RIKI)
% 
% mooring_angle('../wind_tides/tide_moorings.dat','SBE_Mriki_60m.dat' , 270, [datenum(2011, 09, 20 ,0, 0 ,0) datenum(2011, 10, 12, 12, 0, 0)])
%

% ----------------------------- Tides ---------------------------- %
data = load(tidefile);

time_tide = datenum(data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), 0);
tide = data(:,6);
clear data
% --------------------------------------------------------------- %


% -------------------------- Thermistors ------------------------ %
load(rbrfile);

t0 = datenum([RBR.starttime(7:10) '/' RBR.starttime(4:5) '/'  RBR.starttime(1:2) ' ' RBR.starttime(12:22)]);
tf = datenum([RBR.endtime(7:10) '/' RBR.endtime(4:5) '/'  RBR.endtime(1:2) ' ' RBR.endtime(12:22)]);
dt = RBR.sampleperiod/3600/24;
time_rbr = [datenum(t0:dt:tf)];
time_rbr = time_rbr(1:size(RBR.data,1));

data = RBR.data;
T_rbr = data(:,1);
P_rbr = data(:,2);
Z_rbr = data(:,3);

clear data RBR
% --------------------------------------------------------------- %

% $$$ % -------------------------- SeaBird CTD ------------------------ %
% $$$ data = load(rbrfile);
% $$$ datefile = ['date_' rbrfile];
% $$$ fid = fopen(datefile);
% $$$ C = textscan(fid, '%s', 'delimiter', '\n');
% $$$ date = char(C{1});
% $$$ siz = size(date,1);
% $$$ 
% $$$ MM=ones(siz, 1)*9;
% $$$ for i = 1:siz
% $$$     if date(i, 3:5)=='Oct'
% $$$         MM(i)=10;
% $$$     end
% $$$ end
% $$$ 
% $$$ YYYY = str2num(date(:,6:9));
% $$$ DD = str2num(date(:,1:2));
% $$$ hh = str2num(date(:,11:12));
% $$$ mm = str2num(date(:,14:15));
% $$$ ss = str2num(date(:,17:18));
% $$$ 
% $$$ % we keep RBR notation for simplicity
% $$$ time_rbr = datenum(YYYY, MM, DD, hh, mm, ss);
% $$$ 
% $$$ T_rbr = data(:,2)';
% $$$ C_rbr = data(:,3)';
% $$$ Z_rbr = data(:,1)';
% $$$ 
% $$$ clear data date
% $$$ % --------------------------------------------------------------- %

% itp tide to RBR resolution
tide_itp = interp1(time_tide, tide, time_rbr);

% restrict range
I = find(time_rbr>timelim(1) & time_rbr<timelim(2));

mean_tide = mean(tide_itp(I));
rbr_depth = mean(Z_rbr(I));
rbr_displ = [rbr_depth-Z_rbr(I)]';
tide_displ = [mean_tide-tide_itp(I)];

keyboard
plot(time_rbr(I), rbr_displ)                            
hold on
plot(time_rbr(I), tide_displ, 'r')  
plot(time_rbr(I), rbr_displ-tide_displ, 'k')

legend('rbr displ. rel. to equil.', 'tide displ. rel. to equil.', 'difference')
datetick('x', 7)
xlabel('date (sept./oct.)')
ylabel('height(m)')



keyboard
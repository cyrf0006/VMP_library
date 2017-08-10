clear


% Define grid
t1 = datenum(2011, 09, 20, 0, 0, 0);
%t1 = datenum(2011, 10, 12, 0, 0, 0);
t2 = datenum(2011, 10, 12, 0, 0, 0);
dt = 10/60/24;
%dt = 10/86400;

zvec = 50:1:150;
time_vec = [t1:dt:t2];


% Load every thermistor data and average to required resolution
disp('RBR...')
T_rbr70m = getrbr('RBR_mat/019599_70m.mat', time_vec);
T_rbr80m = getrbr('RBR_mat/019600_80m.mat', time_vec);
T_rbr90m = getrbr('RBR_mat/019604_90m.mat', time_vec);



% load SBE
disp('SBE...')
SBE = load('SBE_Mriki_60m.dat');
P_sbe = SBE(:,1);
T_sbe = SBE(:,2);
C_sbe = SBE(:,3);
%S_sbe = 
fid = fopen('date_SBE_Mriki_60m.dat');
C = textscan(fid, '%s', 'delimiter', '\n');
datesbe = char(C{1});
t1 = datenum(datesbe(1,:));
t2 = datenum(datesbe(2,:));
te = datenum(datesbe(end,:));
time_sbe = t1:t2-t1:te;

% load minilog
disp('MINILOG...')
MINI = load('mini_Mriki_100m.dat');
T_mini = MINI;
fid = fopen('date_mini_Mriki_100m.dat');
C = textscan(fid, '%s', 'delimiter', '\n');
datemini = char(C{1});
t1 = datenum(datemini(1,:));
t2 = datenum(datemini(2,:));
te = datenum(datemini(end,:));
time_mini = t1:t2-t1:te;

% ADCP
disp('ADCP...')
load M_RIKI_PTzt.mat
P_adcp = P_adcp/1000;

% tide file for depth correction
tide  = load('../tide_shear/tide_2009-2011.dat');
time_tide = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);


% bin all field to required resolution
for i = 1:length(time_vec)
    % ADCP
    I = find(time_adcp >= time_vec(i)-dt/2 & time_adcp <= time_vec(i)+dt/2);
    T_adcp150m(i) = nanmean(T_adcp(I));
    P_adcp150m(i) = nanmean(P_adcp(I));
    
    %SBE
    I = find(time_sbe >= time_vec(i)-dt/2 & time_sbe <= time_vec(i)+dt/2);
    T_sbe60m(i) = nanmean(T_sbe(I));
    C_sbe60m(i) = nanmean(C_sbe(I));
    P_sbe60m(i) = nanmean(P_sbe(I));
    
    %MINI
    I = find(time_mini >= time_vec(i)-dt/2 & time_mini <= time_vec(i)+dt/2);
    T_mini100m(i) = nanmean(T_mini(I));
    
end

% tide (interpolation from hourly)
tidelevel = interp1(time_tide, level, time_vec);

% ---- Correct vertical position for displacements ---- %
zcorr = tidelevel - nanmean(tidelevel);

Z_adcp = P_adcp150m-zcorr;
Z_sbe = P_sbe60m-zcorr;

% real target depth
pos_sbe = min(Z_sbe);
pos_adcp = min(Z_adcp);

pos_rbr70m = (pos_adcp - pos_sbe)./(150-60)*10+pos_sbe;
pos_rbr80m = (pos_adcp - pos_sbe)./(150-60)*20+pos_sbe;
pos_rbr90m = (pos_adcp - pos_sbe)./(150-60)*30+pos_sbe;
pos_mini100m = (pos_adcp - pos_sbe)./(150-60)*40+pos_sbe;


adcp_displ = Z_adcp - pos_adcp;
sbe_displ = Z_sbe - pos_sbe;
rbr70m_displ = nan(size(sbe_displ));
rbr80m_displ = nan(size(sbe_displ));
rbr90m_displ = nan(size(sbe_displ));
mini100m_displ = nan(size(sbe_displ));
I = find(~isnan(adcp_displ)==1);
rbr70m_displ(I) = sbe_displ(I) - (sbe_displ(I) - adcp_displ(I))./(150-60)*10;
rbr80m_displ(I) = sbe_displ(I) - (sbe_displ(I) - adcp_displ(I))./(150-60)*20;
rbr90m_displ(I) = sbe_displ(I) - (sbe_displ(I) - adcp_displ(I))./(150-60)*30;
mini100m_displ(I) = sbe_displ(I) - (sbe_displ(I) - adcp_displ(I))./(150-60)*40;
I = find(isnan(adcp_displ)==1);
rbr70m_displ(I) = sbe_displ(I) - (nanmean(sbe_displ) - nanmean(adcp_displ))./(150-60)*10;
rbr80m_displ(I) = sbe_displ(I) - (nanmean(sbe_displ) - nanmean(adcp_displ))./(150-60)*20;
rbr90m_displ(I) = sbe_displ(I) - (nanmean(sbe_displ) - nanmean(adcp_displ))./(150-60)*30;
rbr100m_displ(I) = sbe_displ(I) - (nanmean(sbe_displ) - nanmean(adcp_displ))./(150-60)*40;


Z_rbr70m = pos_rbr70m + rbr70m_displ;
Z_rbr80m = pos_rbr80m + rbr80m_displ;
Z_rbr90m = pos_rbr90m + rbr90m_displ;
Z_mini100m = pos_mini100m + mini100m_displ;
% ----------------------------------------------------- %




% ---- 1st method: vertical interp1 ---- %
% low-pass filtering
[b, a] = butter(4, .4/144, 'low');
T_sbe60mf = filtfilt(b, a, T_sbe60m);
T_rbr70mf = filtfilt(b, a, T_rbr70m);
T_rbr80mf = filtfilt(b, a, T_rbr80m);
T_rbr90mf = filtfilt(b, a, T_rbr90m);
T_mini100mf = filtfilt(b, a, T_mini100m);
T_adcp150mf = filtfilt(b, a, T_adcp150m);
% put them togehter (ungrid)
Z_raw = [Z_sbe; Z_rbr70m; Z_rbr80m; Z_rbr90m; Z_mini100m; Z_adcp];
T_raw = [T_sbe60m; T_rbr70m; T_rbr80m; T_rbr90m; T_mini100m; T_adcp150m];
%T_raw = [T_sbe60mf; T_rbr70mf; T_rbr80mf; T_rbr90mf; T_mini100mf; T_adcp150mf];


% vertical interp1
disp('vertical interpolation')
Tgrid = nan(length(zvec), length(time_vec));
for i = 1:length(Z_sbe)
    I = find(~isnan(T_raw(:,i))==1);
    if length(I) > 1
        if sum(isnan(Z_raw(I,i)))==0; % no Nan in Z 
            Tgrid(:,i) = interp1(Z_raw(I,i), T_raw(I,i), zvec);
        end
    end
end
% -------------------------------------- %

% $$$ % ---- 2nd method: griddata ---- %
% $$$ % vertical interp1
% $$$ disp('2D-interpolation')
% $$$ % low-pass filtering
% $$$ [b, a] = butter(4, .4/144, 'low');
% $$$ T_sbe60mf = filtfilt(b, a, T_sbe60m);
% $$$ T_rbr70mf = filtfilt(b, a, T_rbr70m);
% $$$ T_rbr80mf = filtfilt(b, a, T_rbr80m);
% $$$ T_rbr90mf = filtfilt(b, a, T_rbr90m);
% $$$ T_mini100mf = filtfilt(b, a, T_mini100m);
% $$$ T_adcp150mf = filtfilt(b, a, T_adcp150m);
% $$$ 
% $$$ % put them togehter (ungrid)
% $$$ Z_raw = [Z_sbe Z_rbr70m Z_rbr80m Z_rbr90m Z_mini100m Z_adcp];
% $$$ T_raw = [T_sbe60m T_rbr70m T_rbr80m T_rbr90m  T_mini100m T_adcp150m];
% $$$ %T_raw = [T_sbe60mf T_rbr70mf T_rbr80mf T_rbr90mf T_mini100mf T_adcp150mf];
% $$$ time_raw = [time_vec time_vec time_vec time_vec time_vec time_vec];
% $$$ % remove NaNs
% $$$ I = find(isnan(Z_raw)==1);
% $$$ Z_raw(I) = [];
% $$$ T_raw(I) = [];
% $$$ time_raw(I) = [];
% $$$ I = find(isnan(T_raw)==1);
% $$$ Z_raw(I) = [];
% $$$ T_raw(I) = [];
% $$$ time_raw(I) = [];
% $$$ 
% $$$ % griddata
% $$$ [XI, YI] = meshgrid(time_vec, zvec);
% $$$ Tgrid = griddata(time_raw, Z_raw, T_raw, XI, YI, 'linear');
% $$$ % -------------------------------------- %



imagesc(time_vec, zvec, Tgrid)
set(gca, 'ydir', 'reverse')
colorbar
caxis([0 4])
datetick('x', 7)
xlabel('sept. / oct.')
ylabel('z (m)')


keyboard
% change name to avoid confusion
time_vec_RIKI = time_vec;
%save Tfield_M_RIKI.mat Tgrid zvec time_vec_RIKI
save Tfield_M_RIKIcubic.mat Tgrid zvec time_vec_RIKI

save Tseries_M_RIKI.mat T_sbe60m  T_rbr70m  T_rbr80m  T_rbr90m ...
    T_mini100m T_adcp150m tidelevel time_vec 
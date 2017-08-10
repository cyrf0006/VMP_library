clear


% Define grid
t1 = datenum(2011, 09, 20, 0, 0, 0);
t2 = datenum(2011, 10, 12, 0, 0, 0);
dt = 10/60/24;

hab = 0:50;
time_vec = [t1:dt:t2];


% Load every thermistor data and average to required resolution
T_rbr30m = getrbr('RBR_mat/014766_30m.mat', time_vec);
T_rbr40m = getrbr('RBR_mat/014765_40m.mat', time_vec);
T_rbr50m = getrbr('RBR_mat/014764_50m.mat', time_vec);
T_rbr60m = getrbr('RBR_mat/019603_60m.mat', time_vec);
T_rbr65m = getrbr('RBR_mat/019589_65m.mat', time_vec);
T_rbr70m = getrbr('RBR_mat/019602_70m.mat', time_vec);
T_rbr75m = getrbr('RBR_mat/019652_75m.mat', time_vec);
T_rbr79m = getrbr('RBR_mat/019655_79m.mat', time_vec);


% load SBE
% NO SBE!

% ADCP
load M_N080_PTzt.mat

% to required resolution
for i = 1:length(time_vec)
    I = find(time_adcp >= time_vec(i)-dt/2 & time_adcp <= time_vec(i)+dt/2);
    T_adcp55m(i) = nanmean(T_adcp(I));
end


% put them togehter (ungrid)
hab_raw = [50 40 30 25 20 15 10 5 1];
Tgrid_raw = [T_rbr30m; T_rbr40m; T_rbr50m; T_adcp55m; T_rbr60m; T_rbr65m; ...
             T_rbr70m; T_rbr75m;  T_rbr79m];


% Griddata!
[XI, YI] = meshgrid(time_vec, hab);
Tgrid = griddata(time_vec, hab_raw, Tgrid_raw, XI, YI, 'linear');


% display

imagesc(time_vec, hab, Tgrid)
set(gca, 'ydir', 'normal')
colorbar
caxis([0 4])
datetick('x', 7)
xlabel('sept. / oct.')
ylabel('hab (m)')


% change name to avoid confusion
time_vec_N080 = time_vec;
save Tfield_M_N080.mat Tgrid hab time_vec_N080



% added for cubic interp
Tgrid = griddata(time_vec, hab_raw, Tgrid_raw, XI, YI, 'cubic');
save Tfield_M_N080cubic.mat Tgrid hab time_vec_N080
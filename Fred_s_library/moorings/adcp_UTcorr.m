% Find best shift for correlation between T, U and V
clear


% The beginning is similar to tide_shear_MS.m
fid = fopen('/home/cyrf0006/LaTeX/MS/BMix/matlab_files/spectrum_data.mat');
if fid == -1 % doesnt exist yet!
    disp('please run ''mooring_spectrum.m'' first')
else  
   load('/home/cyrf0006/LaTeX/MS/BMix/matlab_files/spectrum_data.mat');
end

% ----------------------------- Spectral Analysis ---------------------------- %


freq = 1/(timeVec(2)-timeVec(1))/86400; % Hz

% TEST !!!!!!
%S2_bin = sqrt(S2_bin);

% last 20m
I = find(zhab<=20);
u0_20m = nanmean(u_bin(I, :));
v0_20m = nanmean(v_bin(I, :));
U0_20m = sqrt(u0_20m.^2 + v0_20m.^2); 
I = find(zhab_S2<=20);
S0_20m = nanmean(S2_bin(I, :));
Su0_20m = nanmean(S2_bin_u(I, :));
Sv0_20m = nanmean(S2_bin_v(I, :));

% last 10m
I = find(zhab<=10);
u0_10m = nanmean(u_bin(I, :));
v0_10m = nanmean(v_bin(I, :));
U0_10m = sqrt(u0_10m.^2 + v0_10m.^2); 
I = find(zhab_S2<=10);
S0_10m = nanmean(S2_bin(I, :));
Su0_10m = nanmean(S2_bin_u(I, :));
Sv0_10m = nanmean(S2_bin_v(I, :));

%last 5m
I = find(zhab<=5);
u0_5m = nanmean(u_bin(I, :));
v0_5m = nanmean(v_bin(I, :));
U0_5m = sqrt(u0_5m.^2 + v0_5m.^2); 
I = find(zhab_S2<=5);
S0_5m = nanmean(S2_bin(I, :));
Su0_5m = nanmean(S2_bin_u(I, :));
Sv0_5m = nanmean(S2_bin_v(I, :));

% 10m hab
[Y, I] = min(abs(10-zhab));
u10m = u_bin(I, :);
v10m = v_bin(I, :);
U10m = sqrt(u10m.^2 + v10m.^2); 
[Y, I] = min(abs(10-zhab_S2));
S10m = S2_bin(I, :);

% 10-20m
I = find(zhab<=20 & zhab>=10);
u10_20m = nanmean(u_bin(I, :));
v10_20m = nanmean(v_bin(I, :));
U10_20m = sqrt(u0_20m.^2 + v0_20m.^2); 
I = find(zhab_S2<=20 & zhab_S2>=10);
S10_20m = nanmean(S2_bin(I, :));
Su10_20m = nanmean(S2_bin_u(I, :));
Sv10_20m = nanmean(S2_bin_v(I, :));
% ----------------------------------- %

%  --------- Temp for M_N080 --------------- %
load /home/cyrf0006/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/Tfield_M_N080_1min.mat
dt = (time_vec_N080(2)-time_vec_N080(1))*86400;

% bin T to velocities
dtVel = timeVec(2)-timeVec(1);
TMat = nan(length(hab), length(timeVec)); 
for i = 1:length(timeVec)
    I = find(time_vec_N080>=timeVec(i)-dtVel & time_vec_N080<timeVec(i)+dtVel);    
    TMat(:,i) = nanmean(Tgrid(:,I), 2);
end
    
Tgrid = TMat;
T10m = Tgrid(11,:); % 10m hab
I = find(hab<=5); 
T0_5m = nanmean(Tgrid(I,:),1); % last 10m
I = find(hab<=20);
T0_10m = nanmean(Tgrid(I,:),1); % last 10m
I = find(hab<=20);
T0_20m = nanmean(Tgrid(I,:),1); % last 20m
I = find(hab<=20 & hab>=10);
T10_20m = nanmean(Tgrid(I,:),1); % last 10-20m
I = find(hab<=50 & hab>=40);
T40_50m = nanmean(Tgrid(I,:),1); % last 10-20m% ----------------------------------- %



% ------------- X-Covariance -------------- %
MAX_LAG = 60*6; %minutes
dt = dtVel*86400/60; % in minutes
maxLag = 72; %int16(MAX_LAG/dt);%


% $$$ % matlab ex:
% $$$ ww = randn(1000,1);     % White Gaussian noise
% $$$ [cov_ww,lags] = xcov(ww,10,'coeff');
% $$$ stem(lags,cov_ww)



% $$$ [cU,lagsU] = xcov(T0_10m, u0_10m, 72, 'coef');
% $$$ [cV,lagsV] = xcov(T0_10m, v0_10m, 72,'coef');
% $$$ [cU2,lagsU2] = xcov(T10_20m, u10_20m, 72,'coef');
% $$$ [cV2,lagsV2] = xcov(T10_20m, v10_20m, 72,'coef');
[cU,lagsU] = xcov(T0_10m, u0_10m, 72, 'coef');
[cV,lagsV] = xcov(T0_10m, v0_10m, 72,'coef');
[cU2,lagsU2] = xcov(T40_50m, u10_20m, 72,'coef');
[cV2,lagsV2] = xcov(T40_50m, v10_20m, 72,'coef');


figure(4)
clf
subplot(211)
stem(lagsU*5/60, cU)
hold on
stem(lagsV*5/60, cV, 'r')
title('0-10m')
legend('u','v')
xlabel('lag (hour)')
ylabel('corr coef (T vs U,V)')

subplot(212)
stem(lagsU2*5/60, cU2)
hold on
stem(lagsV2*5/60, cV2, 'r')
title('10-20m / 40-50m')
legend('u','v')
xlabel('lag (hour)')
ylabel('corr coef (T vs U,V)')


disp('MAX CORR')
[c1, I] = max(cU);
[c2, J] = max(cV);
disp('u,v for 0-10m (hour):')
[lagsU(I) lagsV(J)].*dt/60
disp('corresponding coef:')
[cU(I) cV(J)]
disp('')
[c1, I] = max(cU2);
[c2, J] = max(cV2);
disp('u,v for 10-20m / T for 40-50m(hour):')
[lagsU2(I) lagsV2(J)]*dt/60
disp('corresponding coef:')
[cU2(I) cV2(J)]

disp(' ')
disp(' ')

disp('MIN CORR (ANTI-CORR)')
[c1, I] = min(cU);
[c2, J] = min(cV);
disp('u,v for 0-10m (hour):')
[lagsU(I) lagsV(J)].*dt/60
disp('corresponding coef:')
[cU(I) cV(J)]
disp(' ')
[c1, I] = min(cU2);
[c2, J] = min(cV2);
disp('u,v for 10-20m (hour):')
[lagsU2(I) lagsV2(J)]*dt/60
disp('corresponding coef:')
[cU2(I) cV2(J)]

keyboard

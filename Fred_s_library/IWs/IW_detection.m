function IW = IW_detection(vel_file, time_file, thresh)
    
% ex: >> IW = IW_detection('M_N080_vel.mat', 'M_N080_PTzt.mat', 1e-3)  
        
load(vel_file);
load(time_file);

% params
rho_0 = 1025;

% find bottom ()
for i = 1:size(w, 1)
    I = find(isnan(w(i, :))==1);
    if length(I)>size(w, 2)./4
        bot = i;
        break
    end
end

bottomDepth = z(bot);

% reduce file to desired depth range (1.5m removed!!!)
w = w(1:bot-3, :);
z = z(1:bot-3);
u = east_vel(1:bot-3, :);
v = north_vel(1:bot-3, :);
dz = z(2)-z(1);

% reduce time range
tmin = datenum(2011, 09, 20, 0, 0, 0);
tmax = datenum(2011, 10, 12, 13, 0, 0);
I = find(time_adcp>tmin & time_adcp<tmax);
time_adcp = time_adcp(I);
w = w(:,I);
u = u(:,I);
v = v(:,I);


w_itp = nan(size(w));
u_itp = nan(size(w));
v_itp = nan(size(w));
% interpolation to remove NaNs
for i = 1:size(w, 1)
    I = find(isnan(w(i, :))==1);
    if ~isempty(I) == 1 
        y = w(i, :);
        x = 1:size(w,2);
        y(I) = [];
        x(I) = [];
        w_itp(i, :) = interp1(x, y, 1:size(w, 2));
        
        y = u(i, :);
        x = 1:size(u,2);
        y(I) = [];
        x(I) = [];
        u_itp(i, :) = interp1(x, y, 1:size(u, 2));
        
        y = v(i, :);
        x = 1:size(v,2);
        y(I) = [];
        x(I) = []; 
        v_itp(i, :) = interp1(x, y, 1:size(v, 2));
    else
        w_itp(i, :) = w(i, :);
        u_itp(i, :) = u(i, :);
        v_itp(i, :) = v(i, :);
    end
end

% Save clean velocities for Ea S2 comparison
hab = bottomDepth-z;
save cleanUVW.mat u_itp v_itp w_itp time_adcp z hab


% band pass filter (tides and noise)
dt_adcp = (time_adcp(2) - time_adcp(1)).*86400;
fs = 1/dt_adcp;
freq_low = 1/90; %Hz
freq_high = 1/1800;
Wn_high = freq_high/(fs/2);
Wn_low = freq_low/(fs/2);

%[b,a] = butter(4,[Wn_low Wn_high]);
[b,a] = butter(4,[Wn_high Wn_low]);

disp('filering data ...')
w_filt = nan(size(w_itp));
for i = 1:size(w_itp, 1)
    w_filt(i, :) = filtfilt(b, a, w_itp(i, :));
    u_filt(i, :) = filtfilt(b, a, u_itp(i, :));
    v_filt(i, :) = filtfilt(b, a, v_itp(i, :));
end
disp('done!')



W_itp = nanmean(w_itp, 1);
W_filt = nanmean(w_filt, 1);


% Compute shear...
du = diff(u_filt, 1, 1);
dv = diff(v_filt, 1, 1);
zS2 = z(1:end-1)+dz;
S2 = (du./dz).^2; + (dv./dz).^2;



% Area energy density
Ea = rho_0*dz*nansum((w_filt./100).^2, 1); % in J/m^2 

plot(time_adcp, Ea, 'k')
hold on
plot([time_adcp(1) time_adcp(end)], [thresh thresh], 'r')
hold off

Ea_raw = [Ea; time_adcp];
save Ea_raw.mat Ea_raw

save S2_Ea.mat Ea S2 zS2 time_adcp



% find peaks
count = 1;
for i = 2:length(Ea)-1
    if Ea(i) > Ea(i-1) & Ea(i) > Ea(i+1)
        if Ea(i)>=thresh
            IW(:,count) = [Ea(i); time_adcp(i)];
            count = count+1;
        end
    end
end

disp(sprintf('%d peaks detected'))
save IW_detected.mat IW

% $$$ % isolate IW
% $$$ tmin = datenum(2011, 09, 30, 0, 0, 0);
% $$$ tmax = datenum(2011, 09, 30, 23, 0, 0);
% $$$ tmin = datenum(2011, 10, 1, 0, 0, 0);
% $$$ tmax = datenum(2011, 10, 1, 23, 0, 0);
% $$$ I = find(time_adcp>tmin & time_adcp<tmax);
% $$$ x = time_adcp(I);
% $$$ y = Ea(I);
% $$$ plot(x, y)
% $$$ 
% $$$ set(gca, 'ygrid', 'on')
% $$$ xlabel('time') 
% $$$ ylabel('Ea (J/m^2)')
% $$$ datetick('x')


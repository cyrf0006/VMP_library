function eps_currents(vel_file, time_file, eps_files, dt, zrange)
    
% ex:
% eps_currents('../sept_2011_mission/Mouillages/M_N080/M_N080_vel.mat','../sept_2011_mission/Mouillages/M_N080/M_N080_PTzt.mat','file_names_bndry', 0.5/24, [20 100])

        
% --------------- ADCP velocities ----------------- %
load(vel_file);
load(time_file);

% find bottom ()
for i = 1:size(w, 1)
    I = find(isnan(w(i, :))==1);
    if length(I)>size(w, 2)./4
        bot = i;
        break
    end
end

% reduce file to desired depth range
u = east_vel(1:bot-1, :);
v = north_vel(1:bot-1, :);
w = w(1:bot-1, :);
z = z(1:bot-1);
dz = z(2)-z(1);

% reduce time range
tmin = datenum(2011, 09, 20, 0, 0, 0);
tmax = datenum(2011, 10, 12, 13, 0, 0);
I = find(time_adcp>tmin & time_adcp<tmax);
time_adcp = time_adcp(I);
w = w(:,I);
u = u(:,I);
v = v(:,I);

% time reduction
t1 = round(time_adcp(1)*24)/24;
t2 = round(time_adcp(end)*24)/24;
t_bin = t1:dt:t2;


for i = 1:length(t_bin);
    I = find(time_adcp >= t_bin(i)-dt/2 & time_adcp < t_bin(i)+dt/2);
    u_bin(:,i) = nanmean(u(:,I), 2);
    v_bin(:,i) = nanmean(v(:,I), 2);
    w_bin(:,i) = nanmean(w(:,I), 2);        
end

% vertical average of velocity profiles
I = find(z >= zrange(1) & z <= zrange(2));
mean_u = nanmean(u_bin, 1);
mean_v = nanmean(v_bin, 1);
mean_w = nanmean(w_bin, 1);
% ------------------------------------------------------ %


% ----------------- Deal with epsilon ------------------ %
fid = fopen(eps_files);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});
no_profile = size(epsfiles, 1);

for profile = 1:no_profile

    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
    maxp(profile) = max(p_eps1(end), p_eps2(end));
 
    % ---- average of 2 epsilon profiles ---- %
    if length(p_eps1) == length(p_eps2) % no problem, easy
        p_k = p_eps1;
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        EPS2(:,2) = eps2;
    elseif min(p_eps1)<min(p_eps2) % p1 drive, nomatter size of p2!
        p_k = p_eps1;        
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        ind = find(p_eps2 == p_eps1(1));
        EPS2(ind:length(p_eps2),2) = eps2;
    else % p2 drive, nomatter size of p1!
        p_k = p_eps2;
        p_N = p_eps1;
        EPS2=nan(length(p_eps2),2);
        EPS2(:,2) = eps2;
        ind = find(p_eps1 == p_eps2(1));
        EPS2(ind:length(p_eps1),1) = eps1;
    end
    
    % "selection" average
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
   
    EPS = nanmean(EPS2,2); %MEAN EPSILON
                        
    
    % if only nans
    if sum(~isnan(EPS))==0;
        continue
    end
    
% $$$     % Homestyle despike
% $$$     [Y, No] = Dspike(EPS, 5, 8);
% $$$ 
% $$$     % uses home despike
% $$$     EPS = Y;
% $$$     
    I = find(p_k >= zrange(1) & p_k <= zrange(2));
    epsilon(profile) = nanmean(EPS(I));
    time_eps(profile) = mtime_eps(1);
    
end 
% ------------------------------------------------------------ %



% ----------- find corresponding eps and u,v,u --------------- %
u3 = ((abs(mean_u) + abs(mean_v))/2).^3;
w2 = w.^2;

count = 1;
clear EPS
for i = 1:length(epsilon)
    [Y, I] = min(abs(time_eps(i)-t_bin));  
    if Y <= dt
        EPS(count) = epsilon(i);
        U3(count) = u3(i);
        W2(count) = w2(i);
        count = count+1;
    end
end


figure(1)
loglog(U3, EPS, '.k')
xlabel('U^3');
ylabel('\epsilon(W/kg)')


figure(2)
loglog(W2, EPS, '.k')
xlabel('U^3');
ylabel('\epsilon(W/kg)')



figure(3)
plot(EPS)
hold on
plot(log10(1e-3*U3/20), 'r')
hold off
ylabel('\epsilon')

keyboard
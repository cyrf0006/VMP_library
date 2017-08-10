function ms_transect(adcp_files, simrad_file, figurefile)
    
    
% function ms_transect(adcp_files)
%
% ex: ms_transect('adcp_24sept_fic1', 'transect2_38_Sv.mat', 'adcp_simrad_t2.png') 
% ex: ms_transect('adcp_24sept_fic1', 'transect1_38_Sv.mat', 'adcp_simrad_t1.png') 
% ex: ms_transect('adcp_24sept_fic1', 'transect3_38_Sv.mat', 'adcp_simrad_t3.png') 

    
% default parameters
zmax_adcp = 230; % adcp maximum depth
z0 = 1; %m default ADCP depth if no pressure sensor    
mag_var = 0;
theta = 42.5; % estuary angle (from Coriolis heading)


load(simrad_file)

% load *.P files names (file in which are recorded .P files)
fid = fopen(adcp_files);
C = textscan(fid, '%s', 'delimiter', '\n');

files = char(C{1});

no_files = size(files, 1); %number of eps_files 

count = 1;
disp('merge ADCP files...')
for i = 1:no_files

    load(files(i, :));
    
    % adjust year
    AA = datestr(ADCP.mtime, 0);
    for j = 1:size(AA, 1)
        AA(j,8:11) = '2011';
    end
    ADCP.mtime = datenum(AA);
    
    % -- merge files -- %
    if count == 1 
        east_vel  = ADCP.east_vel;
        north_vel  = ADCP.north_vel;
        w = ADCP.vert_vel;
        error_vel = ADCP.error_vel;
        time_adcp = ADCP.mtime;
        lon = ADCP.nav_slongitude;
        lat = ADCP.nav_slatitude;
        z_adcp = ADCP.config.ranges;
        bt_east = ADCP.bt_vel(1,:);
        bt_north = ADCP.bt_vel(2,:);
        bt_vert = ADCP.bt_vel(3,:);
        heading = ADCP.heading;
        count = length(time_adcp)+1;
    else
        east_vel(:,count:count+length(ADCP.mtime)-1) = ADCP.east_vel;
        north_vel(:,count:count+length(ADCP.mtime)-1) = ADCP.north_vel;
        w(:,count:count+length(ADCP.mtime)-1) = ADCP.vert_vel;
        error_vel(:,count:count+length(ADCP.mtime)-1) = ADCP.error_vel;
        time_adcp(count:count+length(ADCP.mtime)-1) = ADCP.mtime;
        lon(count:count+length(ADCP.mtime)-1) = ADCP.nav_slongitude;
        lat(count:count+length(ADCP.mtime)-1) = ADCP.nav_slatitude;
        bt_east(count:count+length(ADCP.mtime)-1) = ADCP.bt_vel(1,:);
        bt_north(count:count+length(ADCP.mtime)-1) = ADCP.bt_vel(2,:); 
        bt_vert(count:count+length(ADCP.mtime)-1) = ADCP.bt_vel(3,:);
        heading(count:count+length(ADCP.mtime)-1) = ADCP.heading;
        count = count+length(ADCP.mtime);
    end
end


% Correct velocities for bottom track
for i = 1:length(z_adcp)
    east_vel(i, :) = east_vel(i, :) - bt_east./1000;
    north_vel(i, :) = north_vel(i, :) - bt_north./1000;
    w(i, :) = w(i, :) - bt_vert./1000;
end

% select good time vector
time_sim = datenum(2011, 09, 24) + utc_mtime;
I = find(time_adcp>=time_sim(1) & time_adcp<=time_sim(end));
E = east_vel(:,I);
N = north_vel(:,I);
w = w(:,I);
er = error_vel(:,I);
time_adcp = time_adcp(I);
lon_adcp = lon(I);
lat_adcp = lat(I);
heading = heading(I);

% reduce vertical range
I = find(z_adcp >= z(1) & z_adcp <= z(end));
z_adcp = z_adcp(I);
E = E(I, :);
N = N(I, :);
w = w(I, :);
er = er(I, :);

% correction for unreal values
I = find(abs(E)>10);
E(I) = NaN;
N(I) = NaN;
w(I) = NaN;
er(I) = NaN;

keyboard

% appply mask on ADCP velocities
%(remove double echo, make sure no data under max depth)
for i = 1:length(time_adcp)
    [Y, I] = min(abs(time_adcp(i)-time_sim));
    J = find(z_adcp > H(I) | z_adcp > zmax_adcp);
    D(i) = H(I);
    E(J, i) = NaN;
    N(J, i) = NaN;
    w(J, i) = NaN;
end

% match ADCP on simrad grid (interp in both dimension)
E1 = nan(length(z), length(time_adcp));    % 1st itp (horiz)
N1 = nan(length(z), length(time_adcp));
w1 = nan(length(z), length(time_adcp));
EE = nan(size(Sv')); % final matrix
NN = nan(size(Sv'));
W = nan(size(Sv'));

% vertical itp
disp('vertical interpolation (may take some time)...')
for i = 1:length(time_adcp)
    Evec = E(:,i);
    Nvec = N(:,i);
    Wvec = w(:,i);
    % remove NaNs
    I1 = find(isnan(Evec)==1);
    I2 = find(isnan(Nvec)==1);
    I3 = find(isnan(Wvec)==1);
    if length(I1) > length(z_adcp)-3
        continue
    else  
        Evec(I1) = [];
        Nvec(I2) = [];
        Wvec(I3) = [];
        xvec1 = z_adcp;
        xvec2 = z_adcp;
        xvec3 = z_adcp;
        xvec1(I1)=[];
        xvec2(I2)=[];
        xvec3(I3)=[];
        E1(:, i) = interp1(xvec1, Evec, z);
        N1(:, i) = interp1(xvec2, Nvec, z);
        W1(:, i) = interp1(xvec3, Wvec, z);

    end
end

% correction for sidelobe effect
Rmax = D.*cosd(ADCP.config.beam_angle);
for i = 1:size(E1, 2);
    I = find(z > Rmax(i));
    E1(I, i) = NaN;
    N1(I, i) = NaN;
    W1(I, i) = NaN;
end
% $$$ figure(1)
% $$$ clf
% $$$ imagesc(time_adcp, z, E1)
% $$$ caxis([-.5 .5])
% $$$ hold on
% $$$ plot(time_adcp, D, 'r')
% $$$ plot(time_adcp, Rmax, 'm')
% $$$ hold off


disp('horizontal interpolation...')
I = find(isnan(E1)==1);
for i = 1:length(z)
    Evec = E1(i, :);
    Nvec = N1(i, :);
    Wvec = W1(i, :);

    % remove NaNs
    I1 = find(isnan(Evec)==1);
    I2 = find(isnan(Nvec)==1);
    I3 = find(isnan(Wvec)==1);
    if length(I1) > length(time_adcp)-3
        continue
    else    
        Evec(I1) = [];
        Nvec(I2) = [];
        Wvec(I3) = [];
        xvec1 = time_adcp;
        xvec2 = time_adcp;
        xvec3 = time_adcp;
        xvec1(I1)=[];
        xvec2(I2)=[];
        xvec3(I3)=[];
        % interp
        EE(i, :) = interp1(xvec1, Evec, time_sim);
        NN(i, :) = interp1(xvec2, Nvec, time_sim);
        W(i, :) = interp1(xvec3, Wvec, time_sim);
    end
end

% must mask (remove itp through topography)
Rmax = H.*cosd(ADCP.config.beam_angle);
for i = 1:size(EE, 2);
    I = find(z > Rmax(i));
    EE(I, i) = NaN;
    NN(I, i) = NaN;
    W(I, i) = NaN;
end


% ----------- Tranform time data to space data --------------- %
[x, I] = sort(x);
EE = EE(:,I);
NN = NN(:,I);
W = W(:,I);
time_sim = time_sim(I);
H = H(I);
Rmax = Rmax(I);

if exist('Sv') == 1
    Sv = Sv';
    S = Sv(:,I);
elseif exist('Sp') == 1
    Sp = Sp';
    S = S(:,I);
else
    disp('No Sv nor Sp variable?')
    keyboard
end

x_reg = min(x):mean(diff(x)):max(x);
E_xreg = nan(length(z), length(x_reg));
N_xreg = nan(length(z), length(x_reg));
W_xreg = nan(length(z), length(x_reg));
S_xreg = nan(length(z), length(x_reg));

% Itp
disp('time2space...')
for i = 1:size(EE, 1)
    Evec = EE(i, :);
    Nvec = NN(i, :);
    Wvec = W(i, :);
    Svec = S(i, :);

    I1 = find(~isnan(Evec)==1);
    I2 = find(~isnan(Nvec)==1);
    I3 = find(~isnan(Wvec)==1);
    I4 = find(~isnan(Svec)==1);
    if length(I1) > 1
        E_xreg(i, :) = interp1(x(I1), Evec(I1), x_reg);
        N_xreg(i, :) = interp1(x(I2), Nvec(I2), x_reg);
        W_xreg(i, :) = interp1(x(I3), Wvec(I3), x_reg);
    end
    if length(I4) > 1
        S_xreg(i, :) = interp1(x(I4), Svec(I4), x_reg);
    end
end
% ---------------------------------------------------------- %

% must mask again (remove itp through topography)
Rmax = interp1(x, Rmax, x_reg);
H_xreg = interp1(x, H, x_reg);
for i = 1:size(EE, 2);
    I = find(z > Rmax(i));
    E_xreg(I, i) = NaN;
    N_xreg(I, i) = NaN;
    W_xreg(I, i) = NaN;
    I = find(z > H_xreg(i));
    S_xreg(I, i) = NaN;
end


% correct for magnetic declination AND along-cross shore 
%[U, V] = rotate_vecd(EE, NN, (mag_var+theta));
[U, V] = rotate_vecd(E_xreg, N_xreg, (mag_var+theta));


% Smooth filtering data
% $$$ dt_sim = (abs(time_sim(2) - time_sim(1))).*86400;
% $$$ fs = 1/dt_sim;
% $$$ freq_low = 1/300; %Hz
% $$$ Wn_low = freq_low/(fs/2)

fx = 1/5; % about 5m sampling
dsmooth = 1000;%m
Wn_low = 1/dsmooth/(fx/2);

%[b,a] = butter(4,[Wn_low Wn_high]);
[b,a] = butter(4,Wn_low, 'low');

disp('smooth-filering data ...')
U_filt = nan(size(U));
V_filt = nan(size(V));
W_filt = nan(size(W));

% Crazy loop to smooth non-NaN segments
for i = 1:size(U, 1)
    I = find(~isnan(U(i, :))==1);
    if isempty(I)==1
        continue
    elseif max(diff(I)) == 1 % one segment to itp
        U_filt(i, I) = filtfilt(b, a, U(i, I));
        V_filt(i, I) = filtfilt(b, a, V(i, I));
        W_filt(i, I) = filtfilt(b, a, W(i, I));
    else % 2 or more segments to itp
        J = find(diff(I)>1);
        J = [1 J length(I)];
        for j = 1:length(J)-1
            seg = I(J(j)+1):I(J(j+1));
            
            if length(seg) < 20
                continue
            else
                U_filt(i, seg) = filtfilt(b, a, U(i, seg));
                V_filt(i, seg) = filtfilt(b, a, V(i, seg));
                W_filt(i, seg) = filtfilt(b, a, W(i, seg));
            end
        end
    end
end
disp('done!')


% ------------------------- Plot section ------------------------ %

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 12])   


subplot(311)
imagesc(x, z, U_filt)
caxis([-.5 .5])
ylim([0 350])
set(gca, 'xticklabel', [])
hold on
H2 = ones(size(H))*350;
patch([x fliplr(x) x(1)], [H fliplr(H2) H(1)], [.1 .1 .1], 'edgecolor', 'none');
hold off
adjust_space
text(7.5, 315, 'U (m/s)', 'color', [1 1 1]); 

subplot(312)
imagesc(x, z, V_filt)
caxis([-.5 .5])
ylim([0 350])
set(gca, 'xticklabel', [])
ylabel('depth (m)')
hold on
H2 = ones(size(H))*350;
patch([x fliplr(x) x(1)], [H fliplr(H2) H(1)], [.1 .1 .1], 'edgecolor', 'none');
hold off
adjust_space
c = colorbar('position', [.92 .52 .03 .3]);
text(7.5, 315, 'V (m/s)', 'color', [1 1 1]); 


subplot(313)
imagesc(x, z, S_xreg)
caxis([-120 -20])
ylim([0 350])
xlabel('distance from southshore (km)')
adjust_space
c = colorbar('position', [.92 .12 .03 .21]);

hold on
H2 = ones(size(H))*350;
patch([x fliplr(x) x(1)], [H fliplr(H2) H(1)], [.1 .1 .1], 'edgecolor', 'none');
hold off

text(37, 280, datestr(time_adcp(1), 1), 'color', [1 1 1]);
text(37, 315, [datestr(time_adcp(1), 15) ' - ' datestr(time_adcp(end), 15)], 'color', [1 1 1]); 

if exist('Sv')==1
    text(7.5, 315, 'S_v', 'color', [1 1 1]); 
else
        text(7.5, 315, 'S_p', 'color', [1 1 1]); 
end

% ----------------------- Plot section (end) ----------------------- %



%print('-dpng', '-r300', 'adcp_simrad_t2.png')
print('-dpng', '-r300', figurefile)

keyboard

datafile = figurefile(1:end-4);
save datafile
% Could wrote a portion of script to save data here...
function timedepth_eps_shear(eps_files, epsMin, epsMax, varargin)

% function timedepth_eps(no_profile, epsMin, epsMax, varargin)
% was previously isopyc_tidecycle.m
% 
% ex: timedepth_eps_shear('epsProfiles1001.list', 1e-9, 1e-5, 180, 1)

% "ls -1 eps_profile*.mat | sed 's/\.mat//' > eps_files"
%
% 
%



tt0 = datenum(2009, 10, 01, 16, 45, 0);

% Varargin test
if isempty(varargin)==1
    zmax = 324.5;
    zbin = 1;
    zmin = 0.5;
elseif size(varargin,2)==1
    zmax = varargin{1};
    zbin = 1;
    zmin = zbin./2;
elseif size(varargin,2)==2
    zmax = varargin{1};
    zbin = varargin{2};
    zmin = zbin./2;
else
    disp('Wrong input... try "help var_profile_cal"')
    return
end


% ------- for water column stability -------- %
ddz = .25;
P_N2 = ddz/2:ddz:max(zmax);
g = 9.81;
% ------------------------------------------- %


% load eps_files file
fid = fopen(eps_files);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

no_profile = size(epsfiles, 1); %number of *.P files 

% find time min
eps_fname = epsfiles(1, :);
I = find(eps_fname==' '); %remove white space  
eps_fname(I) = [];
data_fname = eps_fname(5:end);
load(data_fname);

T1 = [str2num(time(1, 1:2)) str2num(time(1, 4:5)) str2num(time(1, 7:8))]; % [hh mm ss] (num vector)

% find time max
eps_fname = epsfiles(end, :); 
I = find(eps_fname==' '); %remove white space  
eps_fname(I) = [];
data_fname = eps_fname(5:end);
load(data_fname);
sizeT = size(time);
T2 = [str2num(time(sizeT(1), 1:2)) str2num(time(sizeT(1), 4:5)) str2num(time(sizeT(1), 7:8))]; % [hh mm ss] (num vector
%set(gca, 'FontSize', fs)
%colorbar
dat = date(1,:);
dd = dat(1:2); 
mm = dat(4:5); 
yyyy = dat(7:10);


%zbin=1;
zvec = [zmin:zbin:zmax]';
rho_mat = nan(size(zvec,1), no_profile);
eps_mat = nan(size(zvec,1), no_profile); 
T_mat = nan(size(zvec,1), no_profile);
K_mat = nan(size(zvec,1), no_profile);
F_mat = nan(size(zvec,1), no_profile);
N_mat = nan(size(zvec,1), no_profile);
N2raw_mat = nan(size(P_N2,1), no_profile);

%%%%%%%%%%%%%%%%%%%%
%  discretization  %
%%%%%%%%%%%%%%%%%%%%

% time discretisation
dT = 0.001; %this can be adapted to modify plotting
Tmin = fix((T1(1)+T1(2)/60+T1(3)/60/60)/dT)*dT; %time in decimal round to 1/dT min.
Tmax = ceil((T2(1)+T2(2)/60+T2(3)/60/60)/dT)*dT; 
time_vector = Tmin:dT:Tmax;


%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%

bin=1;

fid = fopen('./eps_plot_info.mat');
if fid == -1 % doesnt exist yet!
% ---  LOOP on profiles ---- %
for profile = 1:no_profile

    eps_fname = epsfiles(profile, :);
    I = find(eps_fname==' '); %remove white space  
    eps_fname(I) = [];
    data_fname = eps_fname(5:end);
    disp(data_fname);
    
    load(data_fname);
    load(eps_fname);
    
    % time of survey
    if profile==1
        t0 = mtime(1);
    elseif profile==no_profile
        tf = mtime(end);
    end
    
    
    % Test for Fred (pas sur que dans tous les fichiers p_eps1 et p_eps2 sont les memes...)
    if length(p_eps1)~=length(p_eps2) | ( p_eps1(1)~=p_eps2(1) & ...
                                          isnan(p_eps2(1))==0 )
        disp('p_eps1 and p_eps2 mismatch')
        pause
    end
        
    % Combination of the 2 profiles
    p_eps = nanmean([p_eps1; p_eps2]);
    eps = nanmean([eps1; eps2]);  
    
    % Diffusivity
    K = 0.2.*eps./(N.^2);
    % density 
    rho = sw_dens(SBS, SBT, P);
    % Raw buoy. freq.
    for i = 1:length(P_N2)
        I = find(P >= (P_N2(i) - ddz/2) & P < (P_N2(i) + ddz/2)); 
        if length(I) > 5
            pol = polyfit(P(I),rho(I),1);
            rho0 = mean(rho(I));
            N2raw_mat(i, profile) = (g/rho0)*pol(1);
        end
    end

    % bin density to epsilon
    I = find(isnan(rho)==1);
    rho(I)=[];
    PP = P;
    PP(I) = [];    
    rhoI = interp1(PP, rho, p_eps);
       
    % bin temperature to epsilon
    I = find(isnan(SBT)==1);
    temp = SBT;
    temp(I)=[];
    PP = P;
    PP(I) = [];    
    tempI = interp1(PP, temp, p_eps);

    % identify bottom
    bot(profile) = max(P);
     
    % time vector
    timevec(profile) = mtime(1);
    
    for i = 1:length(p_eps);
        [Y I] = min(abs(zvec-p_eps(i)));        
        rho_mat(I, profile) = rhoI(i);
        eps_mat(I, profile) = eps(i);
        T_mat(I, profile) = tempI(i);
        K_mat(I, profile) = K(i);
        N_mat(I, profile) = N(i);
        %F_mat(I, profile) = fluoI(i);
        %TR_mat(I, profile) = transI(i);
    end
    
end 

% For 3 profiles we were far from bottom (2009-10-01)
bot(5:7) = bot(4);



% ---- END loop on profiles ---- % 

% clean density field 
% (linear itp in time to remove nans)
for i = 1:size(rho_mat,1)
    raw_vec = rho_mat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        rho_mat(i,I(1):I(end)) = vec_itp;
    end
end

% clean temperature field 
% (linear itp in time to remove nans)
for i = 1:size(T_mat,1)
    raw_vec = T_mat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        T_mat(i,I(1):I(end)) = vec_itp;
    end
end

% clean N2 field 
% (linear itp in time to remove nans)
for i = 1:size(N_mat,1)
    raw_vec = N_mat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        N_mat(i,I(1):I(end)) = vec_itp;
    end
end


% ---- For mixing efficiency calculation ---- %
timevec_eps = timevec;
zvec_eps = zvec;
save MixEff_eps.mat eps_mat N_mat timevec_eps zvec_eps
% ------------------------------------------- %

save eps_plot_info.mat timevec zvec eps_mat K_mat rho_mat bot

else
    load eps_plot_info.mat
end




%% ---- Getting info on shear ---- %
load('~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/ADCP/adcp_2009-10-01.mat');
gpsFile = '~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/GPS/20091001';


east_vel = ADCP.east_vel;
north_vel = ADCP.north_vel;
mtime = ADCP.mtime;
z = ADCP.config.ranges;

% Reduce to timevec
I = find(mtime >= timevec(1) & mtime <= timevec(end));
east_vel = east_vel(:,I);
north_vel = north_vel(:,I);
mtime = mtime(I);

I = find(~isnan(nansum(east_vel,1))==1);
east_vel = east_vel(:,I);
north_vel = north_vel(:,I);
mtime = mtime(I);

[u v] = adcp_velcor(east_vel, north_vel, mtime, gpsFile);


% Filtering
filterTime = 5; % in minutes 
dt = (mtime(2)-mtime(1));
% Low pass 
fs = 1/(dt*86400);
freq_low = 1/(filterTime*60); %Hz
Wn_low = freq_low/(fs/2);
[b,a] = butter(4, Wn_low);


u_filt = nan(size(u));
v_filt = nan(size(v));

for i = 1:size(u, 1)
    % remove NANs
    uvec = u(i,:);
    I = find(isnan(uvec)==1);
    xitp = 1:length(uvec);
    x = xitp;
    uvec(I) = [];
    x(I) = [];
    uvec = interp1(x, uvec, xitp);    
    u_filt(i, :) = filtfilt(b, a, uvec);
    
    vvec = v(i,:);
    I = find(isnan(vvec)==1);
    xitp = 1:length(vvec);
    x = xitp;
    vvec(I) = [];
    x(I) = [];
    vvec = interp1(x, vvec, xitp);    
    v_filt(i, :) = filtfilt(b, a, vvec);
end


du = diff(u_filt,1,1);
dv = diff(v_filt,1,1);
dz = z(2) - z(1);

S2Mat = (du/dz).^2 + (dv/dz).^2;
zS2 = z(1:end-1)+dz/2;


%save ADCP_shear.mat S2Mat zS2 mtime


% $$$ % average shear
% $$$ dtShear = 60/86400;
% $$$ timeVecShear = mtime(1):dtShear:mtime(end);
% $$$ S2_bin = nan(size(S2Mat,1), length(timeVec))
% $$$ 
% $$$ for i = 1:length(timeVecShear)
% $$$     I = find(mtime>=timeVecShear(i)-dtShear/2 & mtime< ...
% $$$              timeVecShear(i)+dtShear/2);
% $$$     S2_bin(i) = 

% $$$ I = find(mtime > datenum(2012, 10, 25, 14, 48, 0) & mtime < datenum(2012, 10, 25, 16, 0, 0) );
% $$$ J = find(zS2<80);
% $$$ 
% $$$ % Must now find the bottom
% $$$ Hobs = [79 79 65]+1; % +1m ADCP depth
% $$$ Xobs = [datenum(2012, 10, 25, 14, 48, 0) datenum(2012, 10, 25, 14, 50, 0) datenum(2012, 10, 25, 16, 0, 0)];


H = interp1(timevec, bot, mtime);
%S2 = S2Mat(:,I);
%S2 = S2Mat;

for i = 1:length(mtime)
    II = find(zS2>H(i));
    S2(II, i) = NaN;
end




%% ---- PLot Section ---- %%
h = figure('Visible', 'off');
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 25 15])

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.05 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.2; % very left of figure
rigs = 0.2; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

z0 = 1.6; zf = 38.6;
t0 = datenum(2009, 10, 01, 17, 30, 0);
tf = datenum(2009, 10, 01, 19, 15, 0);
xRect = [t0 tf tf t0 t0];
yRect = [zf zf z0 z0 zf];


%% -- S1 -- %%
subplot(211)


%imagesc(mtime(I), zS2(J), log10(S2))
%contourf(mtime(I), zS2(J), log10(S2(J,I)))
I = find(zS2<=30); % remove depth >30m (sidelobe)
contourf(mtime, zS2(I), log10(S2Mat(I,:)), 'linestyle', 'none')
caxis([-6 0])
%imagesc(mtime(I), zS2(J), log10(S2(J,:)))
V = 1020:0.5:1028;

hold on
% bottom
% $$$ % from profile end
% $$$ x = [timevec(1) timevec  timevec(end) timevec(1)];
% $$$ y = [zmax bot zmax zmax];
% from ADCP
topo = load('~/WINDEX/data_processing/Mission_tadoussac_2009/HLC_topo_2009-10-01.mat');
x = [topo.timeVec(1) topo.timeVec  topo.timeVec(end) topo.timeVec(1)]; 
y = [zmax topo.bot zmax zmax];


contour(timevec, zvec, rho_mat, V, 'color', [1 1 1]*0) 
patch(x, y, [1 1 1]*.6)
plot(xRect, yRect, '--k', 'linewidth', 2)
hold off

%caxis([-5 -2])
set(gca, 'ydir', 'reverse')
datetick('x', 15) 
xlim([tt0 timevec(end)])
ylim([0 35])
cb = colorbar;
%xlabel(datestr(timevec(1), 1))
ylabel('Depth (m)')

ylabel(cb,'log_{10}(S^2 (s^{-2}))', 'FontSize', 10)
set(gca, 'xticklabel', [])


adjust_space


cbPos = get(cb, 'pos');
%cbPos(1) = cbPos(1)-.01;
cbPos(2) = cbPos(2)+.025;
cbPos(3) = cbPos(3)*.5;
cbPos(4) = cbPos(4)-.05;
set(cb, 'pos', cbPos)



%% -- S2 -- %%
subplot(212)


% -- density with dissipation -- %
epsMin = 1e-9;
epsMax = 1e-5;

contour(timevec, zvec, rho_mat, V, 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
hold on
field=eps_mat;
x=timevec;
y=zvec;
for j = 1:size(field,2) 
    
    %data = eps_mat(:,j);
    data = field(:,j);
    x(1:length(data)) = timevec(j);
    y = zvec;
      
    rectbar = [0.92 0.15 0.01 0.75]; 

    if j == size(field,2) % with cbar
        cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
    end
 
end
contour(timevec, zvec, rho_mat, V, 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
%x = [timevec(1) timevec  timevec(end) timevec(1)];
%y = [zmax bot zmax zmax];
x = [topo.timeVec(1) topo.timeVec  topo.timeVec(end) topo.timeVec(1)]; 
y = [zmax topo.bot zmax zmax];
patch(x, y, [1 1 1]*.6)
plot(xRect, yRect, '--k', 'linewidth', 2)
hold off
datetick('x', 15)

xlim([tt0 timevec(end)])
ylim([0 35])
colorbar(cb, 'delete');
cb = colorbar;
ylabel(cb,'log_{10}(\epsilon (W kg^{-1}))', 'FontSize', 10)
set(gca, 'ydir', 'reverse')
%set(gca, 'xticklabel', [])
xlabel(datestr(timevec(1), 1))

adjust_space

cbPos = get(cb, 'pos');
%cbPos(1) = cbPos(1)-.01;
cbPos(2) = cbPos(2)+.025;
cbPos(3) = cbPos(3)*.5;
cbPos(4) = cbPos(4)-.05;
set(cb, 'pos', cbPos);


print(h, '-dpng', 'eps_shear.png')
set(h, 'renderer', 'painters')
print(h, '-depsc2', 'eps_shear.eps')


keyboard


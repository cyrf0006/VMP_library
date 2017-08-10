function timedepth_eps_shear2(eps_files, epsMin, epsMax, varargin)

% function timedepth_eps(no_profile, epsMin, epsMax, varargin)
% was previously isopyc_tidecycle.m
% 
% ex: timedepth_eps_shear2('epsProfiles1001.list', 1e-9, 1e-5, 180, 1)
%
% Now run from anywhere:
% ex: timedepth_eps_shear2('~/PhD/WINDEX/MissionTadoussac/epsProfiles1001.list', 1e-9, 1e-5, 35, 1)



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
data_fname = strrep(eps_fname, 'eps_', '');
load(data_fname);

T1 = [str2num(time(1, 1:2)) str2num(time(1, 4:5)) str2num(time(1, 7:8))]; % [hh mm ss] (num vector)

% find time max
eps_fname = epsfiles(end, :); 
I = find(eps_fname==' '); %remove white space  
eps_fname(I) = [];
data_fname = strrep(eps_fname, 'eps_', '');
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
NO3_mat = nan(size(zvec,1), no_profile);
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
    data_fname = strrep(eps_fname, 'eps_', '');
    disp(data_fname);
    
    load(data_fname);
    load(eps_fname);
    
    % time of survey
    if profile==1
        t0 = mtime(1);
    elseif profile==no_profile
        tf = mtime(end);
    end
    
     % ---- Process SA and NO3
    if ~exist('SA')
        [SA, in_ocean] = gsw_SA_from_SP(SBS,P,-68.79,48.61);
        CT = gsw_CT_from_t(SA,SBT,P);
    end

% $$$     rho = gsw_rho(SA,CT,P);        
% $$$     [Y, I] = sort(rho);
% $$$     SA = SA(I);
% $$$     %NO3 = 6.7*SA - 207; %60-250m
% $$$     NO3 = 5.7*SA - 173; %50-300m
% $$$     
% $$$     I = find(NO3<=0);
% $$$     NO3(I) = NaN;
% $$$     clear SA;
    
    SA = sort(SA);
    NO3 = 7.3*SA - 227; %32-34.3 SA
    I = find(NO3<=0);
    I = find(SA<32 | SA>34.3);
    NO3(I) = NaN;
    clear SA;
    
    for i = 1:length(zvec) %bin NO3
        I = find(P >= zvec(i)-zbin/2 & P <= zvec(i)+zbin/2);
        NO3_mat(i, profile) = nanmean(NO3(I));
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
    end
end 

% For 3 profiles we were far from bottom (2009-10-01)
bot(5:7) = bot(4);
% ---- END loop on profiles ---- % 


% ---- NO3 fluxes calculation ---- %
dNO3dz = diff(NO3_mat, 1, 1)./zbin;
Pdiff = zvec(1:end-1)+zbin/2;
FNO3 = nan(size(K_mat));

for i = 1:size(dNO3dz,2)
    dNO3itp = interp1(Pdiff, dNO3dz(:,i), zvec);
    FNO3(:,i) = K_mat(:,i).*dNO3itp;
end


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

save eps_plot_info.mat timevec zvec eps_mat K_mat rho_mat FNO3 bot

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
%h = figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 25 30])

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 5; % no. subplot row
dx = 0.05 ; % horiz. space between subplots
dy = 0.035; % vert. space between subplots
lefs = 0.2; % very left of figure
rigs = 0.2; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
xtra_offset = 0.05;
cbar_offset = 0.01;
cbar_width = .02;

z0 = 1.6; zf = 38.6;
t0 = datenum(2009, 10, 01, 17, 30, 0);
tf = datenum(2009, 10, 01, 19, 15, 0);
xRect = [t0 tf tf t0 t0];
yRect = [zf zf z0 z0 zf];


%% -- S1 -- %%
subplot(5,1,[1 2])


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
%plot(xRect, yRect, '--k', 'linewidth', 2)
hold off

%caxis([-5 -2])
set(gca, 'ydir', 'reverse')
datetick('x', 15) 
xlim([tt0 timevec(end)])
ylim([0 35])
cb = colorbar;
%xlabel(datestr(timevec(1), 1))
ylabel('Depth (m)')

ti = ylabel(cb,'log_{10}(S^2 / s^{-2})', 'FontSize', 10);
set(gca, 'xticklabel', [])
tiPos = get(ti, 'pos');
tiPos(1) = tiPos(1)*1.5;
set(ti, 'pos', tiPos)

adjust_space
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = pos2(4)*2+xtra_offset;
set(gca, 'pos', pos2)
drawnow

cbPos = get(cb, 'pos');
cbPos(1) = pos2(1)+pos2(3)+cbar_offset;
cbPos(4) = pos2(4)*.8;
cbPos(2) = pos2(2)+(pos2(4)-cbPos(4))/2;
cbPos(3) = cbar_width;
set(cb, 'pos', cbPos)



%% -- S2 -- %%
subplot(5,1, [3 4])


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
ti = ylabel(cb,'log_{10}(\epsilon / W kg^{-1})', 'FontSize', 10);
set(gca, 'ydir', 'reverse')
tiPos = get(ti, 'pos');
tiPos(1) = tiPos(1)*2;
set(ti, 'pos', tiPos)
ylabel('Depth (m)')

adjust_space
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = pos2(4)*2+xtra_offset;
set(gca, 'pos', pos2)

cbPos = get(cb, 'pos');
cbPos(1) = pos2(1)+pos2(3)+cbar_offset;
cbPos(4) = pos2(4)*.8;
cbPos(2) = pos2(2)+(pos2(4)-cbPos(4))/2;
cbPos(3) = cbar_width;
set(cb, 'pos', cbPos)



%% -- S3 -- %%
subplot(5,1, 5)


load ~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/ADCP/adcp_2009-10-01.mat;

xticks = t0:15/60/24:tf;

intens = squeeze(ADCP.intens(:,4,:));
zVec = ADCP.config.ranges;
timeVec = ADCP.mtime;
zmax = max(zVec);

I = find(timeVec>=t0 & timeVec <=tf);
timeVec = timeVec(I);
intens = intens(:,I);

% find bottom
I = find(zVec>15);
[Y, J] = max(intens(I,:));
J = J+(size(intens,1)-length(I));
I = find(Y<120);
J(I) = length(zVec);
bot = nan(1, length(J));
for i = 1:length(J)
    bot(i) = zVec(J(i));
end
bot = runmean(bot, 3);

% echo attenuation correction
for i = 1:length(timeVec)
    intens(:,i) = intens(:,i).*20.*log10(zVec);
end


I = find(timevec>=t0 & timevec <=tf);
eps_mat = eps_mat(:,I);
K_mat = K_mat(:,I);
FNO3 = FNO3(:,I);
rho_mat = rho_mat(:,I);
timevec = timevec(I);


imagesc(timeVec, zVec, intens)
colormap(flipud(bone))

%caxis([50 150])
caxis([1000 3000])
freezeColors
drawnow
hold on


% dissipation
field=FNO3*86400; %in d-1
I = find(FNO3<0);
field(I) = NaN;
x=timevec;
y=zvec;
epsMin = 1e-5; % for K
epsMax = 1e-1;
epsMin = 1e-3; % for F
epsMax = 1e3;
rectbar = cbPos;

for j = 1:2:size(field,2)     
    data = field(:,j);
    x(1:length(data)) = timevec(j);
    y = zvec;      
    if j == 1 % with cbar
        cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
    end 
end


% isopycnals
%contour(timevec, zvec, rho_mat, 15, 'color', [1 1 1], 'linewidth', 1)
set(gca, 'xtick', xticks)
xlim([t0 tf])
%datetick('x', 15, 'keeplimits')
xlabel(datestr(timevec(1), 1))
set(gca, 'xticklabel', datestr(xticks, 15))
ti = ylabel(cb, 'log_{10}(F_{NO_3} / mmol m^{-2} d^{-1})');
%ti = ylabel(cb, 'log_{10}(K(m^2 s^{-1}))');
tiPos = get(ti, 'pos');
tiPos(1) = tiPos(1)*1.05;
set(ti, 'pos', tiPos)

% bottom
x = [timeVec(1) timeVec  timeVec(end) timeVec(1)];
y = [zmax bot zmax zmax];
patch(x, y, [1 1 1]*.6)
hold off

adjust_space

pos2 = get(gca, 'pos');
cbPos = get(cb, 'pos');
cbPos(1) = pos2(1)+pos2(3)+cbar_offset;
cbPos(4) = pos2(4)*.8;
cbPos(2) = pos2(2)+(pos2(4)-cbPos(4))/2;
cbPos(3) = cbar_width;
set(cb, 'pos', cbPos)



% *** A tricky procedure to add expansion lines out of axes... ***
hAx1 = gca;
XLIM = get(hAx1, 'xlim');
YLIM = [-20 35];

%# create a second axis as copy of first (without its content), 
%# reduce its size, and set limits accordingly
%hAx2 = copyobj(hAx1,gcf);
hAx2 = axes; %fake axes
             % few axes properties (note enlargement of the
             % vertical position by *1.5)
set(hAx2, 'Position',get(hAx1,'Position').*[1 1 1 1.5], ...
    'XLimMode','manual', 'YLimMode','manual', 'color', 'none', 'yticklabel', ' ')
 
plot([t0 datenum(2009, 10, 01, 17, 50,  0)],[-3 -9], 'k', 'linewidth', 2)
hold on
plot([datenum(2009, 10, 01, 18, 40,  0) tf],[-9 -3], 'k', 'linewidth', 2)
hold off
set(gca, 'ydir', 'reverse')  
set(hAx2, 'xlim', XLIM)  
set(hAx2, 'ylim', YLIM)

axis(hAx2,'off')
uistack(hAx1,'top')


print(h, '-dpng', 'eps_shear_KH.png')
set(h, 'renderer', 'painters')
print(h, '-depsc2', 'eps_shear_KH.eps')

keyboard


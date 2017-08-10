function timedepth_map_defense(eps_files, epsMin, epsMax, varargin)

% function timedepth_eps(no_profile, epsMin, epsMax, varargin)
% was previously isopyc_tidecycle.m
% 
% ex: timedepth_map_defense('epsProfiles1001.list', 1e-9, 1e-5, 200, 1)

% "ls -1 eps_profile*.mat | sed 's/\.mat//' > eps_files"
%
% 
%

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




fid = fopen('plotinfo_timedepth_map.mat');
if fid == -1;
    bin=1;
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
    bot(5:7) = bot(8);

    % BAthymetry from ADCP
    topo = load('~/WINDEX/data_processing/Mission_tadoussac_2009/HLC_topo_2009-10-01.mat');
    I = find(timevec>=topo.timeVec(1) & timevec<=topo.timeVec(end));
    I1 = I(1)-1;
    I2 = I(end)+1;
    timeVecTopo = [timevec(1:I1) topo.timeVec timevec(I2:end)];
    botVecTopo = [bot(1:I1)  topo.bot bot(I2:end)];


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

    save plotinfo_timedepth_map.mat

else
    disp(['** Use already built plot info (plotinfo_timedepth_map.mat) **'])
    load plotinfo_timedepth_map
end




%% VERSION 1 %%
% (no temp.) %
disp('  -> plotting ...')
figure('visible', 'off')
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.09; % very right of figure
tops = 0.01; % top of figure
bots = 0.09; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

FONTSIZE = 18;
FONTSIZE2 = 15;

% -- density contourplot -- %
%V = 1010:0.2:1028;
%cb = contour_rho(timevec, zvec, rho_mat, V);
contour(timevec, zvec, rho_mat, [1020:.3:1026.6], 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
hold on % add turbulence
field = eps_mat;
for j = 1:size(field,2) 
    data = field(:,j);
    x(1:length(data)) = timevec(j);
    y = zvec;
      
    rectbar = [0.92 0.15 0.01 0.8]; 

    if j == size(field,2) % with cbar
        cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
    end 
end

contour(timevec, zvec, rho_mat, [1023:.2:1026.6], 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
x = [timeVecTopo(1) timeVecTopo  timeVecTopo(end) timeVecTopo(1)];
y = [zmax botVecTopo zmax zmax];
patch(x, y, [1 1 1]*.6)
hold off

ylabel(cb,'\epsilon (W kg^{-1})', 'FontSize', FONTSIZE)
xlim([t0 tf])
set(gca, 'ydir', 'reverse')
set(gca, 'xtick', datenum(2009, 10, 1, 14:20,0,0))
set(gca, 'xticklabel', [])
ylabel('Depth (m)', 'fontsize', FONTSIZE)
text(datenum(2009, 10, 01, 20, 15, 0), 145, 'a', 'fontweight', 'bold', 'fontsize', 16)
datetick('x',15)
xlabel(datestr(timevec(1), 1), 'Fontsize', FONTSIZE)
xlim([t0 tf])

% Map
load '/home/cyrf0006/data/SHC/500m/CHS_500m_gulf.mat'

lon_min=-69.75;
lon_max=-69.25;
lat_min=48.05;
lat_max=48.33;

%Points for transect (x,y = lon,lat)
%A = [-68-50/60 48+52/60];
%B = [-68-25/60 48+32/60];
A = [-68.422966 48.537659];
B = [-68.837865 48.850145];

I=find(lat<lat_max & lat>lat_min);
latitude=lat(I);
longitude=lon(I);
bathy=z(I);

I=find(longitude<lon_max & longitude>lon_min);
latitude2=latitude(I);
longitude2=longitude(I);
bathy2=bathy(I);

clear latitude longitude I bathy lat lon z

% on the new grid
y = lat_min:(lat_max-lat_min)/110/5:lat_max;
x = lon_min:(lon_max-lon_min)/80/5:lon_max;
[X,Y] = meshgrid(x,y);
[XI,YI,Z] = griddata(longitude2,latitude2,bathy2,X,Y);

%b=gaussfir(0.001, 5, 3);  %Gaussian filter designer ()
b=gaussfir(10);  %Gaussian filter designer ()

Z2=filter2(b,Z);  %Application of the filter

set(gca, 'fontsize', FONTSIZE)
adjust_space

if zmax == 100
    a2 = axes('position',[0.5 0.15 0.4 0.4]) ; % inset
elseif zmax == 200
    a2 = axes('position',[0.5 0.15 0.5 0.5]) ; % inset
elseif zmax == 150
    a2 = axes('position',[0.45 0.15 0.53 0.53]) ; % inset
else
    a2 = axes('position',[0.5 0.15 0.4 0.4]) ; % inset
end         
pause(.5)    

m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
V=[0:10:300];
caxis([V(1) V(end)]) ;

%pack
hold on
%[H, H] = m_contourf(X,Y,Z2,V, 'linestyle', 'none');
[HH, HH] = m_contour(X,Y,Z2, [20 120 200 300], 'color', 'k');

m_gshhs_h('patch',[.7 .7 .7]); %coastlines                               
%m_coast('patch',[.5 .5 .5]); %coastlines

% $$$ m_line(-69.716606, 48.156604,'marker','o','MarkerFaceColor','k','markersize',6,'color','k');
% $$$ m_text(-69.8, 48.156604, 'Tadoussac', 'vertical', 'bottom', 'horizontal', 'center', 'color', 'k','FontSize',10, 'FontWeight', 'bold');

set(gca, 'fontsize', FONTSIZE2)
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])



% ---- Add few features ---- %

% Stat. 25
lon_min=-69.5846;
lon_max=-69.3145;
lat_min=48.1431;
lat_max=48.3229;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'k') 
m_text(lon_max, lat_min, 'Stat. 25 ', 'vertical', 'bottom', 'horizontal','right', 'color', [0 0 0],'FontSize',FONTSIZE2, 'FontWeight', 'bold')
	
% Add profile position
[latVec, lonVec] = whereare(eps_files, '~/WINDEX/data_processing/Mission_tadoussac_2009/mergedGPS');
for i = 1:length(latVec)
    m_line(lonVec(i), latVec(i),'marker','o','color',[1 0 0]*.9 ,'markersize',2);
end
m_text(lonVec(1), latVec(1), datestr(timevec(1),15), 'vertical', ...
       'bottom', 'horizontal','right', 'color', [1 0 0]*.9,'FontSize',FONTSIZE2, 'FontWeight', 'bold')
m_text(lonVec(end), latVec(end), datestr(timevec(end),15), 'vertical', ...
       'bottom', 'horizontal','right', 'color', [1 0 0]*.9,'FontSize',FONTSIZE2, 'FontWeight', 'bold')
m_text(lonVec(60), latVec(60), datestr(timevec(62),15), 'vertical', ...
       'top', 'horizontal','right', 'color', [1 0 0]*.9,'FontSize',FONTSIZE2, 'FontWeight', 'bold')


% Thermograph position
m_line(-69.6751, 48.1197, 'marker','p','color',[0 0 1]*.5 ...
       ,'markersize',12, 'markerFaceColor', [0 0 1]*.5);

m_grid('box','fancy', 'yticklabel', [], 'xticklabel', [])

set(gca, 'fontsize', FONTSIZE2)

print('-dpng', '-r300', 'timedepth_map1001_defense.png');
set(gcf, 'renderer', 'painters')
print('-depsc2', 'timedepth_map1001_defense.eps');


keyboard


%% VERSION 2 %%
% (with temp.) %
disp('  -> plotting ...')
figure('visible', 'off')
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 32])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.09; % very right of figure
tops = 0.01; % top of figure
bots = 0.08; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

FONTSIZE = 18;
FONTSIZE2 = 15;

subplot(211)
% -- density contourplot -- %
%V = 1010:0.2:1028;
%cb = contour_rho(timevec, zvec, rho_mat, V);
contour(timevec, zvec, rho_mat, [1020:.3:1026.6], 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
hold on % add turbulence
field = eps_mat;
for j = 1:size(field,2) 
    data = field(:,j);
    x(1:length(data)) = timevec(j);
    y = zvec;
      
    rectbar = [0.92 0.60 0.01 0.35]; 

    if j == size(field,2) % with cbar
        cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
    end 
end

contour(timevec, zvec, rho_mat, [1023:.2:1026.6], 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
x = [timeVecTopo(1) timeVecTopo  timeVecTopo(end) timeVecTopo(1)];
y = [zmax botVecTopo zmax zmax];
patch(x, y, [1 1 1]*.6)
hold off

ylabel(cb,'\epsilon (W kg^{-1})', 'FontSize', FONTSIZE)
xlim([t0 tf])
set(gca, 'ydir', 'reverse')
set(gca, 'xtick', datenum(2009, 10, 1, 14:20,0,0))
set(gca, 'xticklabel', [])
ylabel('Depth (m)', 'fontsize', FONTSIZE)
text(datenum(2009, 10, 01, 20, 15, 0), 145, 'a', 'fontweight', 'bold', 'fontsize', 16)

% Map
load '/home/cyrf0006/data/SHC/500m/CHS_500m_gulf.mat'

lon_min=-69.75;
lon_max=-69.25;
lat_min=48.05;
lat_max=48.33;

%Points for transect (x,y = lon,lat)
%A = [-68-50/60 48+52/60];
%B = [-68-25/60 48+32/60];
A = [-68.422966 48.537659];
B = [-68.837865 48.850145];

I=find(lat<lat_max & lat>lat_min);
latitude=lat(I);
longitude=lon(I);
bathy=z(I);

I=find(longitude<lon_max & longitude>lon_min);
latitude2=latitude(I);
longitude2=longitude(I);
bathy2=bathy(I);

clear latitude longitude I bathy lat lon z

% on the new grid
y = lat_min:(lat_max-lat_min)/110/5:lat_max;
x = lon_min:(lon_max-lon_min)/80/5:lon_max;
[X,Y] = meshgrid(x,y);
[XI,YI,Z] = griddata(longitude2,latitude2,bathy2,X,Y);

%b=gaussfir(0.001, 5, 3);  %Gaussian filter designer ()
b=gaussfir(10);  %Gaussian filter designer ()

Z2=filter2(b,Z);  %Application of the filter

set(gca, 'fontsize', FONTSIZE)
adjust_space


if zmax == 100
    a2 = axes('position',[0.5 0.15 0.4 0.4]) ; % inset
elseif zmax == 200
    a2 = axes('position',[0.5 0.15 0.5 0.5]) ; % inset
elseif zmax == 150
    % a2 = axes('position',[0.61 0.62 0.21 0.21]) ; % inset
    a2 = axes('position',[0.55 0.555 0.28 0.3]) ; % inset
else
    
    a2 = axes('position',[0.5 0.15 0.4 0.4]) ; % inset
end        
pause(.5)    

m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
V=[0:10:300];
caxis([V(1) V(end)]) ;

%pack
hold on
%[H, H] = m_contourf(X,Y,Z2,V, 'linestyle', 'none');
[HH, HH] = m_contour(X,Y,Z2, [20 120 200 300], 'color', 'k');

m_gshhs_h('patch',[.7 .7 .7]); %coastlines                               
%m_coast('patch',[.5 .5 .5]); %coastlines

% $$$ m_line(-69.716606, 48.156604,'marker','o','MarkerFaceColor','k','markersize',6,'color','k');
% $$$ m_text(-69.8, 48.156604, 'Tadoussac', 'vertical', 'bottom', 'horizontal', 'center', 'color', 'k','FontSize',10, 'FontWeight', 'bold');

set(gca, 'fontsize', FONTSIZE2)
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])



% ---- Add few features ---- %

% Stat. 25
lon_min=-69.5846;
lon_max=-69.3145;
lat_min=48.1431;
lat_max=48.3229;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'k') 
m_text(lon_max, lat_min, 'Stat. 25 ', 'vertical', 'bottom', 'horizontal','right', 'color', [0 0 0],'FontSize',FONTSIZE2, 'FontWeight', 'bold')
	
% Add profile position
[latVec, lonVec] = whereare(eps_files, '~/WINDEX/data_processing/Mission_tadoussac_2009/mergedGPS');
for i = 1:length(latVec)
    m_line(lonVec(i), latVec(i),'marker','o','color',[1 0 0]*.9 ,'markersize',2);
end
m_text(lonVec(1), latVec(1), datestr(timevec(1),15), 'vertical', ...
       'bottom', 'horizontal','right', 'color', [1 0 0]*.9,'FontSize',FONTSIZE2, 'FontWeight', 'bold')
m_text(lonVec(end), latVec(end), datestr(timevec(end),15), 'vertical', ...
       'bottom', 'horizontal','right', 'color', [1 0 0]*.9,'FontSize',FONTSIZE2, 'FontWeight', 'bold')
m_text(lonVec(60), latVec(60), datestr(timevec(62),15), 'vertical', ...
       'top', 'horizontal','right', 'color', [1 0 0]*.9,'FontSize',FONTSIZE2, 'FontWeight', 'bold')


% Thermograph position
m_line(-69.6751, 48.1197, 'marker','p','color',[0 0 1]*.5 ...
       ,'markersize',12, 'markerFaceColor', [0 0 1]*.5);

m_grid('box','fancy', 'yticklabel', [], 'xticklabel', [])

set(gca, 'fontsize', FONTSIZE2)


% --- S2 --- %
subplot(212)


contourf(timevec, zvec, T_mat, 50, 'linestyle', 'none')
hold on
contour(timevec, zvec, T_mat, [1 1], 'k', 'linewidth', 2) 
%contour(timevec, zvec, rho_mat, [1023:.2:1026.6], 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
x = [timeVecTopo(1) timeVecTopo  timeVecTopo(end) timeVecTopo(1)];
y = [zmax botVecTopo zmax zmax];
patch(x, y, [1 1 1]*.6)
hold off

set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.12 0.01 0.35])
ylim([0 max(zvec)])
datetick('x',15)
%ylabel('depth(m)')
xlabel(datestr(timevec(1), 1), 'Fontsize', FONTSIZE)
xlim([t0 tf])
ylabel(cb,'T (^{\circ}C)', 'FontSize', FONTSIZE)
caxis([0 6])
set(gca, 'fontsize', FONTSIZE)
text(datenum(2009, 10, 01, 20, 15, 0), 145, 'b', 'fontweight', 'bold', 'fontsize', FONTSIZE)

adjust_space



set(gcf, 'renderer', 'painters')
print('-depsc2', 'timedepth_pumpPaper_defense.eps');


% ------------------------------ %
% ----- inserted functions ----- %
% ------------------------------ %


function cb = contour_rho(timevec, zvec, rho_mat, V)

contourf(timevec, zvec, rho_mat, 100, 'linestyle', 'none')
hold on
contour(timevec, zvec, rho_mat, V, 'k') 
hold off
set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.15 0.01 0.75])
ylim([0 max(zvec)])
datetick('x',15)
ylabel('depth(m)')


function cb = contour_E(timevec, zvec, E_mat, V)

contourf(timevec, zvec, E_mat, 50, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.15 0.01 0.75])
ylim([0 max(zvec)])
datetick('x',15)
ylabel('depth(m)')

function cb = add_turbulence(timevec, zvec, field, epsMin, epsMax)
colormap(gray)
freezeColors

hold on
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
hold off


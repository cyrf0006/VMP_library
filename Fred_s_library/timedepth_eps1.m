function timedepth_eps1(eps_files, epsMin, epsMax, varargin)

% function timedepth_eps(no_profile, epsMin, epsMax, varargin)
% was previously isopyc_tidecycle.m
% 
% ex: timedepth_eps('eps_files', 1e-9, 1e-5, 90, 1)

% "ls -1 eps_profile*.mat | sed 's/\.mat//' > eps_files"
%
% 
%
%% Modified in March 2016 for CFL data (only use eps1)

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
I = strfind(eps_fname, 'eps_');
data_fname = eps_fname;
data_fname(I:I+3) = [];
load(data_fname);

T1 = [str2num(time(1, 1:2)) str2num(time(1, 4:5)) str2num(time(1, 7:8))]; % [hh mm ss] (num vector)

% find time max
eps_fname = epsfiles(end, :); 
I = find(eps_fname==' '); %remove white space  
eps_fname(I) = [];
I = strfind(eps_fname, 'eps_');
data_fname = eps_fname;
data_fname(I:I+3) = [];
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

% ---  LOOP on profiles ---- %
for profile = 1:no_profile

    eps_fname = epsfiles(profile, :);
    I = find(eps_fname==' '); %remove white space  
    eps_fname(I) = [];
    I = strfind(eps_fname, 'eps_');
    data_fname = eps_fname;
    data_fname(I:I+3) = [];
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
        
    % Combination of the 2 profiles %% Here the change!!!
    p_eps = p_eps1;
    eps = eps1;  
    
    % Diffusivity
    K = 0.2.*eps./(N.^2);
    
    %% Extra cleaning for CFL
    if profile == 53
        I = find(P>100);
        SBT(I) = NaN;
    end
    
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
    I = find(isnan(SBT));
    temp = SBT;
    temp(I)=[];
    PP = P;
    PP(I) = [];    
    tempI = interp1(PP, temp, p_eps);

    % identify bottom
    bot(profile) = max(P);
     
    % bin fluorescence to epsilon
    I = find(isnan(fluoro)==1);
    fluo = fluoro;
    fluo(I)=[];
    
    if ~isempty(fluo)
        PP = p;
        PP(I) = [];  
        fluoI = interp1(PP, fluo, p_eps);
        % bin trans to epsilon
        I = find(isnan(trans)==1);
        trans(I)=[];
        PPP = p;
        PPP(I) = [];  
        transI = interp1(PPP, trans, p_eps);
    else
        fluoI = nan(size(p_eps));
        transI = nan(size(p_eps));
    end

    % time vector
    timevec(profile) = mtime(1);
    
    for i = 1:length(p_eps);
        [Y I] = min(abs(zvec-p_eps(i)));        
        rho_mat(I, profile) = rhoI(i);
        eps_mat(I, profile) = eps(i);
        T_mat(I, profile) = tempI(i);
        K_mat(I, profile) = K(i);
        N_mat(I, profile) = N(i);
        F_mat(I, profile) = fluoI(i);
        TR_mat(I, profile) = transI(i);
    end
    
end 
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

% clean fluorescence 
% (linear itp in time to remove nans)
for i = 1:size(T_mat,1)
    raw_vec = F_mat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        F_mat(i,I(1):I(end)) = vec_itp;
    end
end

% clean trans
% (linear itp in time to remove nans)
for i = 1:size(TR_mat,1)
    raw_vec = TR_mat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        TR_mat(i,I(1):I(end)) = vec_itp;
    end
end


%% Extra cleaning for CFL data



% -------------------- %
% --- Density plot --- %
% -------------------- %
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

% -- density contourplot -- %
V = 1018:0.2:1028;
cb = contour_rho(timevec, zvec, rho_mat, V);
ylabel(cb,'\rho (kg m^{-3})', 'FontSize', 10)
filename = [datestr(timevec(1), 29) '_isopyc.png'];
print('-dpng', '-r300', filename);

% $$$ % Extra isopycnal to highlight displacement
% $$$ hold on 
% $$$ contour(timevec, zvec, rho_mat, [1027 1027], 'linestyle', '-', 'linewidth', 3, 'edgecolor', 'k')
% $$$ hold off
% $$$ filename = [datestr(timevec(1), 29) '_isopyc_xtraiso.png'];
% $$$ print('-dpng', '-r300', filename);



% -- density with dissipation -- %
epsMin = 1e-9;
epsMax = 1e-5;
% $$$ epsMin = 1e-8;
% $$$ epsMax = 1e-5;
colorbar(cb, 'delete');
cb = add_turbulence(timevec, zvec, eps_mat, epsMin, epsMax);
ylabel(cb,'\epsilon (W kg^{-1})', 'FontSize', 10)
xlim([t0 tf])
xlabel(datestr(timevec(1), 1))


hold on
contour(timevec, zvec, rho_mat, V, 'linestyle', '-', 'linewidth', 1, 'edgecolor', 'k')
x = [timevec(1) timevec  timevec(end) timevec(1)];
y = [zmax bot zmax zmax];
patch(x, y, [1 1 1]*.6)


filename = [datestr(timevec(1), 29) '_isopyc_eps_xtraiso.png'];
print('-dpng', '-r300', filename);
% $$$ filename = [datestr(timevec(1), 29) '_isopyc_eps_xtraiso.eps'];
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', filename);

%print('-dpng', '-r300','isopyc_eps_xtraiso.png')
% $$$ 
% $$$ 
% -- density with diffusivity -- %
figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

cb = contour_rho(timevec, zvec, rho_mat, V);

KMin = 1e-7;
KMax = 1e-3;
colorbar(cb, 'delete');
cb = add_turbulence(timevec, zvec, K_mat, KMin, KMax);
ylabel(cb,'K (m^2 s^{-1})', 'FontSize', 10)
xlim([t0 tf])
xlabel(datestr(timevec(1), 1))
filename = [datestr(timevec(1), 29) '_isopyc_K.png'];
print('-dpng', '-r300', filename);
%print('-dpng', '-r300','isopyc_K.png')
% $$$ 
% $$$ 
% $$$ % ----- Raw buoy. freq. ----- %
% $$$ figure(9)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% $$$ 
% $$$ cb = contour_rho(timevec, zvec, rho_mat, V);
% $$$ 
% $$$ KMin = 1e-8;
% $$$ KMax = 1e-1;
% $$$ colorbar(cb, 'delete');
% $$$ cb = add_turbulence(timevec, P_N2, abs(N2raw_mat), KMin, KMax);
% $$$ ylabel(cb,'N^2 (s^{-2})', 'FontSize', 10)
% $$$ xlim([t0 tf])
% $$$ xlabel(datestr(timevec(1), 1))
% $$$ filename = [datestr(timevec(1), 29) '_isopyc_N2raw.png'];
% $$$ print('-dpng', '-r300', filename);
% $$$ 
% $$$ 
% $$$ figure(5)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% $$$ % -- fluorescence contourplot -- %
% $$$ V = 1018:0.2:1028;
% $$$ %cb = contour_F(timevec, zvec, F_mat);
% $$$ imagesc(timevec, zvec, F_mat)
% $$$ cb = colorbar;
% $$$ ylabel(cb,'fluorescence (ppb)', 'FontSize', 10)
% $$$ colormap(jet)
% $$$ xlim([t0 tf])
% $$$ hold on
% $$$ patch(x, y, [1 1 1]*.6)
% $$$ hold off
% $$$ filename = [datestr(timevec(1), 29) '_fluo.png'];
% $$$ print('-dpng', '-r300', filename);

% $$$ figure(5)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% $$$ % -- trans contourplot -- %
% $$$ V = 1018:0.2:1028;
% $$$ cb = contour_F(timevec, zvec, TR_mat);
% $$$ ylabel(cb,'trans (???)', 'FontSize', 10)
% $$$ colormap(jet)
% $$$ xlim([t0 tf])
% $$$ filename = [datestr(timevec(1), 29) '_trans.png'];
% $$$ print('-dpng', '-r300', filename);
% $$$ 
% $$$ % $$$ 
% $$$ % $$$ figure(6)
% $$$ % $$$ clf
% $$$ % $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% $$$ % $$$ 
% $$$ % $$$ % -- epsilon contourplot -- %
% $$$ % $$$ cb = contour_E(timevec, zvec, log10(eps_mat));
% $$$ % $$$ ylabel(cb,'\epsilon (W/kg)', 'FontSize', 10)
% $$$ % $$$ colormap(jet)
% $$$ % $$$ filename = [datestr(timevec(1), 29) '_eps.png'];
% $$$ % $$$ print('-dpng', '-r300', filename);
% $$$ % $$$ 
% $$$ % $$$ figure(7)
% $$$ % $$$ clf
% $$$ % $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% $$$ % $$$ 
% $$$ % $$$ % -- epsilon contourplot -- %
% $$$ % $$$ cb = contour_F(timevec, zvec, TR_mat);
% $$$ % $$$ ylabel(cb,'trans (FTU)', 'FontSize', 10)
% $$$ % $$$ colormap(jet)
% $$$ % $$$ filename = [datestr(timevec(1), 29) '_trans.png'];
% $$$ % $$$ print('-dpng', '-r300', filename);
% $$$ 
% $$$ figure(8)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
% $$$ 
% $$$ % -- epsilon contourplot -- %
% $$$ %cb = contour_E(timevec, zvec, log10(N_mat.^2));
% $$$ %cb = contourf(timevec, zvec, log10(N_mat.^2), 'linestyle', 'none');
% $$$ cb = imagesc(timevec, zvec, log10(N_mat.^2));
% $$$ set(gca, 'ydir', 'reverse')
% $$$ cb = colorbar;
% $$$ set(cb, 'pos',  [0.92 0.15 0.01 0.75]);
% $$$ ylim([0 max(zvec)]);
% $$$ datetick('x',15);
% $$$ ylabel('depth(m)');
% $$$ ylabel(cb,'N^2 (s^{-1})', 'FontSize', 10);
% $$$ colormap(jet)
% $$$ filename = [datestr(timevec(1), 29) '_N2.png'];
% $$$ print('-dpng', '-r300', filename);
% $$$ 
% ------------------------ %
% --- Temperature plot --- %
% ------------------------ %
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

V = [];
cb = contour_T(timevec, zvec, T_mat, V);
xlim([t0 tf])
xlabel(datestr(timevec(1), 1))
ylabel(cb,'T (^{\circ}C)', 'FontSize', 10)
% $$$ load('/home/cyrf0006/PhD/CTD_IML4/BIN/renamed/BRcolormap3.mat')
% $$$ colormap(mycolormap)
hold on    
contour(timevec, zvec, T_mat, [1 1], 'edgecolor', 'k', 'LineStyle', '-','linewidth', 2 ) 
hold off
filename = [datestr(timevec(1), 29) '_temp.png'];
%load BR_colormap
%colormap() 
%caxis([-1 5])
set(gca, 'fontsize', 15)
patch(x, y, [1 1 1]*.6)
colormap(jet)
keyboard
print('-dpng', '-r300', filename);
filename = [datestr(timevec(1), 29) '_temp.eps'];
set(gcf, 'renderer', 'painters')
print('-depsc2', filename);
%print('-dpng', '-r300','temp.png')
% $$$ 
% $$$ % -- Temperature with dissipation -- %
% $$$ epsMin = 1e-9;
% $$$ epsMax = 1e-5;
% $$$ % $$$ epsMin = 1e-8;
% $$$ % $$$ epsMax = 1e-5;
% $$$ colorbar(cb, 'delete');
% $$$ cb = add_turbulence(timevec, zvec, eps_mat, epsMin, epsMax);
% $$$ xlim([t0 tf])
% $$$ xlabel(datestr(timevec(1), 1))
% $$$ cbti = ylabel(cb,'\epsilon (W kg^{-1})', 'FontSize', 15);
% $$$ ypos = get(cbti, 'pos');
% $$$ ypos(1) = 13;
% $$$ set(cbti, 'pos', ypos)  
% $$$ 
% $$$ hold on    
% $$$ contour(timevec, zvec, T_mat, [1 1], 'edgecolor', 'k', 'LineStyle', '-','linewidth', 2 ) 
% $$$ hold off
% $$$ filename = [datestr(timevec(1), 29) '_temp_eps.png'];
% $$$ print('-dpng', '-r300', filename);
% $$$ %print('-dpng', '-r300','temp_eps.png')


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

function cb = contour_T(timevec, zvec, T_mat, V)

contourf(timevec, zvec, T_mat, 50, 'linestyle', 'none')
hold on
contour(timevec, zvec, T_mat, V, 'k') 
hold off
set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.15 0.01 0.75])
ylim([0 max(zvec)])
datetick('x',15)
ylabel('depth(m)')

function cb = contour_F(timevec, zvec, F_mat, V)

contourf(timevec, zvec, F_mat, 50, 'linestyle', 'none')
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




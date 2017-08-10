function timedepth_chi(chi_files, epsMin, epsMax, varargin)

% function timedepth_chi(chi_files, epsMin, epsMax, varargin)
% 
% ex: timedepth_chi('chi_files', 1e-9, 1e-5, 50, 1)

% "ls -1 chi_profile*.mat | sed 's/\.mat//' > chi_files"
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

% load eps_files file
fid = fopen(chi_files);
C = textscan(fid, '%s', 'delimiter', '\n');
chifiles = char(C{1});

no_profile = size(chifiles, 1); %number of *.P files 

% find time min
chi_fname = chifiles(1, :); 
I = find(chi_fname==' '); %remove white space  
chi_fname(I) = [];
data_fname = chi_fname(5:end);
load(data_fname);

T1 = [str2num(time(1, 1:2)) str2num(time(1, 4:5)) str2num(time(1, 7:8))]; % [hh mm ss] (num vector)

% find time max
chi_fname = chifiles(end, :); 
I = find(chi_fname==' '); %remove white space  
chi_fname(I) = [];
data_fname = chi_fname(5:end);
load(data_fname);
sizeT = size(time);
T2 = [str2num(time(sizeT(1), 1:2)) str2num(time(sizeT(1), 4:5)) str2num(time(sizeT(1), 7:8))]; % [hh mm ss] (num vector
dat = date(1,:);
dd = dat(1:2); 
mm = dat(4:5); 
yyyy = dat(7:10);


%zbin=1;
zvec = [zmin:zbin:zmax]';
rho_mat = nan(size(zvec,1), no_profile);
chi_mat = nan(size(zvec,1), no_profile); 
T_mat = nan(size(zvec,1), no_profile);
K_mat = nan(size(zvec,1), no_profile);
F_mat = nan(size(zvec,1), no_profile);


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

    chi_fname = chifiles(profile, :); 
    I = find(chi_fname==' '); %remove white space  
    chi_fname(I) = [];
    data_fname = chi_fname(5:end);
    disp(data_fname);

    load(data_fname);
    load(chi_fname);
    
    % time of survey
    if profile==1
        t0 = mtime(1);
    elseif profile==no_profile
        tf = mtime(end);
    end
    
% $$$     
% $$$     % Test for Fred (pas sur que dans tous les fichiers p_eps1 et p_eps2 sont les memes...)
% $$$     if length(p_eps1)~=length(p_eps2) | ( p_eps1(1)~=p_eps2(1) & ...
% $$$                                           isnan(p_eps2(1))==0 )
% $$$         disp('p_eps1 and p_eps2 mismatch')
% $$$         pause
% $$$     end
        
    % Combination of the 2 profiles
% $$$     p_eps = nanmean([p_eps1; p_eps2]);
% $$$     eps = nanmean([eps1; eps2]);  
% $$$          
% $$$     p_chi = nanmean([p_chi1; p_chi2]);
% $$$     chi = nanmean([chi1; chi2]);  
    chi = chi1;       
    p_chi = p_chi1;
    
    % density 
    rho = sw_dens(SBS, SBT, P);
    % bin density to epsilon
    I = find(isnan(rho)==1);
    rho(I)=[];
    PP = P;
    PP(I) = [];    
    rhoI = interp1(PP, rho, p_chi);
    
    
    % bin temperature to epsilon
    I = find(isnan(SBT)==1);
    temp = SBT;
    temp(I)=[];
    PP = P;
    PP(I) = [];    
    tempI = interp1(PP, temp, p_chi);
    
% $$$     % bin fluorescence to epsilon
% $$$     I = find(isnan(fluoro)==1);
% $$$     fluo = fluoro;
% $$$     fluo(I)=[];
% $$$     PP = p;
% $$$     PP(I) = [];  
% $$$     if length(PP) == length(fluo)
% $$$         fluoI = interp1(PP, fluo, p_chi);
% $$$     else
% $$$         fluoI = nan(size(p_chi));
% $$$     end
% $$$      
% $$$     % bin trans to epsilon
% $$$     I = find(isnan(trans)==1);
% $$$     trans(I)=[];
% $$$     PPP = p;
% $$$     PPP(I) = [];  
% $$$     if length(PPP) == length(trans)
% $$$         transI = interp1(PPP, trans, p_chi);
% $$$     else
% $$$         transI = nan(size(p_eps));
% $$$     end
% $$$     
    
    % time vector
    timevec(profile) = mtime(1);
    
    for i = 1:length(p_chi);
        [Y I] = min(abs(zvec-p_chi(i)));        
        rho_mat(I, profile) = rhoI(i);
        chi_mat(I, profile) = chi(i);
        T_mat(I, profile) = tempI(i);
% $$$ % $$$         K_mat(I, profile) = K(i);
% $$$         N_mat(I, profile) = N(i);
% $$$         F_mat(I, profile) = fluoI(i);
% $$$         TR_mat(I, profile) = transI(i);
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


% ---- For mixing efficiency calculation ---- %
timevec_chi = timevec;
zvec_chi = zvec;
save MixEff_chi.mat chi_mat T_mat timevec_chi zvec_chi
% ------------------------------------------- %


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
%print('-dpng', '-r300', filename);
%print('-dpng', '-r300','isopyc.png')

% Extra isopycnal to highlight displacement
hold on 
contour(timevec, zvec, rho_mat, [1027 1027], 'linestyle', '-', 'linewidth', 3, 'edgecolor', 'k')
hold off
filename = [datestr(timevec(1), 29) '_isopyc_xtraiso.png'];
%print('-dpng', '-r300', filename);
%print('-dpng', '-r300','isopyc_xtrasiso.png')




% -- density with dissipation -- %
% $$$ epsMin = 1e-9;
% $$$ epsMax = 1e-5;
% $$$ epsMin = 1e-8;
% $$$ epsMax = 1e-5;
colorbar(cb, 'delete');
cb = add_turbulence(timevec, zvec, chi_mat, epsMin, epsMax);
ylabel(cb,'\chi (^{\circ}C^2 s^{-1})', 'FontSize', 10)
xlim([t0 tf])
xlabel(datestr(timevec(1), 1))
filename = [datestr(timevec(1), 29) '_isopyc_chi_xtraiso.png'];
print('-dpng', '-r300', filename);
%print('-dpng', '-r300','isopyc_eps_xtraiso.png')


% ------------------------ %
% --- Temperature plot --- %
% ------------------------ %
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 20 10])

V = [];
cb = contour_T(timevec, zvec, T_mat, V);
xlim([t0 tf])
xlabel(datestr(timevec(1), 1))
ylabel(cb,'T (^{\circ}C)', 'FontSize', 10)

hold on    
contour(timevec, zvec, T_mat, [1 1], 'edgecolor', 'k', 'LineStyle', '-','linewidth', 2 ) 
hold off
filename = [datestr(timevec(1), 29) '_temp.png'];
caxis([-1 5])
set(gca, 'fontsize', 15)
filename = [datestr(timevec(1), 29) '_temp.png'];
print('-dpng', '-r300', filename);

% -- Temperature with dissipation -- %
% $$$ epsMin = 1e-9;
% $$$ epsMax = 1e-5;
% $$$ epsMin = 1e-8;
% $$$ epsMax = 1e-5;
colorbar(cb, 'delete');
cb = add_turbulence(timevec, zvec, chi_mat, epsMin, epsMax);
xlim([t0 tf])
xlabel(datestr(timevec(1), 1))
cbti = ylabel(cb,'\chi (^{\circ}C^2 s^{-1})', 'FontSize', 15);
ypos = get(cbti, 'pos');
ypos(1) = 12;
set(cbti, 'pos', ypos)  

keyboard
hold on    
contour(timevec, zvec, T_mat, [1 1], 'edgecolor', 'k', 'LineStyle', '-','linewidth', 2 ) 
hold off
filename = [datestr(timevec(1), 29) '_temp_chi.png'];
print('-dpng', '-r300', filename);



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
    data = field(:,j);
    x(1:length(data)) = timevec(j);

    I = find(~isnan(data)==1); % ignore Nans
    data = data(I);
    y = zvec(I);    
    x = x(I);

    rectbar = [0.92 0.15 0.01 0.75]; 
    if j == size(field,2) % with cbar
        cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
    end

 
end
hold off
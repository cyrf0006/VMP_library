function timedepth_eps_hab(eps_files, epsMin, epsMax, Pavg)

% function timedepth_eps(no_profile, epsMin, epsMax, varargin)
% was previously isopyc_tidecycle.m
% 
% ex: timedepth_eps('eps_files', 1e-9, 1e-5, 50, 1)
%
% "ls -1 eps_profile*.mat | sed 's/\.mat//' > eps_files"
%
% 
%


zmin = 0;
zmax = 120;
zbin = 1;
Pbin = [zmin+zbin/2:zbin:zmax-zbin/2]';
hab = Pbin;
%Pavg = 20;  % last Pavg meters of the water column average

g = 9.81;

% load eps_files file
fid = fopen(eps_files);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

no_profile = size(epsfiles, 1); %number of *.P files 

% find time min
fname = epsfiles(1, :);
I = find(fname==' '); %remove white space  
fname(I) = [];
load(fname)
T1 = mtime_eps(1);

% find time max
fname = epsfiles(end, :); 
I = find(fname==' '); %remove white space  
fname(I) = [];
load(fname)
T2 = mtime_eps(1);

% raw matrix to fill
epsMat = nan(length(Pbin), no_profile);
NMat = nan(length(Pbin), no_profile);
RMat = nan(length(Pbin), no_profile);

fid = fopen('./Matrices.mat');
if fid == -1 % doesnt exist

    for i = 1:no_profile
        
        disp(sprintf('profile %d', i))
        fname = epsfiles(i, :);
        I = find(fname==' ');   
        fname(I) = [];    
        load(fname)
        
        fname = fname(5:end);
        I = find(fname==' ');   
        fname(I) = [];    
        load(fname)
        
        epsilon = nanmean([eps1', eps2'], 2);
        rho = sw_dens(SBS, SBT, P);

        
        % Bin profiles
        epsbin = nan(length(Pbin),1);
        Nbin = nan(length(Pbin),1);
        Rbin = nan(length(Pbin),1);
        
        for j = 1:length(Pbin)
            I = find(p_eps1 >= (Pbin(j) - zbin/2) & p_eps1 <= (Pbin(j) + zbin/2));
            if ~isempty(I) == 1
                epsbin(j) = nanmean(epsilon(I));
                Nbin(j) = nanmean(N(I));
            end
            I = find(P >= (Pbin(j) - zbin/2) & P <= (Pbin(j) + zbin/2));
            if ~isempty(I) == 1
                Rbin(j) = nanmean(rho(I));
            end        
        end
        
        % rename variables and reference relative to hab
        epsilon = epsbin;
        N = Nbin;
        R = Rbin;
        I = find(~isnan(N)==1); % we take N on purpose
        epsilon(1:I(end)) = flipud(epsilon(1:I(end))); % now hab
        N(1:I(end)) = flipud(N(1:I(end))); % now hab
        R(1:I(end)) = flipud(R(1:I(end))); % now hab
        NMat(:, i) = N;
        epsMat(:, i) = epsilon;
        RMat(:, i) = R;
        proftime(i) = mtime(1);    
    end 
    save Matrices.mat epsMat RMat NMat proftime

else
    load Matrices
end



timevec = proftime;
% clean density field 
% (linear itp in time to remove nans)
for i = 1:size(RMat,1)
    raw_vec = RMat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        RMat(i,I(1):I(end)) = vec_itp;
    end
end

% clean N2 field 
% (linear itp in time to remove nans)
for i = 1:size(NMat,1)
    raw_vec = NMat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        NMat(i,I(1):I(end)) = vec_itp;
    end
end


% Reduce vertical size and rename for plotting (copy-paste)
I = find(hab<=Pavg);
zvec = hab(I);
rho_mat = RMat(I,:);
eps_mat = epsMat(I,:);
N_mat = NMat(I,:);



% -------------------- %
% --- Density plot --- %
% -------------------- %
t0 = T1-1/24;
tf = T2+1/24;
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])

% -- density contourplot -- %
V = 1018:1:1028;
cb = contour_rho(timevec, zvec, rho_mat, V);
ylabel(cb,'\rho (kg m^{-3})', 'FontSize', 10)
set(gca, 'ydir', 'normal')
filename = ['hab_' datestr(timevec(1), 29) '_isopyc.png'];
print('-dpng', '-r300', filename);
%print('-dpng', '-r300','isopyc.png')

% $$$ % Extra isopycnal to highlight displacement
% $$$ hold on 
% $$$ contour(timevec, zvec, rho_mat, [1025 1025], 'linestyle', '-', 'linewidth', 3, 'edgecolor', 'k')
% $$$ hold off
% $$$ filename = ['hab_' datestr(timevec(1), 29) '_isopyc_xtraiso.png'];
% $$$ print('-dpng', '-r300', filename);
% $$$ %print('-dpng', '-r300','isopyc_xtrasiso.png')



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
filename = ['hab_' datestr(timevec(1), 29) '_isopyc_eps.png'];
print('-dpng', '-r300', filename);

set(gcf, 'renderer', 'painters')
filename = ['hab_' datestr(timevec(1), 29) '_isopyc_eps.eps'];
print('-depsc2', filename)
%print('-dpng', '-r300','isopyc_eps_xtraiso.png')

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
% $$$ filename = ['hab_' datestr(timevec(1), 29) '_isopyc_N2raw.png'];
% $$$ print('-dpng', '-r300', filename);


% ------------------------------ %
% ----- inserted functions ----- %
% ------------------------------ %


function cb = contour_rho(timevec, zvec, rho_mat, V)

contourf(timevec, zvec, rho_mat, 100, 'linestyle', 'none')
% $$$ hold on
% $$$ contour(timevec, zvec, rho_mat, V, 'k') 
% $$$ hold off
set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.15 0.01 0.75])
ylim([0 max(zvec)])
datetick('x',15)
ylabel('hab(m)')

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
ylabel('hab(m)')

function cb = contour_F(timevec, zvec, F_mat, V)

contourf(timevec, zvec, F_mat, 50, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.15 0.01 0.75])
ylim([0 max(zvec)])
datetick('x',15)
ylabel('hab(m)')

function cb = contour_E(timevec, zvec, E_mat, V)

contourf(timevec, zvec, E_mat, 50, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
cb = colorbar;
set(cb, 'pos',  [0.92 0.15 0.01 0.75])
ylim([0 max(zvec)])
datetick('x',15)
ylabel('hab(m)')

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




function gamma_M2(tfiles, epsfiles, chifiles, Pinfo, zhab, tidefile, varargin)

%function gamma_N2(epsfile, chifiles, zbin, zhab, tidefile, varargin)

% usage ex (in /home/cyrf0006/WINDEX/data_processing/BMix_study): 
%    gamma_M2('hit_bottom_tprofiles', 'hit_bottom_epsprofiles','hit_bottom_chiprofiles',[0 200 5], 20, 'tide_2009-2012.dat')
%    gamma_M2('hit_bottom_tprofiles_with2010',
%    'hit_bottom_epsprofiles_with2010','hit_bottom_chiprofiles_with2010',[0
%    200 1], 20, 'tide_2009-2012.dat', 'nochi')

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deal varargin
if size(varargin, 2) > 2
    disp('Wrong input')
elseif strcmp(varargin, 'nochi') == 1
    skip = 1;
else
    skip = 0;
end

    

% depth range for extracting the average
g = 9.81;
zmin = Pinfo(1); % top of the fisrt bin
zmax = Pinfo(2);
zbin = Pinfo(3);
nboot = 500;    
Pbin = [zmin+zbin/2:zbin:zmax-zbin/2]';
hab = Pbin;
% load *.P files names (file in which are recorded .P files)
fid = fopen(tfiles);
C = textscan(fid, '%s', 'delimiter', '\n');
tfiles = char(C{1});

fid = fopen(epsfiles);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

fid = fopen(chifiles);
C = textscan(fid, '%s', 'delimiter', '\n');
chifiles = char(C{1});


no_profile = size(epsfiles, 1); %number of eps_files 

% raw matrix to fill
gamma_mat = nan(length(Pbin), no_profile);
% ----------- Bring all profile at same vertical resolution --------------- %

for i = 1:no_profile

    
    disp(sprintf('profile %d', i))
    fname = tfiles(i, :);
    I = find(fname==' ');   
    fname(I) = [];    
    load(fname)
    
    if skip == 0
        fname = chifiles(i, :);
        I = find(fname==' ');   
        fname(I) = [];    
        load(fname)
    end
    
    
    % WATCH OUT! 
    % eps1 and eps2 exist in both epsfiles and chifiles!
    % (epsfiles must be loaded after chifiles)
    fname = epsfiles(i, :);
    I = find(fname==' ');   
    fname(I) = [];    
    load(fname)
    
    %maxP(i) = max(P);
    % basic test if a thermistor is missing or broken
    if nochi == 0
        dt1dz = diff(t1)./diff(p);
        dt2dz = diff(t2)./diff(p);
        if nanstd(dt2dz) > 2*nanstd(dt1dz);
            chi2 = chi2*NaN;
        end
        if nanstd(dt1dz) > 2*nanstd(dt2dz);
            chi1 = chi1*NaN;
        end
        
        chi = nanmean([chi1', chi2'], 2);
        I = find(chi<=0);
        chi(I) = NaN;
    else
        chi = eps1*NaN;
    end
    
    epsilon = nanmean([eps1', eps2'], 2);
    
    
    % Bin profiles
    epsbin = nan(length(Pbin),1);
    chibin = nan(length(Pbin),1);
    Nbin = nan(length(Pbin),1);
    Tbin = nan(length(Pbin),1);
    dtdzbin = nan(length(Pbin),1);    
    for j = 1:length(Pbin)
        I = find(p_eps1 >= (Pbin(j) - zbin/2) & p_eps1 <= (Pbin(j) + zbin/2));
        if ~isempty(I) == 1
            epsbin(j) = nanmean(epsilon(I));
            chibin(j) = nanmean(chi(I));
            Nbin(j) = nanmean(N(I));
        end
        
        
        I = find(P >= (Pbin(j) - zbin/2) & P <= (Pbin(j) + zbin/2));
        if ~isempty(I) == 1
            Tbin(j) = nanmean(SBT(I));
            p = polyfit(P(I),SBT(I),1);
            dtdzbin(j) = p(1);
        end
    end
    
  
    % rename variables
    epsilon = epsbin;
    chi = chibin;
    N = Nbin;
    T = Tbin;
    dtdz = dtdzbin;

    % Compute Gamma
    % HERE SHOULD I REORDER THE PROFILE FIRST???
    %dtdz = gradient(T)./dp;
    gamma = N.^2.*chi./(2*epsilon.*dtdz.^2);
    
    I = find(isnan(gamma)==1);
    gamma(I) = [];
    dtdz(I) = [];
    T(I) = [];
    N(I) = [];
    epsilon(I) = [];
    chi(I) = [];
    
    gamma_mat(1:length(gamma), i) = flipud(gamma); %fliplr to reference to hab
    dtdz_mat(1:length(gamma), i) = flipud(dtdz); %fliplr to reference to hab
    T_mat(1:length(gamma), i) = flipud(T); %fliplr to reference to hab
    N_mat(1:length(gamma), i) = flipud(N); %fliplr to reference to hab
    eps_mat(1:length(gamma), i) = flipud(epsilon); %fliplr to reference to hab
    chi_mat(1:length(gamma), i) = flipud(chi); %fliplr to reference to hab
    proftime(i) = mtime(1);
end 

save turbulence_M2.mat gamma_mat dtdz_mat N_mat eps_mat chi_mat proftime


% ---------------------------------------------------- %
keyboard
% -------------- Compute time to high tide ------------ %
tide  = load(tidefile);
    
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);

% find high tide time 
count = 1;
clear T
for i = 2:length(mtime)-1

    if level(i)>level(i-1) & level(i)>level(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        count = count+1;
    end
end
% T contains the hour of each high tide
clear mtime

for i = 1:length(proftime)
    [Y, I] = min(abs(proftime(i)-T));
    A(i) = (proftime(i)-T(I))*24;
    B(i) = L(I); %level of the closest hightide
end
% ------------------------------------------------------ %

time2 = A;

dtide = 0.25;
dtide = 1;
reg_tide = -6:dtide:6;


J = find(hab<=zhab);
gamma_tide = nan(length(hab(J)), length(reg_tide));
dtdz_tide = nan(length(hab(J)), length(reg_tide));
T_tide = nan(length(hab(J)), length(reg_tide));
N_tide = nan(length(hab(J)), length(reg_tide));
eps_tide = nan(length(hab(J)), length(reg_tide));
chi_tide = nan(length(hab(J)), length(reg_tide));


for i = 1: length(reg_tide)
    I = find(time2 > reg_tide(i) - dtide & time2 < reg_tide(i) + dtide);
    gamma_tide(:,i) = nanmean(gamma_mat(J,I), 2);
    dtdz_tide(:,i) = nanmean(dtdz_mat(J,I), 2);
    T_tide(:,i) = nanmean(T_mat(J,I), 2);
    N_tide(:,i) = nanmean(N_mat(J,I), 2);
    eps_tide(:,i) = nanmean(eps_mat(J,I), 2);
    chi_tide(:,i) = nanmean(chi_mat(J,I), 2);
end
    
figure(1)
clf
%V = [-5:0.01:-3];
imagesc(reg_tide, hab(J), gamma_tide);%, V, 'linestyle', 'none');
xlabel('time to high tide (hours)')
ylabel('hab (m)')
c= colorbar;
caxis([0 .2])
title('\Gamma (s^{-2})') 
set(gca, 'ydir', 'normal')
ti = ylabel(c,'{\Gamma}', 'FontSize', 8);

keyboard

% $$$ 
% $$$ % LOESS-type test
% $$$ N2_meanvec = nanmean(N2_mat(J, :));
% $$$ I = find(~isnan(N2_meanvec)==1);
% $$$ smooth_N2 = loess(time2(I), N2_meanvec(I), sort(time2(I)), 0.3, 1);
% $$$ [time_N2, II] = sort(time2(I));
% $$$ 
% $$$ save N2_M2.mat N2_2p5 N2_97p5 N2_ave N2_error smooth_N2 time_N2
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ imagesc(reg_tide, hab, N2_tide)
% $$$ set(gca, 'ydir', 'normal') 
% $$$ ylim([0 50])
% $$$ c = colorbar;
% $$$ caxis([-5 -3])
% $$$ xlabel('time relative to high tide')
% $$$ ylabel('hab')
% $$$ ti = ylabel(c,'log(N^2) (s^{-2})', 'FontSize', 8);
% $$$ 
% $$$ figure(2)
% $$$ clf
% $$$ errorbar(reg_tide, N2_ave, N2_error, '.k')
% $$$ ylabel('N^2 (s^{-2})')
% $$$ xlabel('time relative to high tide')
% $$$ set(gca, 'yscale', 'log') 
% $$$ xlim([-7 7]) 
% $$$ hold on
% $$$ plot(time_N2, smooth_N2, 'k', 'linewidth', 2)
% $$$ 







% bootstrap on Gamma
J = find(hab<=zhab);
[TIME2, HAB] = meshgrid(time2, hab(J));
GAM = gamma_mat(J,:);
I = find(isnan(GAM)==1);
GAM(I)=[];
TIME2(I) = [];
nboot = 1000;
N = length(GAM);

% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    gamma_boot_b(b) = nanmean(log10(GAM(r)));
end

% mean of random sampling
gamma_boot_dot = nanmean(gamma_boot_b);

% compute error (EFRON & GONG 83)
% $$$ eps_error = sqrt(sum((diff([eps_boot_b eps_boot_dot],1, 2)).^2, 2)./(nboot-1));

% Compute 95% confidence interval
gamma_sort = sort(gamma_boot_b);
CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

gamma_2p5 = gamma_sort(CI_2p5);
gamma_97p5 = gamma_sort(CI_97p5);

[10^gamma_boot_dot 10^gamma_2p5 10^gamma_97p5]




% GAMMA relative to M2 (in log domain!)

gamma_tide = nan(3, length(reg_tide));


for i = 1: length(reg_tide)
    I = find(TIME2 > reg_tide(i) - dtide & TIME2 < reg_tide(i) + dtide);
    gamvec = GAM(I);
    N = length(gamvec);

    % create random sampling
    for b = 1:nboot
        r = rand(N,1);
        r = ceil(r*N/1);
        gamma_boot_b(b) = nanmean(log10(gamvec(r)));
    end
    gamma_tide(1,i) = nanmean(gamma_boot_b);
    gamma_sort = sort(gamma_boot_b);
    CI_2p5 = round(2.5/100*nboot);
    CI_97p5 = round(97.5/100*nboot);
    gamma_tide(2,i) = gamma_sort(CI_2p5);
    gamma_tide(3,i) = gamma_sort(CI_97p5);
    
end


figure(2)
clf
errorbar(reg_tide, 10.^gamma_tide(1,:), 10.^gamma_tide(2,:), ...
         10.^gamma_tide(3,:), 'k')
xlabel('Time to hightide (hour)')
ylabel('\Gamma')


keyboard
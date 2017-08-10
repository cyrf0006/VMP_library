function N2_M2composite(file_names, zbin, zhab, tidefile, varargin)

% function N2_M2composite(file_names, zbin, tidefile)
%
%
% where: - 
%
% usage ex: N2_M2composite('hit_bottom_2011_3', 0.25, 20, 'tide_2009-2011.dat')
%                                     or   
%           N2_M2composite('hit_bottom_2011_3', 1, 20, 'tide_2009-2011.dat', 'skipload')  
%          
% in /home/cyrf0006/WINDEX/data_processing/BMix_study:
% N2_M2composite('hit_bottom_tprofiles', 1, 20, 'tide_2009-2012.dat', 'skipload')
% N2_M2composite('hit_bottom_tprofiles_with2010', 1, 20, 'tide_2009-2012.dat')
% N2_M2composite('hit_bottom_tprofiles_slope', 1, 20, 'tide_2009-2012.dat')

%
% author: F. Cyr - 2012/05/09
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deal varargin
if size(varargin, 2) > 2
    disp('Wrong input')
elseif strcmp(varargin, 'skipload') == 1
    skip = 1;
else
    skip = 0;
end

    

% depth range for extracting the average
g = 9.81;
zmin = 0; % top of the fisrt bin
zmax = 350;
nboot = 500;    
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';
hab = P_bin;
% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

no_profile = size(epsfiles, 1); %number of eps_files 

% raw matrix to fill
N2_mat = nan(length(P_bin), no_profile);
Ld_mat = nan(length(P_bin), no_profile);

% ----------- Loop on profiles to get N2(hab) --------------- %
if skip == 0
for profile = 1:no_profile

    disp(sprintf('profile %d', profile))
    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
        
    rho = sw_dens(SBS, SBT, P);
    
    
    % flagg bad values (remove values where W > 10 x std(W))
    I = find(~isnan(rho)==1);
    J = find(P(I) > max(P(I))-0.25);
    rho(I(J)) = NaN;
    I = find(W > nanmean(W)-10*nanstd(W) & W < nanmean(W)+10*nanstd(W));
    rho = rho(I);
    P = P(I);
        
    % Raw buoy. freq.
    clear N2_vec Ld_vec
    
    % N2 calculation
    for i = 1:length(P_bin)
        I = find(P >= (P_bin(i) - zbin/2) & P < (P_bin(i) + zbin/2)); 
        if length(I) > 5
            %pol = polyfit(P(I),rho(I),1); %unsorted
            pol = polyfit(P(I),sort(rho(I)),1); %sorted
            rho0 = mean(rho(I));
            N2_vec(i) = (g/rho0)*pol(1);
        else
            N2_vec(i) = NaN;
        end
    end
    I = find(~isnan(N2_vec)==1);
    N2_mat(1:length(I), profile) = N2_vec(fliplr(I)); %fliplr to reference to hab
    proftime(profile) = mtime(1);
    
    % L_D calculation
    II = find(~isnan(rho)==1);
    % interp to remove NaNs inside profile
    rho = interp1(P(II), rho(II), P);
    I = find(isnan(rho) ==1 );
    rho(I) = [];
    P(I) = [];
    [Y, I] = sort(rho);
    Ld = P - P(I);
    for i = 1:length(P_bin)
        I = find(P >= (P_bin(i) - zbin/2) & P < (P_bin(i) + zbin/2));
        Ld_vec(i) = nanmean(Ld(I));
    end
    I = find(~isnan(Ld_vec)==1);
    Ld_mat(1:length(I), profile) = Ld_vec(fliplr(I)); 
    
    zbottom(profile) = max(P); %used in APEF_M2.m
end 

save N2.mat N2_mat Ld_mat hab zbottom proftime %also used in APEF_M2.m
clear P_bin
else
    disp('Quick method, skip loading each profile')
    load N2.mat2
end
% ---------------------------------------------------- %

% -------------- Compute time to high tide ------------ %
tide  = load(tidefile);
    
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);

% find high tide time 
count = 1;
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
N2_tide = nan(length(hab(J)), length(reg_tide));
Ld_tide = nan(length(hab(J)), length(reg_tide));

clear  N2_2p5  N2_97p5 N2_ave N2_error

for i = 1: length(reg_tide)
    I = find(time2 > reg_tide(i) - dtide & time2 < reg_tide(i) + dtide);
    
    
    N2_tide(:,i) = nanmean(N2_mat(J,I), 2);
    Ld_tide(:,i) = nanmean(Ld_mat(J,I), 2);
end
    
figure(1)
clf
V = [-5:0.01:-3];
contourf(reg_tide, hab(J), log10(N2_tide), V, 'linestyle', 'none');
xlabel('time to high tide (hours)')
ylabel('hab (m)')
colorbar   
title('N2 (s^{-2})') 

figure(2)    
clf
V = 0:.01:.5;
contourf(reg_tide, hab(J), abs(Ld_tide), V, 'linestyle', 'none')
xlabel('time to high tide (hours)')
ylabel('hab (m)')
colorbar
title('Ld (m)')    

keyboard


    % ---- tests for negative N2 ---- %
% $$$     II = find(N2_mat(:,I)<0);
% $$$     N2_tide(:,i) = nanmean(log10(abs(N2_mat(:,I))), 2);
% $$$     
% $$$     test_mat = log10(abs(N2_mat(:,I)));
% $$$     one_mat = ones(size(test_mat));
% $$$     one_mat(II) = -1;
% $$$     test_mat = test_mat.*one_mat; 
% $$$     N2_tide2(:,i) = nanmean(test_mat, 2);
% $$$     
% $$$     % bootstrap 
% $$$     NN2 = N2_mat(J,I);
% $$$     NN2 = NN2(:);
% $$$     I = find(isnan(NN2)==1);
% $$$     NN2(I) = [];
% $$$     I = find(NN2<=0);
% $$$     NN2(I) = [];
% $$$     
% $$$     for b = 1:nboot
% $$$         r = rand(N,1);
% $$$         r = ceil(r*N/1);
% $$$         N2_boot_b(1:length(J),b) = nanmean(NN2(r));
% $$$     end
    
     % ------------------------------ %
    
    
% $$$     N = length(I);
% $$$     for b = 1:nboot
% $$$         r = rand(N,1);
% $$$         r = ceil(r*N/1);
% $$$         N2_boot_b(1:length(J),b) = nanmean(N2_mat(I));
% $$$     end
% $$$     
% $$$     
% $$$    % mean of random sampling
% $$$    N2_boot_dot = nanmean(N2_boot_b);
% $$$ 
% $$$    % compute error (EFRON & GONG 83)
% $$$    N2_error(i) = sqrt(sum((diff([N2_boot_b N2_boot_dot],1, 2)).^2, 2)./(nboot-1));
% $$$ 
% $$$    % Compute 95% confidence interval
% $$$    N2_sort  = sort(N2_boot_b);
% $$$    CI_2p5 = round(2.5/100*nboot);
% $$$    CI_97p5 = round(97.5/100*nboot);
% $$$    N2_2p5(i) = N2_sort(CI_2p5);
% $$$    N2_97p5(i) = N2_sort(CI_97p5);
% $$$    N2_ave(i) = N2_boot_dot;    
% $$$ end

% LOESS-type test
N2_meanvec = nanmean(N2_mat(J, :));
I = find(~isnan(N2_meanvec)==1);
smooth_N2 = loess(time2(I), N2_meanvec(I), sort(time2(I)), 0.3, 1);
[time_N2, II] = sort(time2(I));


save N2_M2.mat N2_2p5 N2_97p5 N2_ave N2_error smooth_N2 time_N2 

figure(1)
clf
imagesc(reg_tide, hab, N2_tide)
set(gca, 'ydir', 'normal') 
ylim([0 50])
c = colorbar;
caxis([-5 -3])
xlabel('time relative to high tide')
ylabel('hab')
ti = ylabel(c,'log(N^2) (s^{-2})', 'FontSize', 8);

figure(2)
clf
errorbar(reg_tide, N2_ave, N2_error, '.k')
ylabel('N^2 (s^{-2})')
xlabel('time relative to high tide')
set(gca, 'yscale', 'log') 
xlim([-7 7]) 
hold on
plot(time_N2, smooth_N2, 'k', 'linewidth', 2)


keyboard

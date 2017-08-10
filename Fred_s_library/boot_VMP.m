function boot_VMP(file_names, zbin, varargin)

% function boot_VMP(file_names, zbin, varargin)
%
% where: - no_profile = number of profile to consider
%        - zbin = bins for the final plot
%        - varargin = name of the file containing K profile
%
% usage ex: boot_VMP('file_names', 5, 'K_all.dat');
%           boot_VMP('file_names', 5);
%           boot_VMP('file_names_veryall', 5,'K_veryall.dat');
%           boot_VMP('file_names', 5,'K_veryall.dat', 'file_names_veryall');
%           boot_VMP('epsProfiles.list', 1,'K_HLC.dat');  (Nutrient Pumping study)
%
% file_names is a file containing eps_profile*.mat files that we want to consider
% In linux, an easy command to do in folder containing *eps*.mat is:
% 
% "ls -1 *eps_profile*.mat | sed 's/\.mat//' > file_names"
%
% will compute a mean diffusion coefficient over the profiles. This
% function needs that shear_analysis had been run and
% eps_profileXXX.mat are in current folder. If varargin is present,
% the function will create a file containing:
%
% - K_mean: a vector containing the mean K
% - P_K: a vector containing the correcponding depths

% author: F. Cyr - 2010/05/25
%
% MODIFICATIONS:
%  - 2011/02/07: Major changes in input parameters. Now uses a list
%  instead a number of profiles.
%  - 2011/02/18: now uses varargin as last argument
%  - 2011/02/21: consider 95% bootstrap CI, instead of standard
%  error
%  - 2011/11/08: Remove the option to choose to skip or not bins in
%  profiles; now the script do bin averages one at a time in the
%  loop on profiles.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with varargin
if isempty(varargin)==1
    sav = 0;
    xtra_contour = 0;
elseif size(varargin,2)==1
    sav = 1;
    xtra_contour = 0;
    outfile = varargin{1};
elseif size(varargin,2)==2
    sav = 1;
    outfile = varargin{1};
    xtra_contour = 1;
    file_names2 = varargin{2};
else
    disp('Wrong input... try "help mean_K"')
    return
end

% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 350;
nboot = 500;    
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';

% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});


siz = size(epsfiles);
no_profile = siz(1); %number of eps_files 

% raw matrix to fill
mat_eps_bin = nan(length(P_bin), no_profile);
mat_K_bin = mat_eps_bin;
mat_N2_bin = mat_eps_bin;

N2_count = 0;
N2_min = -6;

%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%

for profile = 1:no_profile

    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
    maxp(profile) = max(p_eps1(end), p_eps2(end));
    %%%%%%%%%%%%%%%%%%%%%%
    % - EPSILON, N2, K - %
    %%%%%%%%%%%%%%%%%%%%%%
    
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
    
    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8);

    % uses home despike
    EPS = Y;
    
    % ---- mean N2 ---- %
    N2=nan(length(p_k),1);
    if length(p_N)==length(p_k);
        N2(:) = N.^2;
    else
        N2(ind:length(p_N))=N.^2;
    end
% $$$     
% $$$    N2 = 8e-8.*sqrt(N2).^.76;
    
    % ---- compute K_rho ---- %
    K_rho = GAMMA.*EPS./(N2);

   
    % ---- Remove unrealistic diffusivity ---- %
    I = find(log10(N2)<N2_min);
    if ~isempty(I)==1
        K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
        K_rho(I)=NaN;
        N2_count = N2_count+length(I);
    end
    
  
    % ---------- bin profile -------------- %
    for i = 1:length(P_bin)
        I = find(p_k >= P_bin(i)-zbin/2 & p_k <= P_bin(i)+zbin/2);
        mat_eps_bin(i, profile) = nanmean(EPS(I));
        mat_K_bin(i, profile) = nanmean(K_rho(I));
        mat_N2_bin(i, profile) = nanmean(N2(I));
    end
    % ------------------------------------ %
    
    % ------ remove last bins ------------- %
    I = find(~isnan(mat_eps_bin(:,profile))==1);
    
    % keyboard
    no_remove = 0;
    if length(I)<no_remove
        mat_eps_bin(:, profile)=NaN;
        mat_K_bin(:, profile)=NaN;
        mat_N2_bin(:, profile)=NaN;
    elseif no_remove ~= 0
        mat_eps_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_K_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_N2_bin(I(end-no_remove+1:end), profile)=NaN;
    end
    % ------------------------------------ %
end 


%%%%%%%%%%%%%%%%%%%%
% --- bootstrap --- %
%%%%%%%%%%%%%%%%%%%%
% just rename
mat_eps  = mat_eps_bin;
mat_K = mat_K_bin;
mat_N2 = mat_N2_bin;

% -- bootstrap -- %
disp('bootstrap...')
N = size(mat_eps,2);

% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    
    eps_boot_b(:,b) = nanmean(mat_eps(:,r),2);
    K_boot_b(:,b) = nanmean(mat_K(:,r),2);
    N2_boot_b(:,b) = nanmean(mat_N2(:,r),2);
end

% mean of random sampling
eps_boot_dot = nanmean(eps_boot_b,2);
K_boot_dot = nanmean(K_boot_b,2);
N2_boot_dot = nanmean(N2_boot_b,2);

% compute error (EFRON & GONG 83)
eps_error = sqrt(sum((diff([eps_boot_b eps_boot_dot],1, 2)).^2, 2)./(nboot-1));
K_error = sqrt(sum((diff([K_boot_b K_boot_dot],1, 2)).^2, 2)./(nboot-1));
N2_error = sqrt(sum((diff([N2_boot_b N2_boot_dot],1, 2)).^2, 2)./(nboot-1));

% Compute 95% confidence interval
eps_sort = sort(eps_boot_b, 2);
K_sort  = sort(K_boot_b, 2);
N2_sort  = sort(N2_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

eps_2p5 = eps_sort(:,CI_2p5);
K_2p5 = K_sort(:,CI_2p5);
N2_2p5 = N2_sort(:,CI_2p5);
eps_97p5 = eps_sort(:,CI_97p5);
K_97p5 = K_sort(:,CI_97p5);
N2_97p5 = N2_sort(:,CI_97p5);


% $$$ % mean profile
% $$$ eps_ave = nanmean(mat_eps,2);
% $$$ K_ave = nanmean(mat_K,2);
% $$$ N2_ave = nanmean(mat_N2,2);
eps_ave = eps_boot_dot;
K_ave = K_boot_dot;
N2_ave = N2_boot_dot;

% $$$ % remove NaN bins if there is (Camil's CFL profiles)
% $$$ I = find(isnan(eps_ave)==1);
% $$$ if isempty(I)==0;
% $$$     eps_ave(I) = [];
% $$$     K_ave(I) = [];
% $$$     N2_ave(I) = [];
% $$$     eps_error(I) = [];
% $$$     K_error(I) = [];
% $$$     N2_error(I) = [];
% $$$     P_bin(I) = [];
% $$$     
% $$$ end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually remove last bin %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ P_bin(end-1:end) = [];
% $$$ eps_ave(end-1:end) = [];
% $$$ eps_error(end-1:end) = [];
% $$$ N2_ave(end-1:end) = [];
% $$$ N2_error(end-1:end) = [];
% $$$ K_ave(end-1:end) = [];
% $$$ K_error(end-1:end) = [];
% $$$ eps_2p5(end-1:end) = [];
% $$$ K_2p5(end-1:end) = [];
% $$$ N2_2p5(end-1:end) = [];
% $$$ eps_97p5(end-1:end) = [];
% $$$ K_97p5(end-1:end) = [];
% $$$ N2_97p5(end-1:end) = [];

% $$$ P_bin(end) = [];
% $$$ eps_ave(end) = [];
% $$$ eps_error(end) = [];
% $$$ N2_ave(end) = [];
% $$$ N2_error(end) = [];
% $$$ K_ave(end) = [];
% $$$ K_error(end) = [];
% $$$ eps_2p5(end) = [];
% $$$ K_2p5(end) = [];
% $$$ N2_2p5(end) = [];
% $$$ eps_97p5(end) = [];
% $$$ K_97p5(end) = [];
% $$$ N2_97p5(end) = [];


%%%%%%%%%%%%%%%%%%%%
% --- plotting --- %
%%%%%%%%%%%%%%%%%%%%




% ---------- binned ----------- %
figure(2)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 16 10])

subplot(1, 3, 1);
semilogx(eps_ave, P_bin, 'k', 'linewidth', 0.25);

% shade area
hold on
x2 = eps_2p5;
x1 = eps_97p5;
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
semilogx(eps_ave, P_bin, 'k', 'linewidth', 0.25);
hold off
set(gca, 'ydir', 'reverse')
axis([1e-9 5e-4 0 300])
set(gca, 'xtick', [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4])
ylabel('P (dbar)', 'FontSize', 9)
xlabel('\epsilon (W kg^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')

% Write letter identification
%text(2e-6, 295, 'a', ...
%         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');

subplot(1, 3, 2)
semilogx(N2_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
x2 = N2_2p5;
x1 = N2_97p5;
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
semilogx(N2_ave, P_bin, 'k', 'linewidth', 0.25)
hold off
set(gca, 'ydir', 'reverse')
axis([4e-5 5e-3 0 300])
set(gca, 'xtick', [1e-4 1e-3])
xlabel('N^2 (s^{-2})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')

% Write letter identification
%text(3e-3, 176, 'b', ...
%         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');


subplot(1, 3, 3)

semilogx(K_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
%semilogx(K_ave*0+5.5e-5, P_bin, '--k', 'linewidth', 0.25)
x2 = K_2p5;
x1 = K_97p5;
% $$$ I = find(x2<0);
% $$$ if isempty(I)==0    % for CAmil's scripts!!
% $$$ % $$$     %keyboard
% $$$ % $$$     x2(I)=[];
% $$$ % $$$     x1(I)=[];
% $$$ % $$$     P_bin(I)=[];
% $$$ % $$$     K_ave(I)=[];
% $$$ x2(I)=1e-10;
% $$$ end
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
semilogx(K_ave, P_bin, 'k', 'linewidth', 0.25)
hold off
set(gca, 'ydir', 'reverse')
axis([5e-6 1e-2 0 300])
set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
xlabel('K (m^2 s^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')

% Write letter identification
%Ylim = [0 180];
%text(5e-4, 176, 'c', ...
%         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');




% save figure
%print('-dpng', '-r300','ENK_boot_5mbin.png')
print('-dpng', '-r300','ENK_boot_5mbinCI.png')
set(gcf, 'renderer', 'painters')
%print('-depsc2', 'ENK_boot_5mbin.eps')
print('-depsc2', 'ENK_boot_5mbinCI.eps')



% $$$ % figure for error on K
% $$$ figure(3)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 6 6])
% $$$ 
% $$$ plot(K_error./K_ave, P_bin)
% $$$ set(gca, 'ydir', 'reverse')
% $$$ ylabel('P (dbars)', 'FontSize', 10)
% $$$ xlabel('K_{error}/K_{ave}')
% $$$ 
% $$$ %print('-dpng', '-r300','K_booterror_5mbin.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'K_booterror_5mbin.eps')


%%%%%%%%%%%%%%%%%%%%%%
% - within the CIL - %
%%%%%%%%%%%%%%%%%%%%%%
% $$$ zmin_cil = 20;
% $$$ zmax_cil = 180;
% $$$ I = find(P_bin>zmin_cil & P_bin<zmax_cil);
% $$$ 
% $$$ disp('\epsilon_{CIL}...')
% $$$ % $$$ x1 = (eps_ave+eps_error);
% $$$ % $$$ x2 = (eps_ave-eps_error);
% $$$ x2 = eps_2p5;
% $$$ x1 = eps_97p5;
% $$$ [mean(eps_ave(I)) mean(x2(I)) mean(x1(I))]
% $$$ 
% $$$ disp('K_{CIL}...')
% $$$ % $$$ x1 = (K_ave+K_error);
% $$$ % $$$ x2 = (K_ave-K_error);
% $$$ x2 = K_2p5;
% $$$ x1 = K_97p5;
% $$$ [mean(K_ave(I)) mean(x2(I)) mean(x1(I))]



%%%%%%%%%%%%%%%%%%%%%%
% - save variables - %
%%%%%%%%%%%%%%%%%%%%%%

if sav == 1;
    dlmwrite(outfile, [P_bin K_ave K_2p5 K_97p5], 'delimiter',' ','precision',6);
    save boot_VMP.mat P_bin K_ave K_2p5 K_97p5 eps_ave eps_2p5  eps_97p5 N2_ave N2_2p5 N2_97p5
end    
disp(sprintf('  %d profiles used', N))



%%%%%%%%%%%%%%%%%%%%%%
% --- Extra plot --- %
%%%%%%%%%%%%%%%%%%%%%%
if  xtra_contour == 1;
    disp('Extra plot for diffusivity')
    K_ave1 = K_ave;
    N2_ave1 = N2_ave;
    eps_ave1 = eps_ave;
    xtra_bootVMP
end   


% $$$ %%%%%%%%%%%%%%%%%%%%%%%
% $$$ % - histogram of N2 - %
% $$$ %%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ figure(3)
% $$$ clf
% $$$ 
% $$$ I = 1:size(mat_N2,1)*size(mat_N2,2);
% $$$ A = mat_N2(I);
% $$$ V =  [-5.5-0.5/2:0.5:-1.5+0.5/2];
% $$$ hist(log10(A), V);  
% $$$ xlabel('log(N2)')
% $$$ 
% $$$ disp(sprintf('%d diffusivity bins ignored because log(N2)<%d', N2_count, N2_min));
% $$$ percent = round(N2_count/length(I)*100);
% $$$ disp(sprintf('(%d percent of the total)', percent))
% $$$ 

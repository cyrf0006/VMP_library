function boot_VMP_Gustavo(file_names, zbin, zmax)

% function boot_VMP(file_names, zbin, zmax)
%
% where: - no_profile = number of profile to consider
%        - zbin = bins for the final plot
%        - varargin = name of the file containing K profile
%
% usage ex: boot_VMP_Gustavo('20110721_eps.list', .25, 150)

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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
          %zmax = 350;
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
mat_T_bin = mat_eps_bin;
mat_S_bin = mat_eps_bin;
mat_F_bin = mat_eps_bin;

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
    
    %%%%%%%%%%%%%%%%%%%%%%
    % -       T, S     - %
    %%%%%%%%%%%%%%%%%%%%%%
    load(fname(5:end))
    MAXP(profile) = max(p_eps1(end), p_eps2(end));
    
    % ---------- bin profile -------------- %
    for i = 1:length(P_bin)
        I = find(p_k >= P_bin(i)-zbin/2 & p_k <= P_bin(i)+zbin/2);
        mat_eps_bin(i, profile) = nanmean(EPS(I));
        mat_K_bin(i, profile) = nanmean(K_rho(I));
        mat_N2_bin(i, profile) = nanmean(N2(I));
        mat_F_bin(i, profile) = nanmean(FLUORO(I));

        I = find(P >= P_bin(i)-zbin/2 & P <= P_bin(i)+zbin/2);
        mat_T_bin(i, profile) = nanmean(SBT(I)); 
        mat_S_bin(i, profile) = nanmean(SBS(I)); 

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
mat_T = mat_T_bin;
mat_S = mat_S_bin;
mat_F = mat_F_bin;

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
    T_boot_b(:,b) = nanmean(mat_T(:,r),2);
    S_boot_b(:,b) = nanmean(mat_S(:,r),2);
    F_boot_b(:,b) = nanmean(mat_F(:,r),2);

end

% mean of random sampling
eps_boot_dot = nanmean(eps_boot_b,2);
K_boot_dot = nanmean(K_boot_b,2);
N2_boot_dot = nanmean(N2_boot_b,2);
T_boot_dot = nanmean(T_boot_b,2);
S_boot_dot = nanmean(S_boot_b,2);
F_boot_dot = nanmean(F_boot_b,2);

% compute error (EFRON & GONG 83)
eps_error = sqrt(sum((diff([eps_boot_b eps_boot_dot],1, 2)).^2, 2)./(nboot-1));
K_error = sqrt(sum((diff([K_boot_b K_boot_dot],1, 2)).^2, 2)./(nboot-1));
N2_error = sqrt(sum((diff([N2_boot_b N2_boot_dot],1, 2)).^2, 2)./(nboot-1));

% Compute 95% confidence interval
eps_sort = sort(eps_boot_b, 2);
K_sort  = sort(K_boot_b, 2);
N2_sort  = sort(N2_boot_b, 2);
T_sort  = sort(T_boot_b, 2);
S_sort  = sort(S_boot_b, 2);
F_sort  = sort(F_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

eps_2p5 = eps_sort(:,CI_2p5);
K_2p5 = K_sort(:,CI_2p5);
N2_2p5 = N2_sort(:,CI_2p5);
T_2p5 = T_sort(:,CI_2p5);
S_2p5 = S_sort(:,CI_2p5);
F_2p5 = F_sort(:,CI_2p5);

eps_97p5 = eps_sort(:,CI_97p5);
K_97p5 = K_sort(:,CI_97p5);
N2_97p5 = N2_sort(:,CI_97p5);
T_97p5 = T_sort(:,CI_97p5);
S_97p5 = S_sort(:,CI_97p5);
F_97p5 = F_sort(:,CI_97p5);


% $$$ % mean profile
eps_ave = eps_boot_dot;
K_ave = K_boot_dot;
N2_ave = N2_boot_dot;
T_ave = T_boot_dot;
S_ave = S_boot_dot;
F_ave = F_boot_dot;



%%%%%%%%%%%%%%%%%%%%
% --- plotting --- %
%%%%%%%%%%%%%%%%%%%%
% ---------- binned ----------- %
figure(2)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 16 10])

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 5; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.08; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.03; % top of figure
bots = 0.15; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %



subplot(1, 5, 1);
plot(T_ave, P_bin, 'k', 'linewidth', 0.25);
% shade area
hold on
x2 = T_2p5;
x1 = T_97p5;
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
plot(T_ave, P_bin, 'k', 'linewidth', 0.25);
hold off
set(gca, 'ydir', 'reverse')
%axis([0 10 0 zmax])
%set(gca, 'xtick', [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4])
ylabel('P (dbar)', 'FontSize', 9)
xlabel('T (^{\circ}C)', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space

subplot(1, 5, 2)
plot(S_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
x2 = S_2p5;
x1 = S_97p5;
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
plot(S_ave, P_bin, 'k', 'linewidth', 0.25)
hold off
set(gca, 'ydir', 'reverse')
%axis([25 34 0 zmax])
%set(gca, 'xtick', [1e-4 1e-3])
xlabel('S', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


subplot(1, 5, 3)
semilogx(N2_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
x2 = N2_2p5;
x1 = N2_97p5;
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
semilogx(N2_ave, P_bin, 'k', 'linewidth', 0.25)
hold off
set(gca, 'ydir', 'reverse')
axis([4e-5 5e-3 0 zmax])
set(gca, 'xtick', [1e-4 1e-3])
xlabel('N^2 (s^{-2})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space



subplot(1, 5, 4)
plot(F_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
x2 = F_2p5;
x1 = F_97p5;
I = find(~isnan(x1)==1 & ~isnan(x1)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P_bin(I); flipud(P_bin(I)); P_bin(I(1))], [.8 .8 .8], 'edgecolor', 'none');
plot(F_ave, P_bin, 'k', 'linewidth', 0.25)
hold off
set(gca, 'ydir', 'reverse')
%axis([11 18 0 zmax])
%set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
xlabel('F (ppb)', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space



subplot(1, 5, 5);
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
axis([1e-10 5e-5 0 zmax])
set(gca, 'xtick', [1e-9 1e-8 1e-7 1e-6])
set(gca, 'yticklabel', [])
xlabel('\epsilon (W kg^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'xminorgrid', 'off')
set(gca, 'tickdir', 'out')
adjust_space


% save figure
print('-dpng',['TSNFE_' file_names(1:8) '_' num2str(zmax) 'm.png'])
set(gcf, 'renderer', 'painters')
print('-depsc2', ['TSNFE_' file_names(1:8) '_' num2str(zmax) 'm.eps'])



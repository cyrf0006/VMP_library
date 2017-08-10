function boot_VMP_bias(file_names, zbin, varargin)

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
%
% author: F. Cyr - 2010/05/25
%
% MODIFICATIONS:
%  - 2011/02/07: Major changes in input parameters. Now uses a list
%  instead a number of profiles.
%  - 2011/02/18: now uses varargin as last argument
%  - 2011/02/21: consider 95% bootstrap CI, instead of standard error
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
    xtra_contour = 0;
    file_names2 = varargin{2};
elseif size(varargin,3)==3
    sav = 1;
    outfile = varargin{1};
    xtra_contour = 1;
    file_names3 = varargin{3};
else
    disp('Wrong input... try "help mean_K"')
    return
end

% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 200;
nboot = 500;    

%%% - NEAP -- %%

% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles_n = char(C{1});

siz = size(epsfiles_n);
no_profile_n = siz(1); %number of eps_files 

% 1st profile (to compute pressure vector)
fname = epsfiles_n(1, :); 
% remove blank (created if files in list has not the same size)
I = find(fname==' ');   
fname(I) = [];
%keyboard
load(fname);

dp = p_eps1(3)-p_eps1(1); % skip one value out of 2
z1 = zmin+dp/2;                          
p_k = [z1:dp:zmax]';

% raw matrix to fill
mat_eps_n = sparse(length(p_k), no_profile_n);
mat_K_n = sparse(length(p_k), no_profile_n);
mat_N2_n = sparse(length(p_k), no_profile_n);

%%% -- SPRING -- %%%
% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names2);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles_s = char(C{1});

siz = size(epsfiles_s);
no_profile_s = siz(1); %number of eps_files 

% raw matrix to fill
mat_eps_s = sparse(length(p_k), no_profile_s);
mat_K_s = sparse(length(p_k), no_profile_s);
mat_N2_s = sparse(length(p_k), no_profile_s);


N2_count = 0;
N2_min = -6;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on profiles (neap)%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for profile = 1:no_profile_n

    fname = epsfiles_n(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)

    %%%%%%%%%%%%%%%%%%%%%%
    % - EPSILON, N2, K - %
    %%%%%%%%%%%%%%%%%%%%%%
    
    % ---- average of 2 epsilon profiles ---- %
    EPS2=nan(length(p_k),2);
    
    
    % store 1st profile (if not NaNs..)
    if nansum(~isnan(eps1))>0 %0 if eps1 has only Nans
        
        IND1 = []; % IND1 is the first good index of p_eps1
        ind1=zmin; % ind1 is the corresponding index in p_k
        while isempty(IND1)==1;
            ind1 = ind1+1; 
            IND1 = find(p_eps1==p_k(ind1));
        end
        
        IND2 = [];
        ind2=zmax; 
        while isempty(IND2)==1;
            ind2 = ind2-1; 
            IND2 = find(p_eps1==p_k(ind2));
        end
        
        dz = p_eps1(2)-p_eps1(1); %dz in p_eps1/p_eps2 
        EPS2(ind1:ind2, 1) = eps1(IND1:round(1/dz):IND2);
    end
    
     
    % store 2nd profile
    if nansum(~isnan(eps2))>0 %0 if eps1 has only Nans

        IND1 = []; % IND11 is the first good index of p_eps1
        ind1=zmin; % ind1 is the corresponding index in p_k
        while isempty(IND1)==1;
            ind1 = ind1+1; 
            IND1 = find(p_eps2==p_k(ind1));
        end
        
        IND2 = [];
        ind2=zmax; 
        while isempty(IND2)==1;
            ind2 = ind2-1; 
            IND2 = find(p_eps2==p_k(ind2));
        end

        dz = p_eps2(2)-p_eps2(1); %dz in p_eps1/p_eps2 
        EPS2(ind1:ind2, 2) = eps2(IND1:round(1/dz):IND2);        
    end
    
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
    
    EPS = nanmean(EPS2,2); %MEAN EPSILON
    EPS= EPS2(:,1);

    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8); 
    % uses home despike
    EPS = Y;
    
    % ---- mean N2 ---- %
    N2=nan(length(p_k),1);
    N2(ind1:ind2) = N(IND1:round(1/dz):IND2).^2;
    
    % ---- compute K_rho ---- %
    K_rho = GAMMA.*EPS./(N2);

   
    % ---- Remove unrealistic diffusivity ---- %
    I = find(log10(N2)<N2_min);
    if ~isempty(I)==1
        K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
        K_rho(I)=NaN;
        N2_count = N2_count+length(I);
    end
    
    % ---- store all profiles in matrix ---- %
    mat_eps_n(:, profile) = EPS; % those are 1-m binned
    mat_K_n(:, profile) = K_rho;
    mat_N2_n(:, profile) = N2;

end 


%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on profiles (spring)%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for profile = 1:no_profile_s

    fname = epsfiles_s(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)

    %%%%%%%%%%%%%%%%%%%%%%
    % - EPSILON, N2, K - %
    %%%%%%%%%%%%%%%%%%%%%%
    
    % ---- average of 2 epsilon profiles ---- %
    EPS2=nan(length(p_k),2);
    
    
    % store 1st profile (if not NaNs..)
    if nansum(~isnan(eps1))>0 %0 if eps1 has only Nans
        
        IND1 = []; % IND1 is the first good index of p_eps1
        ind1=zmin; % ind1 is the corresponding index in p_k
        while isempty(IND1)==1;
            ind1 = ind1+1; 
            IND1 = find(p_eps1==p_k(ind1));
        end
        
        IND2 = [];
        ind2=zmax; 
        while isempty(IND2)==1;
            ind2 = ind2-1; 
            IND2 = find(p_eps1==p_k(ind2));
        end
        
        dz = p_eps1(2)-p_eps1(1); %dz in p_eps1/p_eps2 
        EPS2(ind1:ind2, 1) = eps1(IND1:round(1/dz):IND2);
    end
    
     
    % store 2nd profile
    if nansum(~isnan(eps2))>0 %0 if eps1 has only Nans

        IND1 = []; % IND11 is the first good index of p_eps1
        ind1=zmin; % ind1 is the corresponding index in p_k
        while isempty(IND1)==1;
            ind1 = ind1+1; 
            IND1 = find(p_eps2==p_k(ind1));
        end
        
        IND2 = [];
        ind2=zmax; 
        while isempty(IND2)==1;
            ind2 = ind2-1; 
            IND2 = find(p_eps2==p_k(ind2));
        end

        dz = p_eps2(2)-p_eps2(1); %dz in p_eps1/p_eps2 
        EPS2(ind1:ind2, 2) = eps2(IND1:round(1/dz):IND2);        
    end
    
    
    % "selection" average (uncomment plotting portion to inspect)
    
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
    
    EPS = nanmean(EPS2,2); %MEAN EPSILON
    EPS= EPS2(:,1);
    
    
    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8);
    % uses home despike
    EPS = Y;
    
    % ---- mean N2 ---- %
    N2=nan(length(p_k),1);
    N2(ind1:ind2) = N(IND1:round(1/dz):IND2).^2;
    
    % ---- compute K_rho ---- %
    K_rho = GAMMA.*EPS./(N2);

   
    % ---- Remove unrealistic diffusivity ---- %
    I = find(log10(N2)<N2_min);
    if ~isempty(I)==1
        K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
        K_rho(I)=NaN;
        N2_count = N2_count+length(I);
    end
    
    % ---- store all profiles in matrix ---- %
    mat_eps_s(:, profile) = EPS; % those are 1-m binned
    mat_K_s(:, profile) = K_rho;
    mat_N2_s(:, profile) = N2;

end 




% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % Center / Border comparison %
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % NOTE: this is weak since it doesnt check profile depth!!
% $$$ hb = 10;
% $$$ CIL = [50 150];
% $$$ count = 1;
% $$$ count2 = 1;
% $$$ for i=1:size(mat_K,2)
% $$$     
% $$$     
% $$$     % cut profile
% $$$     I = find(~isnan(mat_K(:,i))==1);
% $$$     Kprof = mat_K(I,i); % single Kprofile
% $$$     zz = p_k(I); % with its depth vector
% $$$                  
% $$$     % Mean K boundary
% $$$     if  max(zz) > CIL(1) & max(zz)<CIL(2)+hb % at border (50-160m)   
% $$$         J = find(zz > max(zz)-hb & zz<CIL(2) & zz>CIL(1)); %last 10m,but between 50-150
% $$$         Kbndry(count: count+length(J)-1) = Kprof(J);
% $$$         count = count+length(J);
% $$$     end
% $$$     
% $$$     % Mean K out of boundary
% $$$     if max(zz) < CIL(1)+hb % profile too short
% $$$         
% $$$     elseif max(zz) > CIL(2)+hb % deep profile away from boundary
% $$$         J = find(zz > CIL(1) & zz < CIL(2)); % z in bndry layer  
% $$$         Kcenter(count2: count2+length(J)-1) = Kprof(J);
% $$$         count2 = count2+length(J);
% $$$     else % profile finishes in the CIL
% $$$         J = find(zz<max(zz)-hb & zz>CIL(1) & zz<CIL(2)); % z in bndry layer  
% $$$         Kcenter(count2: count2+length(J)-1) = Kprof(J);
% $$$         count2 = count2+length(J);
% $$$     end
% $$$ end



%%%%%%%%%%%%%%%%%
% bin variables %
%%%%%%%%%%%%%%%%%
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';

mat_eps_bin_n = nan(length(P_bin), size(mat_eps_n, 2));
mat_K_bin_n = mat_eps_bin_n;
mat_N2_bin_n = mat_eps_bin_n;

for i = 1:length(P_bin)
    I = find(p_k>(zmin+(i-1)*zbin) & p_k<(zmin+(i)*zbin));
    mat_eps_bin_n(i, :) = nanmean(mat_eps_n(I, :), 1);
    mat_K_bin_n(i, :) = nanmean(mat_K_n(I, :), 1);
    mat_N2_bin_n(i, :) = nanmean(mat_N2_n(I, :), 1);
end
    
mat_eps_bin_s = nan(length(P_bin), size(mat_eps_s, 2));
mat_K_bin_s = mat_eps_bin_s;
mat_N2_bin_s = mat_eps_bin_s;

for i = 1:length(P_bin)
    I = find(p_k>(zmin+(i-1)*zbin) & p_k<(zmin+(i)*zbin));
    mat_eps_bin_s(i, :) = nanmean(mat_eps_s(I, :), 1);
    mat_K_bin_s(i, :) = nanmean(mat_K_s(I, :), 1);
    mat_N2_bin_s(i, :) = nanmean(mat_N2_s(I, :), 1);
end


%%%%%%%%%%%%%%%%%%%%
% --- bootstrap --- %
%%%%%%%%%%%%%%%%%%%%
% just rename
mat_eps_n  = mat_eps_bin_n;
mat_K_n = mat_K_bin_n;
mat_N2_n = mat_N2_bin_n;
mat_eps_s  = mat_eps_bin_s;
mat_K_s = mat_K_bin_s;
mat_N2_s = mat_N2_bin_s;

% -- bootstrap -- %
disp('bootstrap...')
Nn = size(mat_eps_n,2);
Ns = size(mat_eps_s,2);

% create random sampling
for b = 1:nboot
    rn = rand(min([Ns Nn]),1);
    rn = ceil(rn*Nn/1);
   
    rs = rand(min([Ns Nn]),1);
    rs = ceil(rs*Ns/1);
    
    
    mat_eps = [mat_eps_n(:,rn) mat_eps_s(:,rs)];
    mat_K = [mat_K_n(:,rn) mat_K_s(:,rs)];
    mat_N2 = [mat_N2_n(:,rn) mat_N2_s(:,rs)];
    %keyboard
    eps_boot_b(:,b) = nanmean(mat_eps,2);
    K_boot_b(:,b) = nanmean(mat_K,2);
    N2_boot_b(:,b) = nanmean(mat_N2,2);
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

% remove NaN bins if there is (Camil's CFL profiles)
I = find(isnan(eps_ave)==1);
if isempty(I)==0;
    eps_ave(I) = [];
    K_ave(I) = [];
    N2_ave(I) = [];
    eps_error(I) = [];
    K_error(I) = [];
    N2_error(I) = [];
    P_bin(I) = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually remove last bin %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_bin(end-1:end) = [];
eps_ave(end-1:end) = [];
eps_error(end-1:end) = [];
N2_ave(end-1:end) = [];
N2_error(end-1:end) = [];
K_ave(end-1:end) = [];
K_error(end-1:end) = [];
eps_2p5(end-1:end) = [];
K_2p5(end-1:end) = [];
N2_2p5(end-1:end) = [];
eps_97p5(end-1:end) = [];
K_97p5(end-1:end) = [];
N2_97p5(end-1:end) = [];

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
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 14 8])

subplot(1, 3, 1);
semilogx(eps_ave, P_bin, 'k', 'linewidth', 0.25);

% shade area
hold on
x2 = eps_2p5;
x1 = eps_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
hold off
set(gca, 'ydir', 'reverse')
axis([1e-9 5e-6 0 180])
set(gca, 'xtick', [1e-9 1e-8 1e-7 1e-6])
set(gca, 'xminortick', 'off')
ylabel('P (dbar)', 'FontSize', 9)
xlabel('\epsilon (W kg^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'xminortick', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

% Write letter identification
text(2e-6, 176, 'a', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');

subplot(1, 3, 2)
semilogx(N2_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
x2 = N2_2p5;
x1 = N2_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
hold off
set(gca, 'ydir', 'reverse')
axis([5e-5 5e-3 0 180])
set(gca, 'xtick', [1e-4 1e-3])
xlabel('N^2 (s^{-2})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

% Write letter identification
text(3e-3, 176, 'b', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');


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
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
set(gca, 'ydir', 'reverse')
axis([5e-6 1e-3 0 180])
set(gca, 'xtick', [1e-5 1e-4 1e-3])
set(gca, 'xminortick', 'off')
xlabel('K (m^2 s^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')

% Write letter identification
Ylim = [0 180];
text(5e-4, 176, 'c', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%
% --- Adjust space --- %
%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ subplot(1, 3, 1);
% $$$ pos1 = get(gca, 'pos');
% $$$ subplot(1, 3, 2);
% $$$ pos2 = get(gca, 'pos');
% $$$ subplot(1, 3, 3);
% $$$ pos3 = get(gca, 'pos');
% $$$ 
% $$$ dx = 0.05; % blank between figures
% $$$ olddx = pos2(1)-(pos1(1)+pos1(2));
% $$$ oldwidth = pos1(2);
% $$$ newwidth = oldwidth+(olddx/2-dx/2);
% $$$ pos1(2) = newwidth;
% $$$ pos2(1) = pos1(1)+pos1(2)+dx; 
% $$$ pos2(2) = newwidth;
% $$$ pos3(1) = pos2(1)+pos2(2)+dx; 
% $$$ pos3(2) = newwidth;
% $$$ 
% $$$ subplot(1, 3, 1);
% $$$ set(gca, 'pos', pos1)
% $$$ subplot(1, 3, 2);
% $$$ set(gca, 'pos', pos2)
% $$$ subplot(1, 3, 3);
% $$$ set(gca, 'pos', pos3)


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
zmin_cil = 20;
zmax_cil = 180;
I = find(P_bin>zmin_cil & P_bin<zmax_cil);

disp('\epsilon_{CIL}...')
% $$$ x1 = (eps_ave+eps_error);
% $$$ x2 = (eps_ave-eps_error);
x2 = eps_2p5;
x1 = eps_97p5;
[mean(eps_ave(I)) mean(x2(I)) mean(x1(I))]

disp('K_{CIL}...')
% $$$ x1 = (K_ave+K_error);
% $$$ x2 = (K_ave-K_error);
x2 = K_2p5;
x1 = K_97p5;
[mean(K_ave(I)) mean(x2(I)) mean(x1(I))]



%%%%%%%%%%%%%%%%%%%%%%
% - save variables - %
%%%%%%%%%%%%%%%%%%%%%%

if sav == 1;
    dlmwrite(outfile, [P_bin K_ave K_2p5 K_97p5], 'delimiter',' ','precision',6);
end    

disp(sprintf('  %d neap profiles', Nn))
disp(sprintf('  %d spring profiles', Ns))
disp(sprintf('  dataset subsampled with N=%d', min([Ns Nn])))




%%%%%%%%%%%%%%%%%%%%%%
% --- Extra plot --- %
%%%%%%%%%%%%%%%%%%%%%%
%file_names2 = file_names3;
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

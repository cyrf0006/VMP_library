function mean_K(file_names, zbin, varargin)

% function mean_K(file_names, zbin, varargin)
%
% where: - no_profile = number of profile to consider
%        - zbin = bins for the final plot
%        - varargin = 1) name of the file containing K profile 
%                     2)file_names2, the extra file names on plot
%                     (MUST be 2nd varargin!!!)
%
% usage ex: mean_K('file_names', 5, 'K_all.dat');
%           mean_K('file_names', 5);
%           mean_K('file_names', 5, 'K_all.dat', 'file_names_veryall');
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
%  - 2011/04/18: Extra contour boundary/center (as varargin in)
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
zmax = 200;
    

% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

siz = size(epsfiles);
no_profile = siz(1); %number of eps_files 

% 1st profile (to compute pressure vector)
fname = epsfiles(1, :); 
load(fname);

dp = p_eps1(3)-p_eps1(1); % skip one value out of 2
z1 = zmin+dp/2;                          
p_k = [z1:dp:zmax]';

% raw matrix to fill
mat_eps = sparse(length(p_k), no_profile);
mat_K = sparse(length(p_k), no_profile);
mat_N2 = sparse(length(p_k), no_profile);


%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%
for profile = 1:no_profile

    fname = epsfiles(profile, :);     
    load(fname)

    %%%%%%%%%%%%%%%%%%%%%%
    % - EPSILON, N2, K - %
    %%%%%%%%%%%%%%%%%%%%%%
    
    % ---- average of 2 epsilon profiles ---- %
    EPS2=nan(length(p_k),2);
    
    
    % ranging 1st profile (if not NaNs..)
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
    
     
    % ranging 2nd profile
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
    
    
    % Selective average (uncomment plotting portion to inspect)

% $$$     plot(EPS2(:,1), p_k)
% $$$     hold on
% $$$     plot(EPS2(:,2), p_k, 'r')
    
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
    
    EPS = nanmean(EPS2,2); %MEAN EPSILON
    
% $$$     plot(EPS, p_k, 'y')
% $$$     set(gca, 'ydir', 'reverse')
% $$$     title(sprintf('profile %d', profile))
% $$$     hold off
% $$$     pause    
    
    
    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8);
    No; % check plotting section
% $$$  plot(EPS, p_k);
% $$$  hold on
% $$$  plot(Y, p_k, 'r');
% $$$  title(sprintf('profile %d', profile))
% $$$  hold off
% $$$  pause
 
    % uses home despike
    EPS = Y;
    
    % ---- mean N2 ---- %
    N2=nan(length(p_k),1);
    N2(ind1:ind2) = N(IND1:round(1/dz):IND2).^2;
   
    % ---- compute K_rho ---- %
    K_rho = GAMMA.*EPS./(N2);

    % ---- store all profiles in matrix ---- %
    mat_eps(:, profile) = EPS;
    mat_K(:, profile) = K_rho;
    mat_N2(:, profile) = N2;

end 

%%%%%%%%%%%%%%%%%%%%
% loop on each bin %
%%%%%%%%%%%%%%%%%%%%
for bin = 1:length(p_k)
    
    % ---- for epsilon ---- %
    V = mat_eps(bin, :); %vector containing value for this bin from every profile
    I = find(isnan(V)==0);
    V = V(I); % remove NaN
    num = length(V); %number of good data for this bin (number of profile with values is this bin

    if num>round(.1*no_profile); % at least 10% of the number of profile has value for this bin
        V = sort(V);
        eps_5p(bin) = V(ceil(.05*num));
        eps_95p(bin) = V(floor(.95*num));
        %eps_mean(bin) = nanmean(V(ceil(.05*num):floor(.95*num))); %means 5-95%
        eps_mean(bin) = nanmean(V); %real mean           
            
    else
        eps_5p(bin) = NaN;
        eps_95p(bin) = NaN;
        eps_mean(bin) = NaN;
    end

    % ---- for N2 ---- %
    V = mat_N2(bin, :); %vector containing value for this bin from every profile
    I = find(isnan(V)==0);
    V = V(I); % remove NaN
    num = length(V); %number of good data for this bin (number of profile with values is this bin
    if num>round(.1*no_profile); % at least 10% of the number of profile has value for this bin
        V = sort(V);
        N2_5p(bin) = V(ceil(.05*num));
        N2_95p(bin) = V(floor(.95*num));
        %N2_mean(bin) = nanmean(V(ceil(.05*num):floor(.95*num)));
        N2_mean(bin) = nanmean(V); %real mean           
    else
        N2_5p(bin) = NaN;
        N2_95p(bin) = NaN;
        N2_mean(bin) = NaN;
    end   

    % ---- for K ---- %
    V = mat_K(bin, :); %vector containing value for this bin from every profile
    I = find(isnan(V)==0);
    V = V(I); % remove NaN
    num = length(V); %number of good data for this bin (number of profile with values is this bin

    if num>round(.1*no_profile); % at least 10% of the number of profile has value for this bin
        V = sort(V);
        K_5p(bin) = V(ceil(.05*num));
        K_95p(bin) = V(floor(.95*num));
        %K_mean(bin) = nanmean(V(ceil(.05*num):floor(.95*num)));
        K_mean(bin) = nanmean(V); %real mean           
    else
        K_5p(bin) = NaN;
        K_95p(bin) = NaN;
        K_mean(bin) = NaN;
    end          


end


%%%%%%%%%%%%%%%%%%%%
% binned variables %
%%%%%%%%%%%%%%%%%%%%
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';

eps_mean_bin = P_bin*0;
eps_5p_bin = eps_mean_bin;
eps_95p_bin = eps_mean_bin;

K_mean_bin=eps_mean_bin;
K_5p_bin = eps_mean_bin;
K_95p_bin = eps_mean_bin;

N2_mean_bin=eps_mean_bin;
N2_5p_bin = eps_mean_bin;
N2_95p_bin = eps_mean_bin;

for i = 1:length(P_bin)
    I = find(p_k>(zmin+(i-1)*zbin) & p_k<(zmin+(i)*zbin));
    eps_mean_bin(i) = nanmean(eps_mean(I));
    eps_5p_bin(i) = nanmean(eps_5p(I));
    eps_95p_bin(i) = nanmean(eps_95p(I));

    N2_mean_bin(i)=nanmean(N2_mean(I));
    N2_5p_bin(i) = nanmean(N2_5p(I));
    N2_95p_bin(i) = nanmean(N2_95p(I));
    
    K_mean_bin(i)=nanmean(K_mean(I));
    K_5p_bin(i) = nanmean(K_5p(I));
    K_95p_bin(i) = nanmean(K_95p(I));
end


%%%%%%%%%%%%%%
% remove NaN %
%%%%%%%%%%%%%%
        
J = find(isnan(eps_5p)==0);
eps_5p = eps_5p(J);
eps_95p = eps_95p(J);
eps_mean = eps_mean(J);
P_eps = p_k(J)';

J = find(isnan(N2_5p)==0);
N2_5p = N2_5p(J);
N2_95p = N2_95p(J);
N2_mean = N2_mean(J);
P_N2 = p_k(J)';

J = find(isnan(K_5p)==0);
K_5p = K_5p(J);
K_95p = K_95p(J);
K_mean = K_mean(J);  
P_K = p_k(J)';

J = find(isnan(eps_5p_bin)==0);
eps_5p_bin = eps_5p_bin(J);
eps_95p_bin = eps_95p_bin(J);
eps_mean_bin = eps_mean_bin(J);
P_eps_bin = P_bin(J);

J = find(isnan(N2_5p_bin)==0);
N2_5p_bin = N2_5p_bin(J);
N2_95p_bin = N2_95p_bin(J);
N2_mean_bin = N2_mean_bin(J);
P_N2_bin = P_bin(J);

J = find(isnan(K_5p_bin)==0);
K_5p_bin = K_5p_bin(J);
K_95p_bin = K_95p_bin(J);
K_mean_bin = K_mean_bin(J);  
P_K_bin = P_bin(J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually remove last bin %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_5p_bin(end) = [];
eps_95p_bin(end) = [];
eps_mean_bin(end) = [];
P_eps_bin(end) = [];

N2_5p_bin(end) = [];
N2_95p_bin(end) = [];
N2_mean_bin(end) = [];
P_N2_bin(end) = [];

K_5p_bin(end) = [];
K_95p_bin(end) = [];
K_mean_bin(end) = [];
P_K_bin(end) = [];

% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % Epsilon distribution in the CIL %
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
idep = find(P_K_bin>=50 & P_eps_bin<=150);

disp('\epsilon CIL mean: '); mean(eps_mean_bin(idep)) 
disp('5% centile: '); mean(eps_5p_bin(idep)) 
disp('95% centile: '); mean(eps_95p_bin(idep))
disp(' ')

disp('K CIL mean: '); mean(K_mean_bin(idep)) 
disp('5% centile: '); mean(K_5p_bin(idep)) 
disp('95% centile: '); mean(K_95p_bin(idep))


%%%%%%%%%%%%%%%%%%%
% customized plot %
%%%%%%%%%%%%%%%%%%%

% ------ Not binned ------ %

figure(1)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 16 12])

subplot(1, 3, 1);
semilogx(eps_5p, P_eps)
set(gca, 'ydir', 'reverse')
axis([1e-10 1e-5 0 180])
set(gca, 'xtick', [1e-10 1e-8 1e-6])
%set(gca, 'xticklabel', [1e-10 1e-8 1e-6])
hold on
semilogx(eps_95p, P_eps, 'r')
semilogx(eps_mean, P_eps, 'k')
hold off
ylabel('Depth (m)', 'FontSize', 18)
%xlabel('\epsilon (W/kg)', 'FontSize', 18)
title('\epsilon (W/kg)', 'FontSize', 18)
set(gca, 'FontSize', 18)

subplot(1, 3, 2)
plot(N2_5p, P_N2)
set(gca, 'ydir', 'reverse')
axis([0 0.005 0 180])
hold on
plot(N2_95p, P_N2, 'r')
plot(N2_mean, P_N2, 'k')
hold off
%xlabel('N^2 (s^{-2})', 'FontSize', 18)
title('N^2 (s^{-2})', 'FontSize', 18)
set(gca, 'FontSize', 18)
set(gca, 'yticklabel', [])


subplot(1, 3, 3)
semilogx(K_5p, P_K)
%plot(K_5p, P_K)
set(gca, 'ydir', 'reverse')
axis([1e-7 1e-3 0 180])
set(gca, 'xtick', [1e-7 1e-5 1e-3])
%set(gca, 'xticklabel', [10^(-8) 10^(-6) 10^(-4)])
hold on
%plot(K_95p, P_K, 'r')
%plot(K_mean, P_K, 'k')
semilogx(K_95p, P_K, 'r')
semilogx(K_mean, P_K, 'k')
hold off
%xlabel('K_{\rho} (m^2/s)', 'FontSize', 10)
title('K_{\rho} (m^2/s)', 'FontSize', 10)
set(gca, 'FontSize', 10)
set(gca, 'yticklabel', [])

subplot(1, 3, 2)
%L = ['5%';'95%'; 'mean(5%-95%)'];
legend('5%', '95%', 'mean(5-95%)', 'Location','SouthEast')
%legend('5%','95%', 'mean(5%-95%)')

print('-depsc2', 'ENK_1m-bin.eps')
 

% ---------- binned ----------- %
figure(2)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 14 15])

subplot(1, 3, 1);
semilogx(eps_mean_bin, P_eps_bin)
set(gca, 'ydir', 'reverse')
axis([5e-11 5e-6 0 180])
set(gca, 'xtick', [1e-10 1e-8 1e-6])
%set(gca, 'xticklabel', [1e-10 1e-8 1e-6])
% shade area
hold on
x1 = (eps_5p_bin);
x2 = (eps_95p_bin);
patch([x1; flipud(x2); x1(1)], [P_eps_bin; flipud(P_eps_bin); P_eps_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
semilogx(eps_mean_bin, P_eps_bin, 'k', 'linewidth', 1);
hold off
ylabel('depth (m)', 'FontSize', 10)
xlabel('\epsilon (W kg^{-1})', 'FontSize', 10)
%title('\epsilon (W/kg)', 'FontSize', 10)
set(gca, 'FontSize', 10)
%set(gca, 'xgrid', 'on')


subplot(1, 3, 2)
semilogx(N2_mean_bin, P_N2_bin)
set(gca, 'ydir', 'reverse')
axis([1e-5 1e-2 0 180])
set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2])
% shade area
hold on
x1 = N2_5p_bin;
x2 = N2_95p_bin;
patch([x1; flipud(x2); x1(1)], [P_N2_bin; flipud(P_N2_bin); P_N2_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
semilogx(N2_mean_bin, P_N2_bin, 'k', 'linewidth', 1)
hold off
xlabel('N^2 (s^{-2})', 'FontSize', 10)
%title('N^2 (s^{-2})', 'FontSize', 10)
set(gca, 'FontSize', 10)
set(gca, 'yticklabel', [])
%set(gca, 'xgrid', 'on')
set(gca, 'xminortick', 'off')



subplot(1, 3, 3)
semilogx(K_mean_bin, P_K_bin)
set(gca, 'ydir', 'reverse')
axis([1e-7 1e-3 0 180])
set(gca, 'xtick', [1e-7 1e-5 1e-3])
%set(gca, 'xticklabel', [10^(-8) 10^(-6) 10^(-4)])
% shade area
hold on
x1 = (K_5p_bin);
x2 = (K_95p_bin);
patch([x1; flipud(x2); x1(1)], [P_K_bin; flipud(P_K_bin); P_K_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
semilogx(K_mean_bin, P_K_bin, 'k', 'linewidth', 1)
%xtra plot for cst K=5.5e-5
semilogx(K_mean_bin*0+5.5e-5, P_K_bin, '--k', 'linewidth', 0.5)
hold off
xlabel('K (m^2 s^{-1})', 'FontSize', 10)
%title('K (m^2 s^{-1})', 'FontSize', 10)
set(gca, 'FontSize', 10)
set(gca, 'yticklabel', [])
%set(gca, 'xgrid', 'on')
                                                                                                                                                                                                            
%subplot(1, 3, 2)
%L = ['5%';'95%'; 'mean(5%-95%)'];
%legend('5%', '95%', 'mean(5-95%)', 'Location','SouthEast')
%legend('5%','95%', 'mean(5%-95%)')

print('-depsc2', 'ENK_5m-bin.eps')
print('-dpng', '-r300','ENK_5m-bin.png')

%%%%%%%%%%%%%%%%%%%%%%
% - save variables - %
%%%%%%%%%%%%%%%%%%%%%%

if sav == 1;
    dlmwrite(outfile, [P_K_bin K_mean_bin], 'delimiter',' ','precision',6);
end    


%%%%%%%%%%%%%%%%%%%%%%
% --- Extra plot --- %
%%%%%%%%%%%%%%%%%%%%%%
keyboard
if  xtra_contour == 1;
    disp('Extra plot for diffusivity')
    mean_xtraK
end   

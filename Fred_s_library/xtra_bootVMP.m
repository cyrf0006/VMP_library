% load *.P files names (file in which are recorded .P files)

fid = fopen(file_names2);
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


clear eps_boot_b K_boot_b N2_boot_b
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

% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % manually remove last bin %
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- plot ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_gray = [.6 .6 .6];


figure(2)    
subplot(1, 3, 1);
hold on
x2 = eps_2p5;
x1 = eps_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], my_gray, 'edgecolor', 'none');
semilogx(eps_ave, P_bin, 'k', 'linewidth', 0.25);
semilogx(eps_ave1, P_bin, '--k', 'linewidth', 0.25);
%plot([1e-10 1e-3], [20 20], '--k')
hold off
set(gca, 'box', 'on')

subplot(1, 3, 2)
hold on
x2 = N2_2p5;
x1 = N2_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], my_gray, 'edgecolor', 'none');
semilogx(N2_ave, P_bin, 'k', 'linewidth', 0.25)
semilogx(N2_ave1, P_bin, '--k', 'linewidth', 0.25)
hold off
set(gca, 'box', 'on')

subplot(1, 3, 3)
hold on
x2 = K_2p5;
x1 = K_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], my_gray, 'edgecolor', 'none');
semilogx(K_ave, P_bin, 'k', 'linewidth', 0.25)
semilogx(K_ave1, P_bin, '--k', 'linewidth', 0.25)
hold off
set(gca, 'box', 'on')


%%%%%%%%%%%%%%%%%%%%%%%%
% --- Adjust space --- %
%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1, 3, 1);
pos1 = get(gca, 'pos');
subplot(1, 3, 2);
pos2 = get(gca, 'pos');
subplot(1, 3, 3);
pos3 = get(gca, 'pos');

dx = 0.05; % blank between figures
olddx = pos2(1)-(pos1(1)+pos1(2));
oldwidth = pos1(2);
newwidth = oldwidth+(olddx/2-dx/2);
pos1(2) = newwidth;
pos2(1) = pos1(1)+pos1(2)+dx; 
pos2(2) = newwidth;
pos3(1) = pos2(1)+pos2(2)+dx; 
pos3(2) = newwidth;

subplot(1, 3, 1);
set(gca, 'pos', pos1)
subplot(1, 3, 2);
set(gca, 'pos', pos2)
subplot(1, 3, 3);
set(gca, 'pos', pos3)


print('-depsc2', 'ENK_5m-bin_2plots.eps')
print('-dpng', '-r300','ENK_5m-bin_2plots.png')





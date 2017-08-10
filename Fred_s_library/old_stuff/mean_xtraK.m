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

N2_count = 0;
N2_min = -6;
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

% $$$     plot(EPS2(:,1), p_k)
% $$$     hold on
% $$$     plot(EPS2(:,2), p_k, 'r')
    
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
    
    EPS = nanmean(EPS2,2); %MEAN EPSILON
    EPS= EPS2(:,1);
% $$$     plot(EPS, p_k, 'y')
% $$$     set(gca, 'ydir', 'reverse')
% $$$     title(sprintf('profile %d', profile))
% $$$     hold off
% $$$     pause    
    
    
    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8);
    % check plotting section
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

   
    % ---- Remove unrealistic diffusivity ---- %
    I = find(log10(N2)<N2_min);
    if ~isempty(I)==1
        K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
        K_rho(I)=NaN;
        N2_count = N2_count+length(I);
    end
    
    % ---- store all profiles in matrix ---- %
    mat_eps(:, profile) = EPS; % those are 1-m binned
    mat_K(:, profile) = K_rho;
    mat_N2(:, profile) = N2;

end 


%%%%%%%%%%%%%%%%%
% bin variables %
%%%%%%%%%%%%%%%%%
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';

mat_eps_bin = nan(length(P_bin), size(mat_eps, 2));
mat_K_bin = mat_eps_bin;
mat_N2_bin = mat_eps_bin;

for i = 1:length(P_bin)
    I = find(p_k>(zmin+(i-1)*zbin) & p_k<(zmin+(i)*zbin));
    mat_eps_bin(i, :) = nanmean(mat_eps(I, :), 1);
    mat_K_bin(i, :) = nanmean(mat_K(I, :), 1);
    mat_N2_bin(i, :) = nanmean(mat_N2(I, :), 1);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- plot ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_gray = [.8 .8 .8];


figure(2)    
subplot(1, 3, 1);
hold on
x2 = eps_2p5;
x1 = eps_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], my_gray, 'edgecolor', 'none');
%semilogx(eps_ave, P_bin, 'k', 'linewidth', 1);
hold off

subplot(1, 3, 2)
hold on
x2 = N2_2p5;
x1 = N2_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
%semilogx(N2_ave, P_bin, 'k', 'linewidth', 1)
hold off

subplot(1, 3, 3)
hold on
x2 = K_2p5;
x1 = K_97p5;
patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
%semilogx(K_ave, P_bin, 'k', 'linewidth', 1)


print('-depsc2', 'ENK_5m-bin_2plots.eps')
print('-dpng', '-r300','ENK_5m-bin_2plots.png')

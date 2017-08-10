function eps_tide(file_names, time2, zbin, ranges)

% function boot_VMP(file_names, time2, zbin, ranges)
%
% where: - 
%        - zbin = bins for the final plot
%        - 
%
% usage ex:  eps_tide('file_names_riki', 'time_to_hightide_riki.mat', 1, [20 180]) 
% usage ex:  eps_tide('file_names_bndry_neap', 'time_to_hightide_bndry_neap.mat', 1, [20 180]) 
%
% file_names is a file containing eps_profile*.mat files that we want to consider
% In linux, an easy command to do in folder containing *eps*.mat is:
% 
% "ls -1 *eps_profile*.mat | sed 's/\.mat//' > file_names"
%
% author: F. Cyr - 2011/11/08
%
% MODIFICATIONS:
%  - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 180;
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

% ----------- loop on profiles ------------- % 
for profile = 1:no_profile

    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)

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
% ---------------------------------------------------- % 

% ------------------ compute time relative to high tide  --------------------- % 
load(time2);
time2 = A;

I = find(P_bin > ranges(1) & P_bin < ranges(2));
EPS = nanmean(mat_eps_bin(I, :), 1);
N2 = nanmean(mat_N2_bin(I, :), 1);

% Save raw epsilon vs tide
II = find(isnan(EPS)==1);
EPS(II)=[]; time2(II)=[];
save eps_tide_raw EPS time2


% regular tide vector (-6, 6)
reg_tide = unique(round(time2));

reg_eps = nan(round(length(EPS)/2), length(reg_tide)); % /2 its just approx

for i = 1: length(EPS)
    [Y, I] = min(abs(time2(i)-reg_tide));
    J = find(isnan(reg_eps(:, I))==1);
    reg_eps(J(1), I)=EPS(i);
    reg_N2(J(1), I)=N2(i);
end

% ---------------------------------------------------------------------------- % 


% ------------------ Bootstrap --------------------- %
for i = 1:length(reg_tide)
    
    I = find(~isnan(reg_eps(:,i)));
    samp = reg_eps(I,i);
    N = length(samp);
    clear eps_boot_b N2_boot_b
    % create random sampling
    for b = 1:nboot
        r = rand(N,1);
        r = ceil(r*N/1);
        eps_boot_b(b) = nanmean(samp(r));
        N2_boot_b(b) = nanmean(samp(r));
    end

    % mean of random sampling
    eps_boot_dot(i) = nanmean(eps_boot_b);
    N2_boot_dot(i) = nanmean(N2_boot_b);

    
    % compute error (EFRON & GONG 83)
    eps_error(i) = sqrt(sum((diff([eps_boot_b eps_boot_dot(i)],1, 2)).^2, 2)./(nboot-1));
    N2_error(i) = sqrt(sum((diff([N2_boot_b N2_boot_dot(i)],1, 2)).^2, 2)./(nboot-1));

end


errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
%set(gca, 'ytick', [])
set(gca, 'yscale', 'log')

xlim([-6.5 6.5])
ylabel('\epsilon (W kg^{-1})')
xlabel('time versus high tide (hour)')
keyboard

save eps_tide_bndry_spring reg_tide eps_boot_dot eps_error

save N2_tide_riki reg_tide N2_boot_dot N2_error
% ---------------------------------------------------- %

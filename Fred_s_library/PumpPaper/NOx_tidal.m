clear all
close all

% some info
epsFileNames = '~/PhD/WINDEX/MissionTadoussac/epsProfiles.list';
pFileNames = '~/PhD/WINDEX/MissionTadoussac/pProfiles.list';
zbin = 1;

% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 350;
nboot = 1000;    
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';


% load list of files
fid = fopen(epsFileNames);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

fid = fopen(pFileNames);
C = textscan(fid, '%s', 'delimiter', '\n');
pfiles = char(C{1});


siz = size(epsfiles);
no_profile = siz(1); %number of eps_files 

% raw matrix to fill
mat_eps_bin = nan(length(P_bin), no_profile);
mat_K_bin = mat_eps_bin;
mat_N2_bin = mat_eps_bin;
mat_NO3_bin = mat_eps_bin;
mat_SA_bin = mat_eps_bin;
mat_CT_bin = mat_eps_bin;
mat_RHO_bin = mat_eps_bin;
timeVec = nan(1, no_profile);


N2_count = 0;
N2_min = -6;

%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%

for profile = 1:no_profile
   
    % load turbulence
    disp(sprintf('profile %d / %d', profile, no_profile));
    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    load(fname);
    
    % load fine scale
    fname = pfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    load(fname);
      
    timeVec(profile) = mtime(1); 
    % ---- Process SA and NO3
    if ~exist('SA')
        [SA, in_ocean] = gsw_SA_from_SP(SBS,P,-68.79,48.61);
        CT = gsw_CT_from_t(SA,SBT,P);
    end
    
    % ---- Process EPSILON, N2, K ---- %
    maxp(profile) = max(p_eps1(end), p_eps2(end));
    
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
        clear SA CT
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
        I = find(P >= P_bin(i)-zbin/2 & P <= P_bin(i)+zbin/2);
        %mat_NO3_bin(i, profile) = nanmean(NO3(I));
        mat_SA_bin(i, profile) = nanmean(SA(I));
        %mat_CT_bin(i, profile) = nanmean(CT(I));
        %  mat_RHO_bin(i, profile) = nanmean(rho(I));
    end
    % ------------------------------------ %
    
    % ------ remove last bins ------------- %
    I = find(~isnan(mat_eps_bin(:,profile))==1);
   
    no_remove = 0;
    if length(I)<no_remove
        mat_eps_bin(:, profile)=NaN;
        mat_K_bin(:, profile)=NaN;
        mat_N2_bin(:, profile)=NaN;
        mat_NO3_bin(:, profile)=NaN;
    elseif no_remove ~= 0
        mat_eps_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_K_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_N2_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_NO3_bin(I(end-no_remove+1:end), profile)=NaN;
    end
    % ------------------------------------ %
    clear SA CT
end 



%% HERE THE REAL CALCULATION 
% sort salinity
mat_SA_bin = sort(mat_SA_bin,1);


mat_NO3_bin = 7.3.*mat_SA_bin - 227; %32-34.3 SA (no winter)
I = find(mat_SA_bin<32 | mat_SA_bin>34.3);
mat_NO3_bin(I) = NaN;

% ---- NO3 fluxes calculation ---- %
[dNO3dx, dNO3dz] = gradient(mat_NO3_bin);
dNO3dz = dNO3dz./zbin;

FNO3 = mat_K_bin.*dNO3dz;


%% Check tidal modulation
time2 = time2hightide(timeVec, ['/home/cyrf0006/WINDEX/' ...
                    'data_processing/BMix_study/tide_2009-2012.dat']);

I  = find(P_bin>=25 & P_bin<=50);
FVec = nanmean(FNO3(I,:));
EVec = nanmean(mat_eps_bin(I,:));
KVec = nanmean(mat_K_bin(I,:));
dNVec = nanmean(dNO3dz(I,:));


%% now plot.
NOx_tidal_plot
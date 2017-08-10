function ENK_CIL_bottom(profile_names, epsprofile_names, hb)

% function ENK_CIL_bottom(profile_names, epsprofile_names, hb)
%
% To compute mean dissipation and diffusivity within the CIL and
% within the bottom boundary layer (hb). This script is modified
% from ENK_CIL.m.
%    
%
% usage ex: ENK_CIL_bottom('hit_bottom_profiles', 'hit_bottom_eps_profiles', 10)
%
% file_names is a file containing eps_profile*.mat files that we want to consider
% In linux, an easy command to do in folder containing *eps*.mat is:
% 
% "ls -1 *eps_profile*.mat | sed 's/\.mat//' > epsprofile_names"
%
% author: F. Cyr - 2011/05/26
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 200;
nboot = 500;    

% load file names
fid = fopen(profile_names);
C = textscan(fid, '%s', 'delimiter', '\n');
pfiles = char(C{1});

fid = fopen(epsprofile_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

siz = size(epsfiles);
no_profile = siz(1); %number of eps_files 

% 1st profile (to compute pressure vector)
fname = epsfiles(1, :); 
% remove blank (created if files in list has not the same size)
I = find(fname==' ');   
fname(I) = [];
load(fname);

dp = p_eps1(3)-p_eps1(1); % skip one value out of 2
z1 = zmin+dp/2;                          
p_k = [z1:dp:zmax]';

% raw matrix to fill
mat_eps = sparse(length(p_k), no_profile);
mat_K = sparse(length(p_k), no_profile);
mat_N2 = sparse(length(p_k), no_profile);
mat_T = sparse(length(p_k), no_profile);


N2_count = 0;
N2_min = -10;

%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%
for profile = 1:no_profile

    
    % -------------------- - T profile --------------------------- %
% $$$     %first time it is run
% $$$     fname = pfiles(profile, :);
% $$$     I = find(fname==' ');   
% $$$     fname(I) = [];
% $$$     load(fname)
% $$$     tfname = sprintf('%sT.dat', fname);
% $$$     dlmwrite(tfname, [P SBT], 'delimiter',' ','precision',6);
% $$$ 
    % second time it is run
    fname = pfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    tfname = sprintf('%sT.dat', fname);
    T = load(tfname);
    SBT = T(:,2) ;
    P = T(:,1);
    
    % bin to eps profile res.
    zbin = p_k(2)-p_k(1);
    for i = 1:length(p_k)
        I = find(P>(p_k(i)-zbin/2) & P<(p_k(i) + zbin/2));
        mat_T(i, profile) = nanmean(SBT(I, :), 1);
    end
    % ------------------------------------------------------------ %





    % --------------------- eps profile -------------------------- %
    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    load(fname)
    disp(fname)
    
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

   
% $$$     % ---- Remove unrealistic diffusivity ---- %
% $$$     I = find(log10(N2)<N2_min);
% $$$     if ~isempty(I)==1
% $$$         K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
% $$$         K_rho(I)=NaN;
% $$$         N2_count = N2_count+length(I);
% $$$     end
    
    % ---- store all profiles in matrix ---- %
    mat_eps(:, profile) = EPS; % those are 1-m binned
    mat_K(:, profile) = K_rho;
    mat_N2(:, profile) = N2;

    % ---------------------------------------------------------------- %

end 


% $$$ %%%%%%%%%%%%%%%%%
% $$$ % bin variables %
% $$$ %%%%%%%%%%%%%%%%%
% $$$ P_bin = [zmin+5/2:5:zmax-5/2]';
% $$$ zbin = 5;
% $$$ 
% $$$ mat_eps_bin = nan(length(P_bin), size(mat_eps, 2));
% $$$ mat_K_bin = mat_eps_bin;
% $$$ mat_N2_bin = mat_eps_bin;
% $$$ 
% $$$ for i = 1:length(P_bin)
% $$$     I = find(p_k>(P_bin(i)-zbin/2) & p_k<(P_bin(i)+zbin/2));
% $$$     mat_eps_bin(i, :) = nanmean(mat_eps(I, :), 1);
% $$$     mat_K_bin(i, :) = nanmean(mat_K(I, :), 1);
% $$$     mat_N2_bin(i, :) = nanmean(mat_N2(I, :), 1);
% $$$ end
% $$$ % Extract CIL depth interval and compute K and EPS (binned)
% $$$ for i = 1:size(mat_K_bin, 2)
% $$$     
% $$$     % select last 10m bins
% $$$     I=find(~isnan(mat_K_bin(:,i))==1);
% $$$     if max(P_bin(I))<150
% $$$     Kvec = mat_K(I(end-1):I(end),i);
% $$$     epsvec = mat_eps(I(end-1):I(end),i);
% $$$     Tvec = mat_T(I(end-1):I(end),i);
% $$$     
% $$$     K_CIL_bin(i) = nanmean(Kvec);
% $$$     eps_CIL_bin(i) = nanmean(epsvec);
% $$$     else
% $$$         K_CIL_bin(i) = NaN;
% $$$         eps_CIL_bin(i) = NaN;
% $$$     end
% $$$     
% $$$ % $$$     % select which in CIL
% $$$ % $$$     I = find(Tvec<=1);
% $$$ % $$$     K_CIL(i) = nanmean(Kvec(I));
% $$$ % $$$     eps_CIL(i) = nanmean(epsvec(I));
% $$$ end

% Extract CIL depth interval and compute K and EPS
%count=1;
for i = 1:size(mat_T, 2)
    
    % select last 10m bins
    I=find(~isnan(mat_K(:,i))==1);
    if max(I)<150 %(not at Riki)
        Kvec = mat_K(I(end-hb+1):I(end),i);
        epsvec = mat_eps(I(end-hb+1):I(end),i);
        Tvec = mat_T(I(end-hb+1):I(end),i);
        
% $$$         % Keep only CIL?
% $$$         I = find(Tvec>1);
% $$$         Kvec(I)=[];
% $$$         epsvec(I)=[];
% $$$     
% $$$     K_CIL(count:count+length(Kvec)-1) = Kvec;
% $$$     eps_CIL(count:count+length(epsvec)-1) = epsvec;
% $$$     count = count +length(I);
    
        K_CIL(i) = nanmean(Kvec);
        eps_CIL(i) = nanmean(epsvec);
    else
        K_CIL(i) = NaN;
        eps_CIL(i) = NaN;
    end
    
end

% remove NaNs
I = find(isnan(K_CIL)==1);
K_CIL(I)=[];
eps_CIL(I)=[];

I = find(isnan(eps_CIL)==1);
K_CIL(I)=[];
eps_CIL(I)=[];

%%%%%%%%%%%%%%%%%%%%
% --- bootstrap --- %
%%%%%%%%%%%%%%%%%%%%

% -- bootstrap -- %
disp('bootstrap...')
N = length(K_CIL);

% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    
    eps_boot_b(b) = nanmean(eps_CIL(r));
    K_boot_b(b) = nanmean(K_CIL(r));
end

% mean of random sampling
eps_boot_dot = nanmean(eps_boot_b);
K_boot_dot = nanmean(K_boot_b);

% Compute 95% confidence interval
eps_sort = sort(eps_boot_b);
K_sort  = sort(K_boot_b);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

eps_2p5 = eps_sort(CI_2p5);
K_2p5 = K_sort(CI_2p5);
eps_97p5 = eps_sort(CI_97p5);
K_97p5 = K_sort(CI_97p5);

eps_ave = eps_boot_dot;
K_ave = K_boot_dot;


%%%%%%%%%%%%%%%%%%%%%%
% - within the CIL - %
%%%%%%%%%%%%%%%%%%%%%%

disp({'\epsilon_{b}...'})
[eps_ave eps_2p5 eps_97p5]

disp({'K_{b}...'})
[K_ave K_2p5 K_97p5]



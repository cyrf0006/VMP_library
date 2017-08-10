function ENK_CIL(profile_names, epsprofile_names)

% function ENK_CIL(profile_names, epsprofile_names)
%
% To compute mean dissipation and diffusivity within the CIL. This
% is new comprage to the time the average was done by boot_VMP.m on
% a fixed depth range (50-150m). Now the The average is done only
% on bins colder than 1degC
%    
%
% usage ex: ENK_CIL('Tfile_names_veryall', 'file_names_veryall')
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
N2_min = -6;

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

    % ---------------------------------------------------------------- %

end 


% Extract CIL depth interval and compute K and EPS
count=1;
for i = 1:size(mat_T, 2)
    I = find(mat_T(:,i)<=1);
    
    % average on profiles
    K_CIL(i) = nanmean(mat_K(I,i));
    eps_CIL(i) = nanmean(mat_eps(I,i));
    % average on bins
% $$$     K_CIL(count:count+length(I)-1) = mat_K(I,i);
% $$$     eps_CIL(count:count+length(I)-1) = mat_eps(I,i);
% $$$     count = count +length(I);
% $$$     
% for visual inspection
% $$$     subplot(1,2,1)
% $$$     plot(mat_T(I,i), p_k(I))
% $$$     set(gca, 'ydir', 'reverse')
% $$$     subplot(1,2,2)
% $$$     semilogx(mat_K(I,i), p_k(I))
% $$$     set(gca, 'ydir', 'reverse')
% $$$     pause

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

disp('\epsilon_{CIL}...')
[eps_ave eps_2p5 eps_97p5]

disp('K_{CIL}...')
[K_ave K_2p5 K_97p5]

keyboard
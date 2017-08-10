function Ri_casts(adcp_s2, eps_profiles, range, varargin)
    
% function Ri_casts(adcp_s2, eps_profiles, range, varargin)
%
% ex: Ri_casts(merged_adcp_shear.mat, file_names_all_one [20 330], merged_adcp_shear_bis.mat, file_names_all_bis)
%
%
%
%
% Varargin test
if isempty(varargin)==1
    bis = 0;
elseif size(varargin,2) == 2
    bis = 1;
    adcp_s2_bis = varargin{1};
    eps_profiles_bis = varargin{2};        
else
    disp('vargin size must be 2 or nothing')
    disp('extra arguments ignored')
    bis = 0;
end


load(adcp_s2);


% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 350;
nboot = 500;    
P_bin = P_S2;
zbin = P_bin(2)-P_bin(1);

thresh = 1200; %sec

% load *.P files names (file in which are recorded .P files)
fid = fopen(eps_profiles);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

no_profile = size(epsfiles, 1); %number of eps_files 

Ri = [];
epsilon = [];
Re_buoy = [];

% loop on profile
count  = 1;
for profile = 1:no_profile

    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
    
    P_bin = P_S2;  
    [Y, I] = min(abs((mtime_eps(1)-t_S2)));
    if Y > thresh/86400;
        continue
    else
        % only keep wanted depth range
        J = find(isnan(S2(:,I))==1);
        P_bin(J) = [];       
        J = find(P_bin > range(1) & P_bin < range(2));
        P_bin = P_bin(J);
        
        % keep corresponding shear
        s2_bin = S2(J, I);
        
        % raw matrix to fill
        eps_bin = nan(length(P_bin), 1);
        N2_bin = eps_bin;
  
               
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
        
        % ---------- bin profile -------------- %
        for i = 1:length(P_bin)
            I = find(p_k >= P_bin(i)-zbin/2 & p_k <= P_bin(i)+zbin/2);
            eps_bin(i) = nanmean(EPS(I));
            N2_bin(i) = nanmean(N2(I));
        end
        % ------------------------------------ %
        
        % *** Here we have s2_bin, eps_bin and N2_bin *** %        
        Ri_bin = N2_bin./s2_bin;
        Re_bin = eps_bin./(1e-6.*N2_bin);
        %        [length(Ri_bin) length(Re_bin)]
        %pause(.1)
        Ri(count:count+length(s2_bin)-1) = Ri_bin;
        epsilon(count:count+length(s2_bin)-1) = eps_bin;
        Re_buoy(count:count+length(s2_bin)-1) = Re_bin;
        count = count+length(s2_bin);
                        
    end
end 

if bis == 1
    Ri_casts_varargin
end






figure(1)
loglog(Ri, epsilon, '.k')
%semilogx(epsilon, tanh(Ri), '.k')
xlim([0 1e5])

figure(2)
%loglog(Ri, Re_buoy, '.k')


%smooth_Re = loess(Ri, epsilon, sort(Ri), 0.3, 1);
%Ri_sort = sort(Ri)

keyboard
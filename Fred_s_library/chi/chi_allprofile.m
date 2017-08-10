function chi_allprofile(profiles)

% function chi_allprofile(profiles)
%
% usage ex: chi_allprofile('Tprofiles_all')
%  Where 'prof_files' comes from:
%     "ls -1 profile*.mat | sed 's/\.mat//' > prof_files"
%
% This function is similar to eps_profile.m, but compute this time
% 'chi' from micro_temperature sensors. This script can also
% compute epsilon from chi and spectral fit to batchelor
% spectrum. Further improvements are to come
%
% F. Cyr - Feb./March 2012.
% ---------------------------------------------------------- %
tic

fid = fopen(profiles);
C = textscan(fid, '%s', 'delimiter', '\n');
pro_files = char(C{1});

no_profiles = size(pro_files,1);


% load 1st file
fname_in = pro_files(1, :);
I = find(fname_in==' ');
fname_in(I)=[];
load(fname_in)
dt1dz = diff(t1)./diff(p);
dt2dz = diff(t2)./diff(p);
if nansum(abs(dt2dz)) > 1e6 & nansum(abs(dt1dz)) > 1e6 %this is weak!
    disp('both microT absent, nothing to do!')
    return
elseif nansum(abs(dt2dz)) > 10*nansum(abs(dt1dz))
    disp('microT2 absent!')
    probe = 1;
elseif nansum(abs(dt1dz)) > 10*nansum(abs(dt2dz))
    disp('microT1 absent!')
    probe = 2;
else
    disp('both microT1 present!')
    probe = 0;
end


for count = 1:no_profiles
    

    fname_in = pro_files(count, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    fname_out = ['chi_' fname_in];
    
    
% $$$     if count > 999
% $$$         fname_out = sprintf('chi_profile%d', count);
% $$$     elseif count > 99
% $$$         fname_out = sprintf('chi_profile0%d', count);
% $$$     elseif count > 9
% $$$         fname_out = sprintf('chi_profile00%d', count);
% $$$     else
% $$$         fname_out = sprintf('chi_profile000%d', count);
% $$$     end
    
    %    load(['../../' fname_in]);

    %  -- Uncomment to consider only profiles that hit the bottom -- %
    % (Since 'most' profile at border hit the bottom in 2011... )
% $$$     load(['../../' fname_in]);
% $$$     if P(end)>120 
% $$$         continue
% $$$     end
    % -------------------------------------------------------------- %
    
    load(fname_in)
    
    d=sprintf('PROCESSING PROFILE %d ...', count);
    disp(d)
    
    % Angle and viscosity
    ax = ang2acc(pitch); ay = ang2acc(roll);
    DENS = sw_dens(SBS, SBT, P);
    [nu, mu] = viscosity(SBS, SBT, DENS);

    % round frequencies
    fs = round(fs);
    FS = round(FS);

    dt1dz = diff(t1)./diff(p);
    dt2dz = diff(t2)./diff(p);

    % ----- > THIS SHOULD BE THE REAL CALL....
% $$$         [chi1, eps1, p_chi1, ratio] = chi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1);
% $$$         %[chi2, eps2, p_chi2, ratio] = chi_vmp(dt2dz,[ax ay az],P, W, nu, fs, FS, 1);

    % ------ >... BUT FOR THE MOMENT, WE ONLY COMPUTE CHI IN ORDER TO
    % SPEED THE TREATMENT

    if probe == 0
        [chi1, p_chi1] = quickchi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1); 
        [chi2, p_chi2] = quickchi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1); 
    elseif probe == 1;
        [chi1, p_chi1] = quickchi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1); 
        chi2 = nan(size(chi1));
        p_chi2 = p_chi1;
    else
        [chi2, p_chi2] = quickchi_vmp(dt2dz,[ax ay az],P, W, nu, fs, FS, 1); 
        chi1 = nan(size(chi2));
        p_chi1 = p_chi2;
    end
    
    % dont compute epsilon for the moment
    eps1 = nan(size(chi1)); eps2 = nan(size(chi1)); ratio = nan;
    
    mtime_eps = interp1(P,MTIME,p_chi1,'linear','extrap');

    save(fname_out,'p_chi1','p_chi2','eps1','eps2', 'mtime_eps','chi1' ,'chi2', 'ratio');    

    
end

toc
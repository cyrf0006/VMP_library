function chi_profile_count(profiles)

% function chi_profile_count(profiles)
%
% Same as chi_profile.m, but with counter renaming in local folder
%
% F. Cyr - Nov. 2012
% ---------------------------------------------------------- %
tic
fid = fopen(profiles);
C = textscan(fid, '%s', 'delimiter', '\n');
pro_files = char(C{1});

siz = size(pro_files);
no_profiles = siz(1); %number of profile* files 
% $$$ end


for count = 1:no_profiles
    fname_in = pro_files(count, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    %        fname_out = ['chi_' fname_in];
    if count<10
        fname_out = sprintf('chi_profile_hitbottom00%d.mat', count);
    else
        if count<100
            fname_out = sprintf('chi_profile_hitbottom0%d.mat', count);
        else %profile>100
            fname_out = sprintf('chi_profile_hitbottom%d.mat', count);
        end
    end
% $$$     end
    
    load(fname_in);
    
    d=sprintf('TREATMENT OF PROFILE %d ...', count);
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
    
    if sum(~isnan(dt1dz))==0
        continue
    else
        [chi1, p_chi1] = quickchi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1); 
        chi2 = nan(size(chi1)); eps1 = nan(size(chi1)); eps2 = ...
               nan(size(chi1)); ratio = nan;
        p_chi2 = p_chi1; 
    end            
    
    
    mtime_eps = interp1(P,MTIME,p_chi1,'linear','extrap');

    save(fname_out,'p_chi1','p_chi2','eps1','eps2', 'mtime_eps','chi1' ,'chi2', 'ratio');    

    
end


toc
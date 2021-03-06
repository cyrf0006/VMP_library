function chi_profile(profiles)

% function chi_profile(profiles)
%
% usage ex: chi_profile('prof_files')
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
% Choose if the function epsilon_vmp.m is explicit or not
% $$$ answer=0;
% $$$ while answer==0
% $$$     R1 = input('Do you want visual inspection of profiles when suspect? (y/n) ', 's');
% $$$     if strcmp('y',R1)==1 
% $$$         explicit=1;
% $$$         answer=1;
% $$$     elseif strcmp('n',R1)==1
% $$$         explicit=0;
% $$$         answer=1;
% $$$     else
% $$$         disp('bad answer! please type y or n')
% $$$     end
% $$$ end

% $$$ 
% $$$ if ischar(profiles)==0
% $$$     no_profiles = profiles;
% $$$ else
    % load profiles*.mat
    fid = fopen(profiles);
    C = textscan(fid, '%s', 'delimiter', '\n');
    pro_files = char(C{1});

    siz = size(pro_files);
    no_profiles = siz(1); %number of profile* files 
% $$$ end


for count = 1:no_profiles
    
% $$$     if ischar(profiles)==0 % the input is a number
% $$$         if count<10
% $$$             fname_in=sprintf('profile00%d', count);
% $$$             fname_out=sprintf('chi_profile00%d', count);
% $$$         else
% $$$             if count<100
% $$$                 fname_in=sprintf('profile0%d', count);
% $$$                 fname_out=sprintf('chi_profile0%d', count);
% $$$             else %no_files>100
% $$$                 fname_in=sprintf('profile%d', count);
% $$$                 fname_out=sprintf('chi_profile%d', count);
% $$$             end
% $$$         end
% $$$     else % the input is a list.

        fname_in = pro_files(count, :);
        I = find(fname_in==' ');
        fname_in(I)=[];
        fname_out = ['chi_' fname_in];
        
% $$$         if count<10
% $$$         fname_out = sprintf('chi_profile_hitbottom00%d.mat', count);
% $$$     else
% $$$         if count<100
% $$$             fname_out = sprintf('chi_profile_hitbottom0%d.mat', count);
% $$$         else %profile>100
% $$$             fname_out = sprintf('chi_profile_hitbottom%d.mat', count);
% $$$         end
% $$$     end
% $$$      fname_out = ['chi_profile' fname_in];
% $$$     end
    
        load(fname_in);

    % test to see if a shear probe is absent
% $$$     if count==1
% $$$         load(fname_in);
% $$$         
% $$$ % $$$         if nanmean(shear1.^2)>0.001 & mean(shear2.^2)>0.001
% $$$ % $$$             disp('both shear probe present')
% $$$ % $$$         elseif nanmean(shear1.^2)>0.001
% $$$ % $$$             shearprobe = 's1';
% $$$ % $$$             disp('S2 absent!?')
% $$$ % $$$         elseif nanmean(shear2.^2)>0.001
% $$$ % $$$             shearprobe = 's2';
% $$$ % $$$             disp('S1 absent!?')
% $$$ % $$$         else
% $$$ % $$$             disp('Both shear probes absent... nothing to do!')
% $$$ % $$$             return
% $$$ % $$$         end
% $$$     end
    
    
    d=sprintf('TREATMENT OF PROFILE %d ...', count);
    disp(d)

    % call shear analysis
% $$$     if exist('shearprobe')==0
% $$$         shear_analysis(fname_in,fname_out, 1, explicit);
% $$$     else
% $$$         shear_analysis(fname_in,fname_out, 1, explicit, shearprobe);
% $$$     end

    
    %chi_analysis(fname_in,fname_out); % replaced by following
    
    % Angle and viscosity
    ax = ang2acc(pitch); ay = ang2acc(roll);
    DENS = sw_dens(SBS, SBT, P);
    [nu, mu] = viscosity(SBS, SBT, DENS);

    % round frequencies
    fs = round(fs);
    FS = round(FS);

    %   if isempty(varargin)==1 % both shear probes are presents!
                            
        % ---------------------- %
        % preanalysis (cleaning, flagging, etc.)
        % ---------------------- %
        % nopreanalysis required?
% $$$         sh1 = shear_preanalysis(shear1, w, pitch, roll, p, 5);
% $$$         sh2 = shear_preanalysis(shear2, w, pitch, roll, p, 5);

        % ---------------------- %
        % beginning treatment
        % ---------------------- %
% $$$     [eps1,p_eps1] = epsilon_vmp(sh1,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$     [eps2,p_eps2] = epsilon_vmp(sh2,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
        
        % compute micro dt/dz
        % disp('here!')
        %keyboard
        dt1dz = diff(t1)./diff(p);
        dt2dz = diff(t2)./diff(p);

        
        % ----- > THIS SHOULD BE THE REAL CALL....
% $$$         [chi1, eps1, p_chi1, ratio] = chi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1);
% $$$         %[chi2, eps2, p_chi2, ratio] = chi_vmp(dt2dz,[ax ay az],P, W, nu, fs, FS, 1);

        % ------ >... BUT FOR THE MOMENT, WE ONLY COMPUTE CHI IN ORDER TO
        % SPEED THE TREATMENT
 
        if sum(~isnan(dt1dz))==0
            continue
        else
            [chi1, p_chi1] = quickchi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1); 
            chi2 = nan(size(chi1)); eps1 = nan(size(chi1)); eps2 = ...
                   nan(size(chi1)); ratio = nan;
            p_chi2 = p_chi1; 
        end            
        
        
        
        
        
% $$$         disp('rendu ici!')
% $$$         keyboard
        %[eps2,p_eps2] = epsilon_vmp_fc(sh2,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);

% $$$     elseif strcmp(varargin{1}, 's1')==1;%one shear probe is absent!
% $$$         sh1 = shear_preanalysis(shear1, w, pitch, roll, p, 5);
% $$$ % $$$     [eps1,p_eps1] = epsilon_vmp(sh1,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$         [eps1,p_eps1] = epsilon_vmp_fc(sh1,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$         eps2 = nan(1,length(eps1));
% $$$         p_eps2 = eps2;
% $$$     else
% $$$         sh2 = shear_preanalysis(shear2, w, pitch, roll, p, 5);
% $$$ % $$$     [eps2,p_eps2] = epsilon_vmp(sh2,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$         [eps2,p_eps2] = epsilon_vmp_fc(sh2,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$         eps1 = nan(1,length(eps2));
% $$$         p_eps1 = eps1;
% $$$     end


% $$$     % ---------------------- %
% $$$     % post-treatment
% $$$     % ---------------------- %
% $$$     if isempty(eps1)==1 & isempty(eps2)==1
% $$$         disp('Shear profile empty, skip profile');
% $$$         disp('Press any key')
% $$$         pause
% $$$         return
% $$$     end
% $$$ 
% $$$     % mean profile
% $$$     eps = nanmean([eps1; eps2]);
% $$$     p_eps = nanmean([p_eps1; p_eps2]);
% $$$ 
% $$$     % buoyancy freq.
% $$$     RHO = sw_dens(SBS,SBT,P);
% $$$     N = buoy_freq(RHO,P,p_eps1,zbin);
% $$$ 
% $$$     % Bin the fluoro and trans data if they exist.
% $$$     if exist('fluoro') & exist('trans')
% $$$         for k = 1:length(p_eps1)
% $$$             I = find(p >= (p_eps1(k) - zbin/2) & p < (p_eps1(k) + zbin/2));  
% $$$             TRANS(k) = nanmean(trans(I));
% $$$             FLUORO(k) = nanmean(fluoro(I)); 
% $$$         end
% $$$     end
% $$$ 
% $$$ 
% $$$     % Calculate the Ozmidov scale
% $$$     Lo = sqrt(eps ./ (N.^3));
% $$$ 
% $$$     nu = 10^(-6);
% $$$     isotropy_index = eps./(nu.*(N.^2));
% $$$ 
% $$$     % Interpolate time at p_eps
    mtime_eps = interp1(P,MTIME,p_chi1,'linear','extrap');

    save(fname_out,'p_chi1','p_chi2','eps1','eps2', 'mtime_eps','chi1' ,'chi2', 'ratio');    

    
end


toc
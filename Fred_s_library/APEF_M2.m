function APEF_M2(N2_matrix, tprofiles, ovt_files, tide_file, zhab)
    
% function Lo_Lt(prof_files, ovt_files)  
%
% (to be run in ~/WINDEX/data_processing/BMix_study)
%
% Inspired from Lo_Lt('Tprofiles_riki', 'fineovt_riki')
% and N2_M2composite(file_names, zbin, tidefile)
%
% usage ex: 
% APEF_M2('N2.mat', 'hit_bottom_tprofiles', 'ovt_files', 'tide_2009-2012.dat', 20)
% APEF_M2('N2_with2010.mat', 'hit_bottom_tprofiles_with2010', 'ovt_files_with2010', 'tide_2009-2012.dat', 20)
% APEF_M2('N2_slope.mat', 'hit_bottom_tprofiles_slope', 'ovt_files_slope', 'tide_2009-2012.dat', 20)
% APEF_M2('N2_20110928.mat', 'hit_bottom_tprofiles_20110928', 'ovt_files_20110928', 'tide_2009-2012.dat', 20)
% 
% **User must run N2_M2composite.m to generate N2_mat.m   
% **User must also run mat2reorange.m to generate 'ovt_files_slope' (see help menu)
%    and in a shell do:
%   for i in fine*.dat; do 
%   reorange $i; 
%   done
%
% F. Cyr, Nov. 2012
    
% few params:
g = 9.81; %m/s2
nu = 1e-6;% m2/s
Kt = 1.4e-7; %m2/s

% profiles
fid = fopen(tprofiles);
C = textscan(fid, '%s', 'delimiter', '\n');
list = char(C{1});

% overturns
fid = fopen(ovt_files);
C = textscan(fid, '%s', 'delimiter', '\n');
list_ovt = char(C{1});

% load stratification (just for proftime!)
load(N2_matrix)
no_profiles = length(proftime);
if no_profiles ~= size(list_ovt, 1)
    disp('Error, no. profiles is not equal to no. overturns')
    exit
end

% initialize OVT matrix
OVT = [];

count=1;
for i = 1:no_profiles

    disp(sprintf('profile %d', i));
    % ---- Work on VMP profiles ---- %
    fname_in = list(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    load(fname_in)   
    ax = ang2acc(pitch); ay = ang2acc(roll);
    % --------------------------------- %
    
    % -------- Work on overturns ------ %
    fname_ovt = list_ovt(i, :); % remove white space
    I = find(fname_ovt==' ');
    fname_ovt(I)=[];

    % get overturns for this profile
    command = sprintf('cut -d" " -f2 %s > /tmp/tmp', fname_ovt);
    system(command);
    ind1 = load('/tmp/tmp');

    command = sprintf('cut -d" " -f3 %s > /tmp/tmp', fname_ovt);
    system(command);
    ind2 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f4 %s > /tmp/tmp', fname_ovt);
    system(command);
    z1 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f5 %s > /tmp/tmp', fname_ovt);
    system(command);
    z2 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f6 %s > /tmp/tmp', fname_ovt);
    system(command);
    L_t = load('/tmp/tmp'); 

% $$$     command = sprintf('cut -d" " -f7 %s > /tmp/tmp', fname_ovt);
% $$$     system(command);
% $$$     ThFluct = load('/tmp/tmp'); % <kg/m3>, cf. rho'
% $$$     
% $$$     command = sprintf('cut -d" " -f9 %s > /tmp/tmp', fname_ovt);
% $$$     system(command);
% $$$     L_t_on_H = load('/tmp/tmp');
% $$$     H = L_t./L_t_on_H;
    
    command = sprintf('cut -d" " -f11 %s > /tmp/tmp', fname_ovt);
    system(command);
    APEF = load('/tmp/tmp');
    
    if isempty(L_t)==1
        continue
    end
      
    no_ovt = length(L_t);
    for j = 1:no_ovt
 
        %if z2(j)-z1(j) < .1 | (zbottom(i)-z2(j) > zhab)% overturn at least 10cm
        if (zbottom(i)-z2(j) > zhab)% overturn at least 10cm
            continue
        end
        
        
        % Few calculations (for Ra & Jb)
        Hz = P(ind2(j))-P(ind1(j)); % overturn vertical size
        rho = sw_dens(SBS(ind1(j):ind2(j)), SBT(ind1(j):ind2(j)), P(ind1(j):ind2(j)));
        rho_sort = sort(rho);
        rho_p = rho_sort(end)-rho_sort(1); % density anomaly
        pp = polyfit(P(ind1(j):ind2(j)),rho_sort,1);
        rho0 = nanmean(rho);
        N = sqrt((g/rho0)*pp(1)); % buoy freq on overturn
        Jb = APEF(j).*N;
        Ra = (g*rho_p*Hz.^3)./(rho0*Kt*nu);
        Ts = (Hz.^2./Jb).^(1/3);
        
        if z1(j) < p(1) % top of watercolumn
            epsilon = NaN;
            gamma = NaN;
        else
            % epsilon calculation on patch
            if  z2(j)-z1(j) < .1% overturn at least 10cm
                epsilon = NaN;
                gamma = NaN;
            else
                I = find(p>=z1(j) & p<=z2(j));
                Ifine = find(P>=z1(j) & P<=z2(j)); 
                eps1 = patch_epsilon(shear1(I), [ax(I), ay(I), az(I)], nanmean(W(Ifine)), nu);
                eps2 = patch_epsilon(shear2(I), [ax(I), ay(I), az(I)], ...
                                     nanmean(W(Ifine)), nu);
            end
            
            if abs(eps1-eps2)<100
                epsilon = nanmean([eps1 eps2]);
                gamma = Jb/epsilon;
            else
                disp('big difference between eps1 and eps2 in APEF.m [K>>]')
                keyboard
            end
        end
        
        %OVT contains: [profindex hab1 hab2, L_t APEF, L_t ovtSize, N(brunt-V), Jb(buoFlux), Raleigh epsilon gamma]
        if isempty(OVT)==1
            %zbottom come from 'N2_matrix', thus zbottom(i)-z1(j) is hab
            OVT = [i zbottom(i)-z1(j) zbottom(i)-z2(j) L_t(j) APEF(j), ...
                   Hz, N, Jb, Ra, epsilon, gamma, Ts];   
        else
            OVT = [OVT; [i zbottom(i)-z1(j) zbottom(i)-z2(j) L_t(j) ...
                         APEF(j), Hz, N, Jb, Ra, epsilon, gamma, Ts]]; 
        end
    end                
    % ------------------------------- % 
end

P_HR = 0:.05:zhab; % already hab since OVT = [i zbottom(i)-z1(j) zbottom(i)-z2(j)  ...                               
APEF_mat = nan(length(P_HR), length(proftime));
Hz_mat = nan(length(P_HR), length(proftime));
N_mat = nan(length(P_HR), length(proftime));
Jb_mat = nan(length(P_HR), length(proftime));
Ra_mat = nan(length(P_HR), length(proftime));


for i = 1:size(OVT, 1)
    I = find(P_HR>OVT(i, 3) & P_HR<OVT(i,2));    
    APEF_mat(I,OVT(i,1)) =  OVT(i,4);
    Hz_mat(I,OVT(i,1)) =  OVT(i,6);
    N_mat(I,OVT(i,1)) =  OVT(i,7);
    Jb_mat(I,OVT(i,1)) =  OVT(i,8);
    Ra_mat(I,OVT(i,1)) =  OVT(i,9);   
end


% Here we save matrices so ~/LaTeX/MS/BMix/matlab_files/graphics.m
% will be able to use it.
proftime_ovt = proftime;
%save partII.mat APEF_mat Hz_mat N_mat Jb_mat Ra_mat P_HR proftime_ovt
%save partII_20110928.mat proftime_ovt OVT
save partII.mat proftime_ovt OVT


% Above is needed by graphics.m
% Below is for APEF_M2.m use alone

% $$$ 
% $$$ % -------------- Compute time to high tide ------------ %
% $$$ tide  = load(tide_file);
% $$$     
% $$$ mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
% $$$ level = tide(:,6);
% $$$ 
% $$$ % find high tide time 
% $$$ count = 1;
% $$$ for i = 2:length(mtime)-1
% $$$ 
% $$$     if level(i)>level(i-1) & level(i)>level(i+1)
% $$$         T(count) = mtime(i); % high tide time
% $$$         L(count) = level(i); % high tide level
% $$$         count = count+1;
% $$$     end
% $$$ end
% $$$ % T contains the hour of each high tide
% $$$ clear mtime
% $$$ 
% $$$ for i = 1:length(proftime)
% $$$     [Y, I] = min(abs(proftime(i)-T));
% $$$     A(i) = (proftime(i)-T(I))*24;
% $$$     B(i) = L(I); %level of the closest hightide
% $$$ end
% $$$ % ------------------------------------------------------ %
% $$$ 
% $$$ time2 = A;
% $$$ 
% $$$ dtide = 1;
% $$$ reg_tide = -6:dtide:6;
% $$$ 
% $$$ 
% $$$ APEF_tide = zeros(length(P_HR), length(reg_tide));
% $$$ Hz_tide = zeros(length(P_HR), length(reg_tide));
% $$$ N_tide = zeros(length(P_HR), length(reg_tide));
% $$$ Jb_tide = zeros(length(P_HR), length(reg_tide));
% $$$ Ra_tide = zeros(length(P_HR), length(reg_tide));
% $$$ 
% $$$ for i = 1: length(reg_tide)
% $$$     I = find(time2 > reg_tide(i) - dtide & time2 < reg_tide(i) + dtide);  
% $$$     APEF_tide(:,i) = nanmean(APEF_mat(:,I), 2);
% $$$     Hz_tide(:,i) = nanmean(Hz_mat(:,I), 2);
% $$$     N_tide(:,i) = nanmean(N_mat(:,I), 2);
% $$$     Jb_tide(:,i) = nanmean(Jb_mat(:,I), 2);
% $$$     Ra_tide(:,i) = nanmean(Ra_mat(:,I), 2);
% $$$ end
% $$$ 
% $$$ 
% $$$ % $$$ figure(1)
% $$$ % $$$ clf
% $$$ % $$$ imagesc(reg_tide, P_HR, Jb_tide)
% $$$ % $$$ set(gca, 'ydir','normal')
% $$$ % $$$ xlabel('Time to hightide (hours)')
% $$$ % $$$ ylabel('hab (m)')
% $$$ % $$$ title('Available potential energy of fluctuations (m^2 s^{-2})')
% $$$ % $$$ colorbar
% $$$ % $$$ caxis([.25 .5])
% $$$ % $$$ keyboard
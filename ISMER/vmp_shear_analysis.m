function vmp_shear_analysis(p_fname_in, p_fname_out, zbin, explicit, varargin);

% function vmp_shear_analysis(p_fname, p_fname_out, zbin)
%
% where 
% - p_fname: the input profile name 
% - p_fname_out: the output epsilon profile name
% - zbin: bin size for epsilon average (will be computed every
% zbin/2). !!!!WE RECOMMAND TO USE zbin=1!!!, because otherwise bins are
% not beginning at the same depth for every profile...
% - explicit: Parameter to choose if you can explicitly choose
% which are the bins to reject. Can be 1 or 0. 
%    1: you will be able to reject suspicious bins after visual
%    inspection
%    0: you'll be able to have a visual inspection, but you won't
%    be able to reject bins.
%
%ex: shear_analysis('profile001', 'eps_profile001', 1, 1)
%
% This function save in a .mat both shear1 and shear2 with their
% pressure [ p1 shear1 p2 shear2 ]
%           ...  ...  ...   ... 
%
%
% Author: Frederic Cyr - 2009/12/17
%
%     Modifications: 
%         F. Cyr (2010/05/26)
%             -> now zbin is a parameter
%         D. Bourgault - (2010?)
%             -> this functions returns now more variables (N,
%             MTIME, Lo, isotropy_index, fluoro...)
%         F. Cyr (2010-05-31) 
%             -> add param "explicit"
%             -> add explicit param to epsilon_vmp.m
%         F. Cyr (2010-10-04)
%             -> test to see if both shear probes were present
%             (eliminate unusefull treatment when S2 absent)
%             Uses varargin.
%             -> modified Dan's mean epsilon with nanmean
%         F. Cyr (2010-11-23)
%             -> Check if eps1 and eps2 are not empty, otherwise
%             skip (in posttreatment section)
%              
%
%% -------------------------------------------------------------- %%

% load selected profile
load(p_fname_in);

% ---------------------- %
% preliminary options
% ---------------------- %

% Options for Dewey's scripts
%zbin = 1.0; % bin size
%zbin = 4.0;
%zbin = 5.0;
%iplt = [1 2 3 0]; %iplt==[figNum1 figNum2 figNum3 spause]
iplt = 0; %   iplt==0: No plotting done.

% Angle and viscosity
ax = ang2acc(pitch); ay = ang2acc(roll);

DENS = sw_dens(SBS, SBT, P);
[nu, mu] = viscosity(SBS, SBT, DENS);

% round frequencies
fs = round(fs);
FS = round(FS);

if isempty(varargin)==1 % both shear probes are presents!
    % ---------------------- %
    % preanalysis (cleaning, flagging, etc.)
    % ---------------------- %
    sh1 = vmp_shear_preanalysis(shear1, w, pitch, roll, p, 5);
    sh2 = vmp_shear_preanalysis(shear2, w, pitch, roll, p, 5);
    % ---------------------- %
    % beginning treatment
    % ---------------------- %
% $$$     [eps1,p_eps1] = epsilon_vmp(sh1,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$     [eps2,p_eps2] = epsilon_vmp(sh2,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$     [eps1,p_eps1] = epsilon_vmp_fc(sh1,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% $$$     [eps2,p_eps2] = epsilon_vmp_fc(sh2,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
    [eps1,p_eps1] = vmp_spectral_integration(sh1,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
    [eps2,p_eps2] = vmp_spectral_integration(sh2,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
elseif strcmp(varargin{1}, 's1')==1;%one shear probe is absent!
    sh1 = vmp_shear_preanalysis(shear1, w, pitch, roll, p, 5);
% $$$     [eps1,p_eps1] = epsilon_vmp(sh1,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
    [eps1,p_eps1] = vmp_spectral_integration(sh1,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
    eps2 = nan(1,length(eps1));
    p_eps2 = eps2;
else
    sh2 = vmp_shear_preanalysis(shear2, w, pitch, roll, p, 5);
% $$$     [eps2,p_eps2] = epsilon_vmp(sh2,[ax ay az],fs,SBT,P,FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
    [eps2,p_eps2] = vmp_spectral_integration(sh2,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
    eps1 = nan(1,length(eps2));
    p_eps1 = eps1;
end


% ---------------------- %
% post-treatment
% ---------------------- %
if isempty(eps1)==1 & isempty(eps2)==1
    disp('Shear profile empty, skip profile');
    disp('Press any key')
    pause
    return
end

% mean profile
eps = nanmean([eps1; eps2]);
p_eps = nanmean([p_eps1; p_eps2]);

% buoyancy freq.
RHO = sw_dens(SBS,SBT,P);
N = buoy_freq(RHO,P,p_eps1,zbin);

% Bin the fluoro and trans data if they exist.
if exist('fluoro') & exist('trans')
  for k = 1:length(p_eps1)
    I = find(p >= (p_eps1(k) - zbin/2) & p < (p_eps1(k) + zbin/2));  
    TRANS(k) = nanmean(trans(I));
    FLUORO(k) = nanmean(fluoro(I)); 
  end
end


% Calculate the Ozmidov scale
Lo = sqrt(eps ./ (N.^3));

nu = 10^(-6);
isotropy_index = eps./(nu.*(N.^2));

% Interpolate time at p_eps
mtime_eps = interp1(P,MTIME,p_eps,'linear','extrap');

if exist('fluoro') & exist('trans')
  save(p_fname_out,'p_eps1','p_eps2','eps1','eps2','N','mtime_eps','Lo','isotropy_index','FLUORO','TRANS');    
else 
  save(p_fname_out,'p_eps1','p_eps2','eps1','eps2','N','mtime_eps','Lo','isotropy_index');
end

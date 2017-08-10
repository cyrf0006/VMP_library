function N = buoy_freq(rho,z,Z,varargin);
% function N = buoy_freq(rho,z,Z,varargin);
%  
% usage: whether N = buoy_freq(rho,z,Z)
%       or N = buoy_freq(rho,z,Z, zbin) (old version)    
% 
% This function returns a background buoyancy frequency profile centered at 
% Z given a raw non-uniform density profile rho(z).
%
% The profile is first sorted. The derivative d(rho)/dz over the desired 
% interval zbin is calculated with a linear fit.
%

% Author: Daniel Bourgault - 2010/04/16
%      
%  modifs: 
%      - 2013-04-04 (F. Cyr): works if Pbin spans depth range not in P
%      - 2014-01-08 (F. Cyr): remove default zbin input and add
%              varargin instead (can compute zbin from Z)
%      
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin)==1
    zbin = Z(2)-Z(1);
elseif size(varargin,2)==1
    zbin = varargin{1};
else
    disp('Wrong input... try "help buoy_freq"')
    return
end

g = 9.81;

if nanmean(rho) < 1000
  display('It is likely that the function buoy_freq was fed with sigma instead of rho.');
  display('N values may be off by a factor of 1000. Double check this.');
end

rho = sort(rho);

for k = 1:length(Z)
    I = find(z >= (Z(k) - zbin/2) & z <= (Z(k) + zbin/2));  
 
    rho0 = mean(rho(I));
    
    if size(I) < 2 | isnan(rho0) == 1 
        N(k) = NaN;
    else        
        rho0 = mean(rho(I));
        p = polyfit(z(I),rho(I),1);
        N(k) = sqrt((g/rho0)*p(1));
    end
end

function indices = vmp_splitProfiles(MAT_fname);

% function indices = split_profiles(MAT_fname);
%
% This script splits profiles from a .mat file containing successive
% profiles.
%
% The parameter is the name of the .MAT file without the extension (ex.
% 'DAT000')
%
% The function returns a matrix containing in first column the profile
% number and in 2nd and 3rd columns the indice of the first and last point
% of the profile. [no_profile i_min i_max  flagg]
%                   ...        ...    ...   ...
% 
% **** Note that the returned indices are the one of the slow matrix ****
% **** Note also that a more restrictive function is written as
% split_profiles2.m. The probes speed perturbations are smaller, but we
% losse many meters at the beginning and end of the profile. ****
%
% Author: Frédéric Cyr - Feb. 2013
%
% Note that this is the new version of the former
% 'split_profiles.m'. History of this function can be found on the
% SVN deposit (under history of vmp_splitProfiles.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wmax = 1.0;
Wmin = 0.3;
Pmin = 1.0;

load(MAT_fname);

I=find(W > Wmax | W < Wmin); % Remove speed which are not betw. 0.4-0.8 m/s
P(I)=NaN;
J=find(P < Pmin); % Remove data where pressure is less than 2dBar (2m) 
P(J)=NaN;

I = find(isnan(P)==1);

count = 1;
i = 1;
indices = [];
while i < length(I)
    if  I(i+1)-I(i) > 1280 % good profile, at least 20secs.
        indices(count, :)= [count I(i)+1 I(i+1)-1 0]; 
        count=count+1;
        i = i+1;
    else
        i = i+1;
    end    
end

% Quick check the result
% $$$ plot(MTIME, P)
% $$$ hold on
% $$$ for i = 1:size(indices,1)
% $$$     plot(MTIME(indices(i,2):indices(i,3)), P(indices(i,2):indices(i,3)), ...
% $$$          'r');
% $$$ end

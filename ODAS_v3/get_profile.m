%% get_profile
% Extract indices to where an instrument is moving up or down
%%
% <latex>\index{Type A!get\_profile}</latex>
%
%%% Syntax
%   profile = get_profile( P, W, pMin, wMin, direction, minDuration, fs )
%
% * [P] Vector of pressure.
% * [W] Vector of rate of change of pressure - make sure it has no flyers.
% * [pMin] Minimum pressure for finding profiles.
% * [wMin] Minimum rate of change of pressure for finding profiles.
% * [direction] Direction of profile, either 'up' or 'down'.
% * [minDuration] Minimum duration (s) for up/down to be considered a profile.
% * [fs] Sampling rate of P and W in samples per second.
% * []
% * [profile] 2 X N matrix where each column is the start and end indexes for 
%         each profile.  Empty if no profiles detected.
%
%%% Description
% Extract the slow sampling indices to the sections of a profile where the 
% instrument is steadily moving up or down for at least min-duration seconds.
%
% Call this function separately for ascents and descents.
%
% Developed for finding sections in Slocum profiles but it can also be used
% for vertical profilers that have numerous profiles in a single file.
%

% *Version History:*
%
% * 2011-12-26 (RGL) initial
% * 2012-04-30 (WID) comments added to allow for Matlab publishing

function profile = get_profile(P, W, P_min, W_min, direction, min_duration, fs)

min_samples = min_duration*fs; %The minimum number of contiguous samples to be declared a profile.
if strcmpi(direction,'up'); W = -W; end

n = find((P > P_min) & (W >= W_min));
if (length(n) < min_samples)
    profile = [];
    return
end
    
diff_n = diff(n);
diff_n = [diff_n(1); diff_n]; % use diff function to find gaps in profiles

m = find(diff_n > 1);% This locates the breaks bewtween profiles

if isempty(m)
    profile = [n(1) ; n(end)]; % There is only a single profile
else
    profile = zeros(2,length(m)+1);
    profile(:,1) = [n(1); n(m(1)-1)]; % This is the first profile
    
    for index = 2:length(m)
        profile(:,index) = [n(m(index-1));n(m(index)-1)];
    end
    profile(:,end) = [n(m(end));n(end)]; % This is the last profile
end
%Now check the length of each profile
profile_length = profile(2,:) - profile(1,:); % the length of each profile

mm = find(profile_length >= min_samples); % find the profiles longer than min_samples
profile = profile(:,mm); % keep only those profiles

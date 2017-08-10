function view_profiles(profiles, var1, var2)

% function view_profile(profiles, var1, var2)
%
% Where 
% - profiles: the number of profile to look at
%        or: a list of profile from  
%        "ls -1 profile*.mat | sed 's/\.mat//' > profiles"
% - var1: variable to plot in x
% - var2: variable to plot in y
% 
% ex: view_profile(50, 'SBT', 'P')
%          OR
%     view_profile('profiles', 'SBT', 'P')
%
% will plot 50 profiles of SBT against P and will ask you to hit any key
% between 2 plot. This function could be use to have a quick look of
% profiles and may be use to detect profile splitting errors...
%
% Author: Frederic Cyr - 2009/12/18
%
% MODIFICATIONS:
%
% - F. Cyr (july 2011): 
%       Add a new type of input (list). The script can now be used
%       as a number or a list as 1st input.
% ------------------------------------------------------------------- %

% Input a number or a list
if ischar(profiles)==0
    no = profiles;
else
    % load profiles*.mat
    fid = fopen(profiles);
    C = textscan(fid, '%s', 'delimiter', '\n');
    pro_files = char(C{1});

    no = size(pro_files, 1);
end

% Main loop
for profile = 1:no
    
    % Input a number
    if ischar(profiles)==0 % the input is a number
      
        if profile<10
            data_fname = sprintf('profile00%d', profile);
        else
            if profile<100
                data_fname = sprintf('profile0%d', profile);
            else %profile>100
                data_fname = sprintf('profile%d', profile);
            end
        end
   
    else % Input a list
        data_fname = pro_files(profile, :);
    end
    
    % Load and Plot
    load(data_fname);

    figure(1)
    clf
    
    eval(sprintf('plot(%s, %s)', var1, var2))
    set(gca, 'ydir', 'reverse');
    title(sprintf('profile %d',profile))
    disp('Press any key to continue');

    pause
    
end

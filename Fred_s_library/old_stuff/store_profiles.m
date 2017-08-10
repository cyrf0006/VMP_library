function store_profiles(profileList, outfilePrefix, outFolder)

% usage ex in /home/cyrf0006/WINDEX/data_processing:
% >> store_profiles('profile_list/rikiProfiles.list', 'profile', './Robin') 
% to copy all profile*.mat in Robin to give them to him.  
%
% F. Cyr, Jan. 2013
    
fid = fopen(profileList);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

noProfiles = size(files, 1);

%loop on profiles
for i = 1:noProfiles

    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    if i>999
        outfile = [outfilePrefix sprintf('%d.mat', i)];
    elseif i>99
        outfile = [outfilePrefix sprintf('0%d.mat', i)];
    elseif i>9
        outfile = [outfilePrefix sprintf('00%d.mat', i)];
    else
        outfile = [outfilePrefix sprintf('000%d.mat', i)];
    end
    

    
    
    command = sprintf('cp %s %s/%s', fname, outFolder, outfile);
    
     system(command);
    
    
end


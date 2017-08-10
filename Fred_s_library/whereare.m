function [latVec, lonVec] = whereare(pfileList, gps_file)
   
% function [latVec, lonVec] = whereare(pfileList, gps_file)
%
% Similar to whereis.m, but take a profile list as input instead of
% a single file. Returns the coordinates of the profiles (computed from the time at
% which the profile have been taken and the time in GPS file)
%
% usage ex.: [latVec, lonVec] = whereare('epsProfiles.list', 'mergedGPS')
%  in ~/WINDEX/data_processing/Mission_tadoussac_2009
%
% Frederic Cyr - August 2013
%
% ------------------------------------------------------------ %
      
% xtract GPS infos
[latGPS lonGPS nGPS] = xtract_track(gps_file);
lonGPS = lonGPS*-1;
 
% Load profiles
fid = fopen(pfileList);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

noProfiles = size(files,1);

latVec = nan(noProfiles,1);
lonVec = nan(noProfiles,1);
for i = 1:noProfiles

    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
    % find a variable named mtime*
    struc = whos('mtime*');
    varname = struc.name;
    command = (['n = ' varname '(1);']);
    eval(command);
    
    % find lat lon
    [mi I] = min(abs(nGPS-n));
    if mi > 1/24
        latVec(i) = NaN;
        lonVec(i) = NaN;
    else
        latVec(i) = latGPS(I);
        lonVec(i) = lonGPS(I);
    end
end

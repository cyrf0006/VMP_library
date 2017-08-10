function [lat_prof, lon_prof] = whereis(pfile, gps_file)
   
% function [lat_prof, lon_prof] = whereis(pfile, gps_file)
%
% returns the coordinates of a profile (computed from the time at
% which the profile have been taken and the time in GPS file)
%
% usage ex.: [lat_prof, lon_prof] = whereis('profile001.mat', 'gps_2011-02-07')
%
% Frederic Cyr (2011-02-07)
%
% ------------------------------------------------------------ %
      
% xtract GPS infos
[lat_gps lon_gps n_gps] = xtract_track(gps_file);
lon_gps = lon_gps*-1;
    
load(pfile);
% profile coord.
[mi n_ind] = min(abs(n_gps-mtime_eps(1)));
lat_prof = lat_gps(n_ind);
lon_prof = lon_gps(n_ind);




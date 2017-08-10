function [latitude longitude n_GPS] = xtract_track(GPSfile)

% [latitude longitude date time] = function xtract_track(GPSfile)
% tentative de plotter la position a la TCL a chaque Xmin.


track = load(GPSfile); % [LAT LATmin LON LONmin DD MM YYYY HH MM SS]


latitude = track(:, 1) + track(:,2)/60;  
longitude = track(:, 3) + track(:,4)/60; 




n_GPS = datenum(track(:,7), track(:,6), track(:,5), track(:,8), track(:,9), track(:,10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% !!!! 20090721 is different!!!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n_GPS = datenum(track(:,5), track(:,6), track(:,7), track(:,8), track(:,9), track(:,10));

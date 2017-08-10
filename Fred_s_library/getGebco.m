function [latVec lonVec z] = getGebco(ncFile, ncField, lims, varargin)


% usage ex:
% [lat lon z] = getGebco('~/data/GEBCO/GEBCO_08.nc', 'z', [30 65 -80 -35]);
% where lims = [lat0 latE lon0 lonE]
%  
% Will return lat (1xm vector), lon (1xn vector) and z (nxm
% matrix).     
    
% Extract from GEBCO documentation:

% The gridded data are stored in a netCDF data file.
%
% Within the GEBCO_08 Grid and GEBCO_08 SID Grid netCDF files, the
% grids are stored as one dimensional arrays of 2-byte siged
% integer values.

%The complete data sets provide global coverages. Each data set
%consists of 21,600 rows x 43,200 columns, resulting in a total of
%933,120,000 data points. The data start at the Northwest corner of
%the file, i.e. for the global file, position 89°59'45"N,
%179°59'45"W, and are arranged in latitudinal bands of 360 degrees
%x 120 points/degree = 43,200 values. The daa range eastward from
%179°59'45"W to 179°59'45"E. Thus, the first band contains 43,200
%values for 89°59'45"N, then followed by a band of 43,200 values at
%89°59'15"N and so on at 30 arc-second latitude intervals down to
%89°59'45"S.

%The data values are pixel centred registered i.e. they refer to
%elevations at the centre of grid cells. For example, the data
%point at 179°59'45"W, 89°59'45"N represents the data value at the
%centre of a grid cell of dimension 30 arc-seconds of latitude and
%longitude centred on this position i.e. 179°59'30"W - 180°W;
%89°59'30"N - 90°N. 
    
% F. Cyr - January 2013
%
%  modifications:
%   - July 2014 (F. Cyr): Add "varargin" to read 1-min resolution
%   as well (**this modif might not work...).


%% preamble
if isempty(varargin) % default 30-sec res.
    dlatsec = 30; %seconds
    dlonsec = 30; %seconds
else
    dlatsec = varargin{1}; %seconds
    dlonsec = varargin{1}; %seconds
end
lat0 = lims(1);
latE = lims(2);
lon0 = lims(3);
lonE = lims(4);

dlat = dlatsec/(60*60);
dlon = dlonsec/(60*60);

topLat = 90-dlat/2; % Positive North hemisphere
botLat = -90+dlat/2; % Positive East
leftLon = -180+dlon/2; 
rightLon = 180-dlon/2;

latVec = topLat:-dlat:botLat;
lonVec = leftLon:dlon:rightLon;

% Identify desired range
I = find(latVec >= lat0 & latVec <= latE);
J = find(lonVec >= lon0 & lonVec <= lonE);
z = nan(length(I), length(J));

LineSkip = (I(1)-1);
colSkip = J(1)-1;

% initial count
start = LineSkip*length(lonVec)+J(1);
count = length(J);

% loop to load and store data
for i = 1:size(z, 1)
    z(i,:) = nc_varget(ncFile, ncField, start, count);
    start = start+length(lonVec);
end


% Restrict Lat Lon
latVec = latVec(I);
lonVec = lonVec(J);


function [Xfield, xbin, zbin, dist]=getTransect(station_file)
    
    
fid = fopen(station_file);
    
C = textscan(fid, '%s', 'delimiter', '\n');
stat = char(C{1});

no_station = size(stat,1); 

% distance from head (stat16 -> 25)
dist = [698.2 512.3 422.3 350.6 286.9 191.5 128.4 79.8 34.9 0];
dist = round(dist);

% Define regular grid
zbin = 1:5:450;
xbin = 0:1:700;
xbin = dist;

X_stations = nan(length(zbin), no_station);

for i = 1:no_station
    
    file = load(stat(i,:)); % [P X X2.5 X97.2 X2.5_raw X97.5_raw]
                            % or [P X X2.5 X97.2]
    
    % Eliminate bins where number of profile is small (no statistics)
    P = file(:,1);
    X = file(:,2);
    X2p5 = file(:,3);
    X97p5 = file(:,4);
    
    I = find(~isnan(X)==1);
    I1 = find(~isnan(X2p5)==1);
    I2 = find(~isnan(X97p5)==1);
 
    % interp. to grid resolution
    Xitp = interp1(P(I), X(I), zbin);
    X2p5itp = interp1(P(I1), X2p5(I1), zbin);
    X97p5itp = interp1(P(I2), X97p5(I2), zbin);

    X_stations(:,i) = Xitp;

end


[XI, YI] = meshgrid(dist,zbin);
ZI = griddata(xbin,zbin,X_stations,XI,YI, 'linear');
Xfield = ZI;

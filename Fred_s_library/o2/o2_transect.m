clear


%
% 
%
% "ls -1 *.bin > ctd_files"
%
% Author: Frederic Cyr, Oct. 2011
% 
% ------------------------------------------------------------- %
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

% -- Getting infos on profiles -- %
% load file names
fid = fopen('station_file');
C = textscan(fid, '%s', 'delimiter', '\n');
stat = char(C{1});

no_station = size(stat,1); 

% distance from head (stat16 -> 23)
dist = [825 698.2 512.3 422.3 350.6 286.9 191.5 128.4 79.8 34.9 0];
dist = round(dist);

% Define regular grid
zbin = 1:1:450;
xbin = 0:1:825;



T_stations = nan(length(zbin), no_station);
O_stations = nan(length(zbin), no_station);

T_transect = nan(length(zbin), length(xbin));
O_transect = nan(length(zbin), length(xbin));


for i = 1:no_station
    fname = stat(i,:);
    I = find(fname==' ');   
    fname(I) = [];
    file = load(fname); % [P T T2.5 T97.2 T2.5_raw T97.5_raw O ....]

    % Eliminate bins where number of profile is small (no statistics)
    T = file(:,2);
    T2p5 = file(:,3);
    T97p5 = file(:,4);
    I = find(~isnan(T)==1);
    I1 = find(~isnan(T2p5)==1);
    I2 = find(~isnan(T97p5)==1);
    
    if length(I)==length(I1) & length(I)==length(I2)
        T_stations(1:length(I),i) = T(I);
    else        
        T_stations(1:length(I2),i) = T(I2);
    end
        
    O = file(:,7);
    O2p5 = file(:,8);
    O97p5 = file(:,9);
    I = find(~isnan(O)==1);
    I1 = find(~isnan(O2p5)==1);
    I2 = find(~isnan(O97p5)==1);
    
    if length(I)==length(I1) & length(I)==length(I2)
        O_stations(1:length(I),i) = O(I);
    else        
        O_stations(1:length(I2),i) = O(I2);    
    end
    
end

% test griddata

[XI, YI] = meshgrid(xbin,zbin);
ZI = griddata(dist,zbin,O_stations,XI,YI, 'cubic');
O_transect = ZI;
ZI = griddata(dist,zbin,T_stations,XI,YI, 'cubic');
T_transect = ZI;

% $$$ 
% $$$ for i = 1:length(zbin)
% $$$     
% $$$     
% $$$     % T
% $$$     I = find(~isnan(T_stations(i, :)) == 1);
% $$$     if length(I) > 2
% $$$         
% $$$         I1 = find(xbin==min(dist(I)));
% $$$         I2 = find(xbin==max(dist(I)));
% $$$         T_transect(i, I1:I2) = interp1(dist(I), T_stations(i, I), xbin(I1:I2));
% $$$      end
% $$$     
% $$$     
% $$$     % O_2
% $$$     I = find(~isnan(O_stations(i, :)) == 1);
% $$$     if length(I) > 2
% $$$         
% $$$         I1 = find(xbin==min(dist(I)));
% $$$         I2 = find(xbin==max(dist(I)));
% $$$         O_transect(i, I1:I2) = interp1(dist(I), O_stations(i, I), xbin(I1:I2));
% $$$     end
% $$$     
% $$$ end

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

V=[0:5:450];
contourf(xbin, zbin, O_transect, V, 'linestyle', 'none')
hold on
contour(xbin, zbin, O_transect, [62.5 75 150 250 350], 'color', 'k', 'ShowText','on');


for i = 1:length(dist)
    plot([dist(i) dist(i)], [0 450], '--k')
end

title('Oxygen concentration (\mu mol/l)')
xlabel('distance from Tadoussac (km)')
ylabel('depth (m)')
ylim([0 375])
set(gca, 'ydir', 'reverse')


save transect_griddata_cub.mat xbin zbin T_transect O_transect
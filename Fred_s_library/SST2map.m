function SST_map(mapData, LatLon, latLims, lonLims, months)
    
% usage ex: 
% SST_map('sstMap.list', '~/data/front_data/AtlanticLatLon.mat', [45 52], [-70 -55], [1:12])
% SST_map('sstMap.list', '~/data/front_data/AtlanticLatLon.mat', [47.5 49.5], [-70 -67.5], [1:12])
% SST_map('sstMapHDF.list', '~/data/front_data/AtlanticLatLon.mat', [47.5 49.5], [-70 -67.5], [1:12])

% !ls /media/Seagate1TB/Front_data/SST_atlantic/Overlay*.mat > sstMap.list
% 
% Lat-Lon corresponding to maps in the list are provided in
% AtlanticLatLon.mat
 
    
fid = fopen(mapData);
C = textscan(fid, '%s', 'delimiter', '\n');
mapFiles = char(C{1});

% Lat Lon of corresponding SST
load(LatLon);
latVec = lat(:,1);
lonVec = lon(1,:);
Ilat = find(latVec>=latLims(1) & latVec<= latLims(2));
Ilon = find(lonVec>=lonLims(1) & lonVec<= lonLims(2));

% Get time of all maps
command = sprintf(['cat %s | sed -s "s/\\(^.*Overlay_\\)\\([0-9].......\\)\\(.*.mat\\)/\\2/g" > /tmp/tmp'], mapData);
system(command);
theDate = load('/tmp/tmp');
n = datenum(num2str(theDate), 'yyyymmdd');
    
% clear few variables
clear lon lat theDate


for imonth = 1:length(months)
    
    % Preselect dates to avoid extra loading
    I = find(str2num(datestr(n, 5)) == months(imonth));
    ImapFiles = mapFiles(I,:);
    noFiles = size(ImapFiles,1);

    disp(['Proceeding with ' datestr(n(I(imonth)),3) ' - ' sprintf(['%d Files'], noFiles)])
    
    % Load and store all maps for this month
    sstCube = nan(length(Ilat), length(Ilon), noFiles);
    for i = 1:noFiles

        fname = ImapFiles(i, :);
        I = find(fname==' ');   
        fname(I) = [];
        
        load(fname)
        image_overlay=image_overlay(Ilat, Ilon);
        I = find(image_overlay==-4 | image_overlay==-5 |  image_overlay==-2.5);
        image_overlay(I) = NaN;
        
        sstCube(:,:,i) = image_overlay;
    end
       
    SST = nanmean(sstCube,3);
    outfile = ['SSTclim_' datestr(n(I(imonth)),3) '.mat'];
    
    save(outfile, 'SST', 'latVec', 'lonVec');
   
   
   
    
    start = [0 0];
stride = [];
end = [100 50];
data = hdfsd('readdata',dataref,start,stride,end);
    
% $$$    
% $$$    outfile = fname(1:end-4);
% $$$    
% $$$    figure(1)
% $$$    clf
% $$$    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])
% $$$    m_proj('mercator', 'long',[lonLims(1) lonLims(2)],'lat',[latLims(1) latLims(2)]);
% $$$    m_grid('box','fancy')
% $$$    hold on
% $$$    %m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
% $$$    m_pcolor(lonVec, latVec,TMat)     
% $$$    shading interp
% $$$    m_gshhs_h('patch',[.5 .5 .5]); %coastlines
% $$$    hold off
% $$$       
% $$$    colorbar
% $$$    caxis([-2 16])
% $$$    
% $$$    print('-dpng', ['T_' outfile '.png']);  
% $$$    set(gcf, 'renderer', 'painters'); % vectorial figure
% $$$    print('-depsc2', ['T_' outfile '.eps']);                                    
% $$$ 
% $$$    pause(1)

end

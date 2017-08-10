function SST2map_fromHDF(mapData, LatLon, latLims, lonLims, months)
    
% usage ex: 
% SST2map_fromHDF('sstMapHDF.list', '~/data/front_data/AtlanticLatLon.mat', [47.5 49.5], [-70 -67.5], [1:12])
% SST2map_fromHDF('sstMapHDF.list', '~/data/front_data/AtlanticLatLon.mat', [47 50], [-72 -65], [1:12])
% SST2map_fromHDF('sstMapHDF.list', '~/data/front_data/AtlanticLatLon.mat', [45 52], [-72 -55], [1:12])
% !ls /media/Seagate1TB/FrontsDisk1/matlab_atlantic/traitement_INPUT/*/*Cut.hdf > sstMapHDF.list
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
latVec = latVec(Ilat);
lonVec = lonVec(Ilon);

% Get time of all maps
command = sprintf(['cat %s | sed -s "s/\\(.*\\)\\([0-9]\\{8\\}\\)\\(.*.hdf\\)/\\2/g" > /tmp/tmp'], mapData);
system(command);
theDate = load('/tmp/tmp');
n = datenum(num2str(theDate), 'yyyymmdd');
    
% clear few variables
clear lon lat theDate

for imonth = 1:length(months)
    tstart = tic;
    % Preselect dates to avoid extra loading
    I = find(str2num(datestr(n, 5)) == months(imonth));
    ImapFiles = mapFiles(I,:);
    noFiles = size(ImapFiles,1);

    disp(['Proceeding with ' datestr(n(I(imonth)),3) ' - ' sprintf(['%d Files'], noFiles)])
    if imonth>1
        disp([' -> Estimated remaining time: ' sprintf('~%d min.', round(telapsed/60*(length(months)-imonth+1)))]);
    end
    outfile = ['SSTclim_LSLELarge_' datestr(n(I(imonth)),5) '.mat'];
    
    % Load and store all maps for this month
    sstCube = nan(length(Ilat), length(Ilon), noFiles);
    for i = 1:noFiles
        
        fname = ImapFiles(i, :);
        I = find(fname==' ');   
        fname(I) = [];
        
        data = hdfread(fname,'<CutImg>');
        data=data(Ilat, Ilon);
        
        % convert to temperature
        Ipos = find(data>0);
        Ineg = find(data<-1);
        
        dataConverted = nan(size(data));
        
        dataConverted(Ipos) = .15*data(Ipos) - 3;
        dataConverted(Ineg) = .15*(256+data(Ineg)) -3;

        sstCube(:,:,i) = dataConverted;
    end
       
    SST = nanmean(sstCube,3);
    save(outfile, 'SST', 'latVec', 'lonVec');
    telapsed = toc(tstart);

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

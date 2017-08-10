function SST = flag_SST(SST, lonVec, latVec, latLims, lonLims, varargin)
    

    
% usage ex:
% SST = flag_SST(SST, lonVec, latVec, [46.8 50],[-71 -66.8])
% SST = flag_SST(SST, lonVec, latVec, [46.8 50],[-71 -66.8], 'plot')  
    
% used to remove temperature calculated in lakes and Saguenay

if isempty(varargin) == 1
    withplot = 0;
elseif strcmp(varargin{1}, 'plot') == 1
    withplot = 1;
else
    withplot = 0;
end

 % Saguenay
 lon_max1=-69.7;
 lon_min1=-71;
 lat_min1=48.12;
 lat_max1=48.5;
 
 % Northern lakes (1)
 lon_max2=-68.2;
 lon_min2=-71;
 lat_min2=49.25;
 lat_max2=50; 
 % Northern lakes (2)
 lon_max3=-69;
 lon_min3=-71;
 lat_min3=49;
 lat_max3=50;

if withplot
    figure(1)
    clf
    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 12])
    m_proj('mercator', 'long',[lonLims(1) lonLims(2)],'lat',[latLims(1) latLims(2)]);
    m_grid('box','fancy')
    hold on
    %m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
    m_pcolor(lonVec,latVec,SST)     
    shading flat
    hold off
    colorbar
    
    % Sagnuenay
    [lomi, lami] = m_ll2xy(lon_min1,lat_min1);
    [loma, lama] = m_ll2xy(lon_max1,lat_max1);

    hold on
    rectangle('Position', [lomi lami loma-lomi lama-lami], 'linewidth', ...
              2, 'edgecolor', 'r') 
    hold off

    % Northern lakes (1)
    [lomi, lami] = m_ll2xy(lon_min2,lat_min2);
    [loma, lama] = m_ll2xy(lon_max2,lat_max2);

    hold on
    rectangle('Position', [lomi lami loma-lomi lama-lami], 'linewidth', ...
              2, 'edgecolor', 'r') 
    hold off
    
    % Northern lakes (2)
    [lomi, lami] = m_ll2xy(lon_min3,lat_min3);
    [loma, lama] = m_ll2xy(lon_max3,lat_max3);

    hold on
    rectangle('Position', [lomi lami loma-lomi lama-lami], 'linewidth', ...
              2, 'edgecolor', 'r') 
    hold off
end


[lonMat, latMat] = meshgrid(lonVec, latVec); 

I = find(latMat>=lat_min1 & lonMat<=lon_max1);
SST(I) = NaN;
I = find(latMat>=lat_min2 & lonMat<=lon_max2);
SST(I) = NaN;
I = find(latMat>=lat_min3 & lonMat<=lon_max3);
SST(I) = NaN;


if withplot
    figure(2)
    clf
    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 12])
    m_proj('mercator', 'long',[lonLims(1) lonLims(2)],'lat',[latLims(1) latLims(2)]);
    m_grid('box','fancy')
    hold on
    %m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
    m_pcolor(lonVec,latVec,SST)     
    shading flat
    hold off
    colorbar
end
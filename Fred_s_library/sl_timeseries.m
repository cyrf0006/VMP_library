function sl_timeseries(cnvFiles)

%usage ex: >> sleiwex_transect('20130609T1')

fid = fopen(cnvFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

siz = size(files);
noProfile = siz(1); %number of eps_files 

Pbin = 0.5:1:100;
dp = 1;
    
moorLon = -69-54.148/60;
moorLat = 48+12.756/60;

densMat = nan(length(Pbin), noProfile);
for i = 1:noProfile
    
    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    disp(fname)
    
    [lat,lon,gtime,data,names,sensors]=cnv2mat(fname);
    
    
    P = data(:,2);
    T = data(:,3);
    S = data(:,4);
    
    Tbin = nan(size(Pbin));    
    Sbin = nan(size(Pbin));
    for j = 1:length(Pbin)

        I = find(P>=Pbin(j)-dp/2 & P<Pbin(j)+dp/2);
        if isempty(I) ~= 1
            Tbin(j) = nanmean(T(I));
            Sbin(j) = nanmean(S(I));
        end
        
    end
    
    densMat(:,i) = sw_dens(Sbin,Tbin,Pbin);
    TMat(:,i) = Tbin;
    SMat(:,i) = Sbin;

    Lat(i) = lat;
    Lon(i) = lon;
    mtime(i) = datenum(gtime(1));
end

figure(3)
pcolor(Lat, Pbin, SMat)   
shading interp   
set(gca, 'ydir', 'reverse')


x = 0;
for i = 2:length(Lat)
    x =  [x m_lldist([Lon(1) Lon(i)], [Lat(1) Lat(i)])];
end


figure(1)
clf
pcolor(x, Pbin, densMat)   
shading interp
hold on
contour(x, Pbin, densMat, 'linecolor', 'k')
set(gca, 'ydir', 'reverse')
xmoor = m_lldist([Lon(1) moorLon], [Lat(1) moorLat]);
plot([xmoor xmoor], [max(Pbin) min(Pbin)], 'color', [1 1 1], 'linestyle', '--')
   
% Map Figure
figure(5)
clf
%figure dimension
paperwidth = 16;%cm
paperheight = 14;%cm
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; % offset between heigth of colorbar and heigth of figure

lon_min=-69.9167;
lon_max=-69.8333;
lat_min=48.20;
lat_max=48.25;

moorLat = 48+12.642/60;
moorLon = -69-53.719/60;


% $$$ 
% $$$ I=find(lat<lat_max & lat>lat_min);
% $$$ latitude=lat(I);
% $$$ longitude=lon(I);
% $$$ bathy=z(I);
% $$$ 
% $$$ I=find(longitude<lon_max & longitude>lon_min);
% $$$ latitude2=latitude(I);
% $$$ longitude2=longitude(I);
% $$$ bathy2=bathy(I);
% $$$ 
% $$$ clear latitude longitude I bathy lat lon z
% $$$ 
% $$$ % on the new grid
% $$$ y = lat_min:(lat_max-lat_min)/110/5:lat_max;
% $$$ x = lon_min:(lon_max-lon_min)/80/5:lon_max;
% $$$ [X,Y] = meshgrid(x,y);
% $$$ [XI,YI,Z] = griddata(longitude2,latitude2,bathy2,X,Y);
% $$$ 
% $$$ %b=gaussfir(0.001, 5, 3);  %Gaussian filter designer ()
% $$$ b=gaussfir(10);  %Gaussian filter designer ()
% $$$ 
% $$$ Z2=filter2(b,Z);  %Application of the filter

m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
%m_gshhs_h('patch',[.5 .5 .5]); %coastlines
%m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy')
m_gshhs_h('patch',[.5 .5 .5]); %coastlines



for i = 1:length(Lat)
    m_line(Lon(i),Lat(i), 'marker','p','MarkerFaceColor','k', ...
           'markersize',6,'color','k');
end
m_line(moorLon,moorLat, 'marker','o','MarkerFaceColor','r', 'markersize',6,'color','r');

% -- custom colormap! -- %
%c=colormap(gray);
%load bathy_blue2;



function plottrack(logged_lat, logged_lon, n_GPS, no_profile)

% new_version: plottrack(logged_lat, logged_lon, time_GPS, no_profile)
% more comments are needed
% old_version: plottrack(logged_lat, logged_lon, date, time_GPS, T1, T2, timestep)

%%%%%%%%%%%%%%%%%
% build the map %
%%%%%%%%%%%%%%%%%


load '/home/cyrf0006/data/Topo_Gulf_St_Lawrence.mat'


% %%%% LCH %%%%%
% lon_min=-70;
% lon_max=-69;
% lat_min=47.83333;
% lat_max=48.5;
% %%%%%%%%%%%%%%%

% %%%% IML4 %%%%%%
% lon_min=-70;
% lon_max=-67;
% lat_min=48;
% lat_max=49.5;
% %%%%%%%%%%%%%%%%

%%%% IML4 - zoom %%%%%%
lon_min=-68.75;
lon_max=-68.25;
lat_min=48.4;
lat_max=48.8;
%%%%%%%%%%%%%%%%%%%%%%%%%

I=find(lat<lat_max & lat>lat_min);
latitude=lat(I);
longitude=lon(I);
bathy=z(I);

I=find(longitude<lon_max & longitude>lon_min);
latitude2=latitude(I);
longitude2=longitude(I);
bathy2=bathy(I);

clear latitude longitude I bathy lat lon z

% on the new grid
y = lat_min:(lat_max-lat_min)/110/5:lat_max;
x = lon_min:(lon_max-lon_min)/80/5:lon_max;
[X,Y] = meshgrid(x,y);
[XI,YI,Z] = griddata(longitude2,latitude2,bathy2,X,Y);
%ZI = interp2(longitude2, latitude2, bathy2, X, Y, 'linear');

b=gaussfir(0.001, 5, 3);  %Gaussian filter designer ()
Z2=filter2(b,Z);  %Application of the filter


clf
m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
%m_gshhs_h('patch',[.5 .5 .5]); %coastlines
%m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('FontSize', 14,'box','fancy');

V=[0:50:300];

%pack
hold on
[H, H] = m_contourf(X,Y,Z2,V);
%[H, H] = m_contourf(X,Y,Z,V);
set(H, 'LineStyle', 'none')
c = colorbar;
ylabel(c,'profondeur (m)', 'FontSize', 14)

%title('Position de chaque profil', 'FontSize', 14);
%lt = xlabel('Longitude', 'FontSize', 14);
ylabel('Latitude','FontSize', 14 )
m_gshhs_h('patch',[.5 .5 .5]); %coastlines



%%%%%%%%%%%%%%%%
% choose track %
%%%%%%%%%%%%%%%%

%! 1st version

% time_vector = T1:timestep:T2;
% %keyboard
% for t = 1:length(time_vector)
%     ind = dsearchn(time, time_vector(t));
%     track_time = time(ind);
%     track_lat = logged_lat(ind);
%     track_lon = logged_lon(ind);
% 
%     hh=round(track_time);
%     mm=round((track_time-hh)*60);
%     tt=datestr(sprintf('%d:%d', hh, mm),15);
% 
% %    m_line(-track_lon, track_lat,'marker','-','MarkerFaceColor','r','markersize',10,'color','r');
%     m_line(-track_lon, track_lat, 'color','r');
%     if mm==0
%         m_text(-track_lon, track_lat, tt, 'vertical', 'bottom', 'horizontal', 'left', 'color', 'k','FontSize',13);
%         m_line(-track_lon, track_lat,'marker','.','MarkerFaceColor','k','markersize',20,'color','k');
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%
% extract time profile %
%%%%%%%%%%%%%%%%%%%%%%%%
%! 2nd version
for profile = 1:no_profile
    
    % the name of the profile
    if profile<10
        data_fname = sprintf('profile00%d', profile);
    else
        if profile<100
            data_fname = sprintf('profile0%d', profile);
        else %profile>100
            data_fname = sprintf('profile%d', profile);
        end
    end

   
    load(data_fname);
    
    dat = date(1,:);    
    dd = dat(1:2); 
    mm = dat(4:5); 
    yyyy = dat(7:10);
    
    T = [str2num(time(1,1:2)) str2num(time(1, 4:5)) str2num(time(1,7:8))];

    t = (T(1)+T(2)/60+T(3)/60/60)'; %time in decimal round to 1/10 sec.

    n_profile = datenum(str2num(yyyy),str2num(mm),str2num(dd), T(1), T(2), T(3));
    
    %ind = dsearchn(time_GPS, t);
    ind = dsearchn(n_GPS, n_profile);
%    track_time = time_GPS(ind);
    track_lat = logged_lat(ind);
    track_lon = logged_lon(ind);

%    keyboard
    
     %plot a point for each profile
    m_line(-track_lon, track_lat,'marker','.','MarkerFaceColor','k','markersize',10,'color','k');

    % for special plot every hour
    if profile==1
        hour_pt1 = ceil(t);
    end
    if profile==no_profile
        hour_pt2 = floor(t);
    end
    
end

% 
%     %special plot every hour
%         hour_vector = hour_pt1:hour_pt2;
%     
%     for i = 1:length(hour_vector)   
%         
%         ind = dsearchn(n_GPS, hour_vector(i));
%         %track_time = n_GPS(ind);
%         track_lat = logged_lat(ind);
%         track_lon = logged_lon(ind);
% 
%         %hh=round(track_time);
%         %mm=round((track_time-hh)*60);
%         %tt=datestr(sprintf('%d:%d', hh, mm),15); 
%         tt = datestr(n_GPS(ind), 15);
%         
%         m_text(-track_lon, track_lat, tt, 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',13, 'FontWeight', 'bold');
%         m_line(-track_lon, track_lat,'marker','.','MarkerFaceColor','k','markersize',20,'color','k');
%     end
    
    set(gca, 'FontSize', 14)
    
    % Other pts of interest
  
    m_line(-68.522186, 48.482571, 'marker','square','MarkerFaceColor','r','markersize',8,'color','k');
    m_text(-68.522186, 48.482571, 'Rimouski', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',13, 'FontWeight', 'bold');
    m_line(-68.583333, 48.666666,'marker','square','MarkerFaceColor','b','markersize',10,'color','k');
    m_text(-68.583333, 48.666666, 'IML4', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',13, 'FontWeight', 'bold');

    
% %m_text(-69.554157, 48.069431, 'Ile-Rouge','vertical','bottom','FontSize',12);
% m_line(-69.4333, 48.25, 'marker','square','MarkerFaceColor','k','markersize',6,'color','k');
% m_text(-69.4333, 48.25, 'B','vertical','bottom','FontSize',12);
% m_line(-69.5500, 48.1667, 'marker','square','MarkerFaceColor','k','markersize',6,'color','k');
% m_text(-69.5500, 48.1667, 'A','vertical','bottom','FontSize',12);
% hold off

% Location of the stations for the transect












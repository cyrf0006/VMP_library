disp('compute currents...')
bottom_stress


disp('map...')

%figure dimension
paperwidth = 16;%cm
paperheight = 14;%cm

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; % offset between heigth of colorbar and heigth of figure

lon_min=-69.3333;
lon_max=-68;
lat_min=48.33333;
lat_max=49;


% load tides
tide = load('~/WINDEX/data_processing/very_all_profiles/tide_2009-2011.dat');
n = datenum(tide(:,1), tide(:,2),tide(:,3),tide(:,4),tide(:,5),0);
h = tide(:,6);
h_tide = interp1(n, h, ave_time300);

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])

% SUBPLOT 1
s1 = subplot('position',[0.05 0.82 0.9 0.17]);
plot(ave_time300, h_tide); hold on
set(gca, 'xticklabel',[ave_time300(1):1:ave_time300(end)])
datetick('x', 7)
xlim([ave_time300(1) ave_time300(end)])
xlabel('sept. / oct.')

% SUBPLOT 2
s2 = subplot('position',[0.05 0.08 .9 0.65]);

m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('box','fancy')

xlabel('Longitude', 'FontSize', 10)
ylabel('Latitude', 'FontSize', 10)
m_gshhs_h('patch',[.7 .7 .7]); %coastlines                               
hold on
m_line(-68.522186, 48.482571, 'marker','p','MarkerFaceColor','k','markersize',6,'color','k');
m_text(-68.522186, 48.482571, 'Rimouski', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',10, 'FontWeight', 'bold');

disp('loop for arrows...')
% mooring coordinates
c_riki = [-68.583333, 48.666666];
c_N080 = [-68.8091,  48.8089];

for i = 1:length(ave_E300)-8
    
    subplot(s1)
    h0 = plot(ave_time300(i), h_tide(i), '.r');
    
    subplot(s2)
    h1 = m_vec(0.1, c_riki(1), c_riki(2), ave_E300(i), ave_N300(i), 'b');
    h2 = m_vec(0.1, c_N080(1), c_N080(2), ave_E600(i), ave_N600(i), 'b');
    if i == 1
       [hpv5, htv5] = m_vec(0.1, -68.25, 48.4166667, 0.05, 0, ...
               'b', 'key', '5 cm s^{-1}');
       set(htv5,'FontSize',8);
    end
    
    %h1 = m_quiver(c_riki(1), c_riki(2), ave_E300(i), ave_N300(i), 'k', 'linewidth', 1);
    %h2 = m_quiver(c_N080(1), c_N080(2), ave_E600(i), ave_N600(i), 'k', 'linewidth', 1);
    
    if i > 999
        filename = sprintf('/tmp/currents_%d.png', i);
    elseif i >99
        filename = sprintf('/tmp/currents_0%d.png', i);
    elseif i > 9
        filename = sprintf('/tmp/currents_00%d.png', i);
    else
        filename = sprintf('/tmp/currents_000%d.png', i);       
    end
    
    print('-dpng', '-r100', filename);
    
    delete(h0);
    delete(h1);
    delete(h2);
    
end


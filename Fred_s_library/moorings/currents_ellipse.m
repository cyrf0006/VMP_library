spring = 3.5;%m

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

spring = 4;
I = find(h_tide > spring);
Ispring = I(1):I(end);
Ineap = 1:length(h_tide);
Ineap(Ispring) = [];



% $$$ % ----- NEAP TIDE --- %
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
% $$$ 
% $$$ % SUBPLOT 1
% $$$ s1 = subplot('position',[0.05 0.82 0.9 0.17]);
% $$$ plot(ave_time300, h_tide, 'k');
% $$$ set(gca, 'xticklabel',[ave_time300(1):1:ave_time300(end)])
% $$$ datetick('x', 7)
% $$$ xlim([ave_time300(1) ave_time300(end)])
% $$$ xlabel('sept. / oct.')
% $$$ 
% $$$ hold on
% $$$ x1 = ave_time300(Ineap(1):Ispring(1)); 
% $$$ x2 = ave_time300(Ispring);
% $$$ x3 = ave_time300(Ispring(end):Ineap(end));
% $$$ 
% $$$ y1 = ones(size(x1));
% $$$ y2 = ones(size(x2));
% $$$ y3 = ones(size(x3));
% $$$ 
% $$$ 
% $$$ patch([x1 fliplr(x1) x1(1)], [y1*0 y1*5 y1(1)*0], [.85 .85 .85], 'edgecolor', 'none');
% $$$ patch([x2 fliplr(x2) x2(1)], [y2*0 y2*5 y2(1)*0], [1 1 1], 'edgecolor', 'none');
% $$$ patch([x3 fliplr(x3) x3(1)], [y3*0 y3*5 y3(1)*0], [.85 .85 .85], 'edgecolor', 'none');
% $$$ 
% $$$ plot(ave_time300, h_tide, 'k');
% $$$ hold off
% $$$ 
% $$$ 
% $$$ % Neap tide
% $$$ % SUBPLOT 2
% $$$ s2 = subplot('position',[0.05 0.08 .9 0.65]);
% $$$ 
% $$$ m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
% $$$ m_grid('box','fancy')
% $$$ 
% $$$ xlabel('Longitude', 'FontSize', 10)
% $$$ ylabel('Latitude', 'FontSize', 10)
% $$$ m_gshhs_h('patch',[.7 .7 .7]); %coastlines                               
% $$$ hold on
% $$$ m_line(-68.522186, 48.482571, 'marker','p','MarkerFaceColor','k','markersize',6,'color','k');
% $$$ m_text(-68.522186, 48.482571, 'Rimouski', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',10, 'FontWeight', 'bold');
% $$$ 
% $$$ disp('loop for arrows...')
% $$$ % mooring coordinates
% $$$ c_riki = [-68.583333, 48.666666];
% $$$ c_N080 = [-68.8091,  48.8089];
% $$$ 
% $$$ for i = 1:length(Ineap)
% $$$     
% $$$     h1 = m_vec(0.4, c_riki(1), c_riki(2), ave_E300(Ineap(i)), ave_N300(Ineap(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
% $$$     h2 = m_vec(0.4, c_N080(1), c_N080(2), ave_E600(Ineap(i)), ave_N600(Ineap(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
% $$$     
% $$$     if i == 1
% $$$        [hpv5, htv5] = m_vec(0.4, -68.25, 48.4166667, 0.2, 0, ...
% $$$                'k', 'key', '20 cm s^{-1}', 'headwidth', 0, 'shaftwidth', ...
% $$$                             0.25);
% $$$        set(htv5,'FontSize',8);
% $$$        hold on
% $$$     end
% $$$    
% $$$     
% $$$ end
% $$$ print('-dpng', '-r300', 'neaptide.png')
% $$$ print('-depsc', 'neaptide.eps')
% $$$ 
% $$$ 

% $$$ 
% $$$ % ----- SPRING TIDE --- %
% $$$ figure(1)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
% $$$ 
% $$$ % SUBPLOT 1
% $$$ s1 = subplot('position',[0.05 0.82 0.9 0.17]);
% $$$ plot(ave_time300, h_tide, 'k');
% $$$ set(gca, 'xticklabel',[ave_time300(1):1:ave_time300(end)])
% $$$ datetick('x', 7)
% $$$ xlim([ave_time300(1) ave_time300(end)])
% $$$ xlabel('sept. / oct.')
% $$$ 
% $$$ hold on
% $$$ x1 = ave_time300(Ineap(1):Ispring(1)); 
% $$$ x2 = ave_time300(Ispring);
% $$$ x3 = ave_time300(Ispring(end):Ineap(end));
% $$$ 
% $$$ y1 = ones(size(x1));
% $$$ y2 = ones(size(x2));
% $$$ y3 = ones(size(x3));
% $$$ 
% $$$ 
% $$$ patch([x1 fliplr(x1) x1(1)], [y1*0 y1*5 y1(1)*0], [1 1 1], 'edgecolor', 'none');
% $$$ patch([x2 fliplr(x2) x2(1)], [y2*0 y2*5 y2(1)*0], [.85 .85 .85], 'edgecolor', 'none');
% $$$ patch([x3 fliplr(x3) x3(1)], [y3*0 y3*5 y3(1)*0], [1 1 1], 'edgecolor', 'none');
% $$$ 
% $$$ plot(ave_time300, h_tide, 'k');
% $$$ hold off
% $$$ 
% $$$ % SUBPLOT 2
% $$$ s2 = subplot('position',[0.05 0.08 .9 0.65]);
% $$$ 
% $$$ m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
% $$$ m_grid('box','fancy')
% $$$ 
% $$$ xlabel('Longitude', 'FontSize', 10)
% $$$ ylabel('Latitude', 'FontSize', 10)
% $$$ m_gshhs_h('patch',[.7 .7 .7]); %coastlines                               
% $$$ hold on
% $$$ m_line(-68.522186, 48.482571, 'marker','p','MarkerFaceColor','k','markersize',6,'color','k');
% $$$ m_text(-68.522186, 48.482571, 'Rimouski', 'vertical', 'top', 'horizontal', 'left', 'color', 'k','FontSize',10, 'FontWeight', 'bold');
% $$$ 
% $$$ disp('loop for arrows...')
% $$$ % mooring coordinates
% $$$ c_riki = [-68.583333, 48.666666];
% $$$ c_N080 = [-68.8091,  48.8089];
% $$$ 
% $$$ for i = 1:length(Ispring)
% $$$     
% $$$     h1 = m_vec(0.4, c_riki(1), c_riki(2), ave_E300(Ispring(i)), ave_N300(Ispring(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
% $$$     h2 = m_vec(0.4, c_N080(1), c_N080(2), ave_E600(Ispring(i)), ave_N600(Ispring(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
% $$$     
% $$$     if i == 1
% $$$        [hpv5, htv5] = m_vec(0.4, -68.25, 48.4166667, 0.2, 0, ...
% $$$                'k', 'key', '20 cm s^{-1}', 'headwidth', 0, 'shaftwidth', ...
% $$$                             0.25);
% $$$        set(htv5,'FontSize',8);
% $$$        hold on
% $$$     end
% $$$    
% $$$     
% $$$ end
% $$$ 
% $$$ print('-dpng', '-r300', 'springtide.png')
% $$$ print('-depsc', 'springtide.eps')




% ---- Part two: average over tidal cycles ---- %

% Identify high tides
II = [];
count = 1;
for i = 2:length(h_tide)-1
    if h_tide(i) > h_tide(i-1) & h_tide(i) > h_tide(i+1)
        II(count) = i;
        count = count+1;
    end
end


for i = 1:length(II)-1
    E300_cycle(i) = nanmean(ave_E300(II(i):II(i+1)));    
    N300_cycle(i) = nanmean(ave_N300(II(i):II(i+1)));    
    E600_cycle(i) = nanmean(ave_E600(II(i):II(i+1)));
    N600_cycle(i) = nanmean(ave_N600(II(i):II(i+1)));
end

IIspring = find(h_tide(II)>=spring);
IIspring = IIspring(1):IIspring(end);
IIspring(end) = [];
IIneap = [1:IIspring(1)-1 IIspring(end)+1:length(II)-1];

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth-5 paperheight])

subplot(211)
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

for i = 1:length(IIneap)
        
    %h1 = m_vec(0.2, c_riki(1), c_riki(2), E300_cycle(i), N300_cycle(i), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
    %h2 = m_vec(0.2, c_N080(1), c_N080(2), E600_cycle(i), N600_cycle(i), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
    h1 = m_vec(0.2, c_riki(1), c_riki(2), E300_cycle(IIneap(i)), N300_cycle(IIneap(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
    h2 = m_vec(0.2, c_N080(1), c_N080(2), E600_cycle(IIneap(i)), N600_cycle(IIneap(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
    
    if i == 1
       [hpv5, htv5] = m_vec(0.2, -68.25, 48.4166667, 0.05, 0, ...
               'k', 'key', '5 cm s^{-1}', 'headwidth', 0, 'shaftwidth', ...
                            0.25);
       set(htv5,'FontSize',8);
       hold on
    end   
end
% letter I.d.
%m_text(-68.05, 48.35, 'a','verticalalignment', 'bottom', 'fontsize', 10, 'fontweight', 'bold','BackgroundColor',[1 1 1]);

subplot(212)
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

for i = 1:length(E300_cycle(IIspring))
    
    h1 = m_vec(0.2, c_riki(1), c_riki(2), E300_cycle(IIspring(i)), N300_cycle(IIspring(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
    h2 = m_vec(0.2, c_N080(1), c_N080(2), E600_cycle(IIspring(i)), N600_cycle(IIspring(i)), 'k', 'shaftwidth', 0.25, 'headwidth', 0);
    
    if i == 1
       [hpv5, htv5] = m_vec(0.2, -68.25, 48.4166667, 0.05, 0, ...
               'k', 'key', '5 cm s^{-1}', 'headwidth', 0, 'shaftwidth', ...
                            0.25);
       set(htv5,'FontSize',8);
       hold on
    end    
end

% letter I.d.
%m_text(-68.05, 48.35, 'b','verticalalignment', 'bottom', 'fontsize', 10, 'fontweight', 'bold','BackgroundColor',[1 1 1]);



keyboard

print('-dpng', '-r300', 'M2cycle_ave.png')
print('-depsc', 'M2cycle_ave.eps')
function aestus_plot_tide(C, x, z, stamp, figPrefix)

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 3; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.05; % very left of figure
rigs = 0.01; % very right of figure
tops = 0.02; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])

textpos = [x(1)+.8*(x(end)-x(1)) z(1)+.9*(z(end)-z(1))];

for i = 1:size(C, 3)
    subplot(3,3,i)
    %imagesc(x/1000, z, C(:,:,i)) 
    contourf(x/1000, z, C(:,:,i), 100, 'linestyle', 'none');
    text(textpos(1)./1000, textpos(2), sprintf('%dh to HT', stamp(i)), ...
         'color', [0 0 0]);
    set(gca, 'ydir', 'reverse')
    adjust_space
    
    clim = [-.1 .1];
    clim = [24 27];
    hold on
    contour(x/1000, z, C(:,:,i), [25 26 26.5 27 27.5], 'k', 'linewidth', 2);
    hold off
% $$$     if i == 1
% $$$         clim = get(gca, 'clim');
% $$$     else
set(gca, 'clim', clim);
% $$$     end
% $$$     
    hold on
    contour(x/1000, z, C(:,:,i), [1 1], 'color', [1 1 1]);
    hold off
    
    if i == 4
        ylabel('z (m)')
    end
    if i == 7
        xlabel('x (km)')
        h = colorbar('Eastoutside');
        PO = [.36, 0.1 0.0130 0.2];
        set(h, 'pos', PO)
    end
end

% Save figure
%set(gcf, 'renderer', 'painters')
%print('-depsc2', 'CIL_plots.eps') % no colorbar
print('-dpng', '-r300',  [figPrefix '.png']) % no colorbar


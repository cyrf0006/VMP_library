function aestus_plot_tide_rho(rho, U, V, x, z, stamp, figPrefix)

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

for i = 1:size(rho, 3)
    
    xi = min(x):500:max(x);
    zi = [min(z):10:max(z)]';
    [xi, zi] = meshgrid(xi, zi);
    u = interp2(x, z, U(:,:,i), xi, zi);
    v = interp2(x, z, V(:,:,i), xi, zi);

    subplot(3,3,i)
    contourf(x/1000, z, rho(:,:,i), 100, 'linestyle', 'none') 
    text(textpos(1)./1000, textpos(2), sprintf('%dh to HT', stamp(i)), ...
         'color', [0 0 0]);
    hold on
    contour(x/1000, z, rho(:,:,i), [1022:1027], 'k', 'linewidth', 2)
    caxis([1022 1027])
    
    set(gca, 'ydir', 'reverse')
    adjust_space
    
    if i == 4
        ylabel('z (m)')
    end
    if i == 7
        xlabel('x (km)')
        h = colorbar('Eastoutside');
        PO = [.36, 0.1 0.0130 0.2];
        set(h, 'pos', PO)
    end
    
    posi = get(gca, 'pos');
    % quiver
    %quiver(xi/1000, zi, u, v, 0.2,'k', 'filled')
    p1=[xi(:)/1000,zi(:)];
    u = u(:); v = v(:);
    I = find(isnan(u)==1);
    u(I) = []; v(I) = []; p1(I,:) = [];
    %m=abs(u+i*v); 
    scale = .8;
    dar = daspect;
    daspect([dar]);%, set(gca,'color',0.3*[1 1 1])
    arrow3(p1,p1+scale*[u,v],'k', 0.8)
    hold off
    set(gca, 'pos', posi)
    
end

% Save figure
%set(gcf, 'renderer', 'painters')
%print('-depsc2', 'CIL_plots.eps') % no colorbar
print('-dpng', '-r300',  [figPrefix '.png']) % no colorbar


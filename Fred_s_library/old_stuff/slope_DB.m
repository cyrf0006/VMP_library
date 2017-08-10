
% - preamble - %
% output figure
figw = 14; %cm width
figh = 4; %cm height

% Parameters for figure costumization
offset1 = 0.01; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 0.0005; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom
offset5 = 0; % Xtra offset top of figure
offset6 = 0; % Xtra offset right of figure
Ylab_offset = -14; %offset between yaxis and ylabel
Ylab_tightoffset = 0.09; 


% - load data - %
load slope_data

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 figw figh])


plot(year, slope, 'k.')
hold on

% shade area
y1 = (slope'./slope')*mean(slope) + std(slope);
y2 = (slope'./slope')*mean(slope) - std(slope);
patch([year; flipdim(year,1); year(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');

plot(year, slope, 'k.')
plot(year, slope, 'k')
plot(year, (slope./slope)*mean(slope), '--k')

ylim([0.1 0.35])  
xlim([1992 2010]);

%%%%%%%%%%%%%%%%%%%%
% - ICI DANIEL ! - %
%%%%%%%%%%%%%%%%%%%%
% interchanger les commentaires des 2 lignes suivantes
%ylab = ylabel('$\dot{T}_{min} (^{\circ}C mo^{-1})','Interpreter','LaTex', 'fontsize', 10, 'VerticalAlignment', 'top');
ylab = ylabel('T_{min} (^{\circ}C mo^{-1})', 'fontsize', 10, 'VerticalAlignment', 'top');

%set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', 1993:2:2009)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')

Ypos = get(ylab, 'position');
set(ylab, 'position', [Ypos(1)+Ylab_offset Ypos(2) Ypos(3)]);

% -- Output figure costumization -- %
Out=get(gca, 'outer'); %get figure box properties
Pos=get(gca, 'position');
Tig=get(gca, 'tight');

% Increase TightInset for colorbar and Xlabel (forced)
Tig = Tig + [Ylab_tightoffset 0 0 0]; 
Tig = Tig + [0 offset4+0.1 offset5 offset6];
% reduce TightInset to reduce white space around figure
set(gca, 'Position', Out - Tig * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

set(gcf, 'renderer', 'painters'); % vectorial figure                                  
                                  
print('-deps2', 'inter_CIL_slope_xlabel.eps')
print('-dpng', '-r300',  'inter_CIL_slope_xlabel.png')

set(gca, 'xticklabel', [], 'fontsize', 10)
print('-deps2', 'inter_CIL_slope.eps')
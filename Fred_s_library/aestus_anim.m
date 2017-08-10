function aestus_anim(x, t0, tf, outfile, caxes)
        
% function aestus_anim(x, t0, tf, outfile, caxes)
%
% Where 
% - x: the field to plot
% - t0 & tf: beginning and end of time vector desired
% - outfile: variable to plot in y
% - caxes: fixed limits of the colorbar
%    
% This script uses viz.m from D. Bourgault aestus viz_toolbox
%    
% ex: aestus_anim('S', 1, 635, 'Sanim.avi', [25 35])
%    
% Will create Sanim.avi of 635 frame, 4 frame-per-second (this can
% be modifed directly in this script)
%
% Author: Frederic Cyr - 2011/07/01
%
% ------------------------------------------------------------------- %


t = t0:tf;

! rm /tmp/*.png

for i = 1:tf-t0;
    
    figure(1)
    clf
    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 10])

    viz(x, t(i));
    caxis(caxes);
    
    fig = sprintf('/tmp/fig%0.4d', i);
    print('-dpng', '-r300', fig)

end
 
!mencoder mf:///tmp/*.png -mf w=1000:h=1000:fps=4:type=png -ovc copy -oac copy -o output.avi
eval(['! mv output.avi ' outfile])



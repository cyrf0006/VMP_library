function TSG_map(mapData, latLims, lonLims)
    
% usage ex: TSG_map('tsgMap.list', [45 52], [-70 -55])
% usage ex: TSG_map('tsgMap.list', [47.5 49.5], [-70 -67.5])

    
fid = fopen(mapData);
C = textscan(fid, '%s', 'delimiter', '\n');
mapFiles = char(C{1});

noFiles = size(mapFiles,1);


for i = 1:noFiles

   fname = mapFiles(i, :);
   disp(fname)
   I = find(fname==' ');   
   fname(I) = [];
       
   load(fname)
   outfile = fname(1:end-4);
   
   figure(1)
   clf
   set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])
   m_proj('mercator', 'long',[lonLims(1) lonLims(2)],'lat',[latLims(1) latLims(2)]);
   m_grid('box','fancy')
   hold on
   %m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
   m_pcolor(lonVec, latVec,TMat)     
   shading interp
   m_gshhs_h('patch',[.5 .5 .5]); %coastlines
   hold off
      
   colorbar
   caxis([-2 16])
   
   print('-dpng', ['T_' outfile '.png']);  
   set(gcf, 'renderer', 'painters'); % vectorial figure
   print('-depsc2', ['T_' outfile '.eps']);                                    

   pause(1)

   
% $$$    figure(2)
% $$$    clf
% $$$    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])
% $$$    m_proj('mercator', 'long',[lonLims(1) lonLims(2)],'lat',[latLims(1) latLims(2)]);
% $$$    m_grid('box','fancy')
% $$$    %m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
% $$$    m_pcolor(lonVec, latVec,SMat)     
% $$$    shading interp
% $$$    
% $$$    hold on
% $$$    m_gshhs_h('patch',[.5 .5 .5]); %coastlines
% $$$    hold off
% $$$ 
% $$$    print('-dpng', ['TSGclim_S_' outfile '.png']);  
% $$$    set(gcf, 'renderer', 'painters'); % vectorial figure
% $$$    print('-depsc2', ['TSGclim_S_' outfile '.eps']);       
% $$$    
end

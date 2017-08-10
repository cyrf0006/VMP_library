top_vec = 0:200;
bot_vec = 100:300;
% usage: heat_matrix = heat_layers(top_vec, bot_vec)
%
% This script has been written to compute the heat content on
% different layer thickness incompassing the CIL 
%
% author: F. Cyr, oct. 2010
%
% ---------------------------------------------------------------------- %

% some constant
g = 9.81;
rho_0 = 1.025e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
T_core = 1; %degC, threshold for the CIL core
freeze_pt = -1.8;


dz = 1;
P = [1:dz:300]';

Tmat = load('Tclim_matrix')';
Smat = load('Sclim_matrix')';



% heat_matrix
heat_mat = zeros(length(top_vec), length(bot_vec));



for i=1:length(top_vec) 
    i
    A = top_vec(i);
    
    for j=1:length(bot_vec)
        
        B = bot_vec(j);    
        I = find(P>=A & P<=B);
        
        if ~isempty(I)==1
            for t = 1:size(Tmat,2) % loop on month
            
                DENS = sw_dens(Smat(I,t), Tmat(I,t), P(I)); %density for bin into CIL core
                HH(t) = cp*nansum(DENS.*(Tmat(I, t)-freeze_pt))*dz/length(I);
                        
            end
        
            [Pfit,S]=polyfit(1:length(HH), HH, 1);
            heat_mat(i,j)=Pfit(1);
        else
            heat_mat(i,j)=NaN; 
        end
        
    end
end

load heat_matrix

load BR_colormap; 
imagesc(bot_vec, top_vec, heat_matrix)
caxis([-1000 1000])
colormap(mycolormap)
xlabel('bottom limit')
ylabel('top limit')   
title('Heat_content slope')
title('Heat content slope (kJ m^{-3} mo^{-1})')

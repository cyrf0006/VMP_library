clear

tic
load grid

% set up time vector
load('time.out');
mtime = time(:,3);
clear time

% find closest high tide
Tperiod = 12.42;
hightide = min(mtime)+Tperiod*3/4/24*0:Tperiod/24:max(mtime);

for i = 1:length(mtime)
    [Y, closestHT] = min(abs(mtime(i)-hightide));
    time2(i) = (mtime(i)-hightide(closestHT))*24;
end

% $$$ % Ignore first 2 days
I = find(mtime<mtime(1)+2);
mtime(I) = [];
time2(I) = [];
frameNeglected = length(I);

% Tidal average of fields....
dtide = 2;
reg_tide = -6:dtide:6;



U_tide = nan(length(z), length(x), length(reg_tide));
V_tide = nan(length(z), length(x), length(reg_tide));
W_tide = nan(length(z), length(x), length(reg_tide));
T_tide = nan(length(z), length(x), length(reg_tide));
S_tide = nan(length(z), length(x), length(reg_tide));
rho_tide = nan(length(z), length(x), length(reg_tide));

for i = 1: length(reg_tide)
    % I is the frame index that must be loaded and averaged
    I = find(time2 > reg_tide(i) - dtide & time2 < reg_tide(i) + dtide);
    U_cube = nan(length(z), length(x), length(I));
    V_cube = nan(length(z), length(x), length(I));
    W_cube = nan(length(z), length(x), length(I));
    T_cube = nan(length(z), length(x), length(I));
    S_cube = nan(length(z), length(x), length(I));

    
    
    % Only one cycle
    j = 100;
    U = getfield('u', I(j)+frameNeglected);
    U = putNaN(U, flag_u);
    U_tide(:,:, i) = U';
    
    V = getfield('v', I(j)+frameNeglected);
    V = putNaN(V, flag_u);
    V_tide(:,:,i) = V';
    
    W = getfield('w', I(j)+frameNeglected);
    W = putNaN(V, flag_s);
    W_tide(:,:,i) = W';
    
    R = getfield('sig', I(j)+frameNeglected);
    R = putNaN(R, flag_s);
    R_tide(:,:,i) = R';
    
    % Many cycle average
% $$$     for j = 1:length(I)
% $$$ % $$$         T = getfield('T', I(j)+frameNeglected);
% $$$ % $$$         T = putNaN(T, flag_s);
% $$$ % $$$         T_cube(:,:,j) = T';
% $$$ % $$$         
% $$$ % $$$         S = getfield('S', I(j)+frameNeglected);
% $$$ % $$$         S = putNaN(S, flag_s);
% $$$ % $$$         S_cube(:,:,j) = S';
% $$$ % $$$         
% $$$         U = getfield('u', I(j)+frameNeglected);
% $$$         U = putNaN(U, flag_u);
% $$$         U_cube(:,:,j) = U';
% $$$         
% $$$         V = getfield('v', I(j)+frameNeglected);
% $$$         V = putNaN(V, flag_u);
% $$$         V_cube(:,:,j) = V';
% $$$         
% $$$         W = getfield('w', I(j)+frameNeglected);
% $$$         W = putNaN(V, flag_s);
% $$$         W_cube(:,:,j) = W';
% $$$         
% $$$         R = getfield('sig', I(j)+frameNeglected);
% $$$         R = putNaN(R, flag_s);
% $$$         R_cube(:,:,j) = R';
% $$$     end    
% $$$     T_tide(:,:, i) = nanmean(T_cube,3);
% $$$     S_tide(:,:, i) = nanmean(S_cube,3);
% $$$     U_tide(:,:, i) = nanmean(U_cube,3);
% $$$     V_tide(:,:, i) = nanmean(V_cube,3);    
% $$$     W_tide(:,:, i) = nanmean(W_cube,3);
% $$$     R_tide(:,:, i) = nanmean(R_cube,3);

end
    
clear T_cube S_cube U_cube V_cube W_cube


keyboard

[Xi, Zi] = meshgrid(x, z);
for i = 1:size(T_tide, 3)
    rho_tide(:,:,i) = sw_dens(S_tide(:,:,i), T_tide(:,:,i), Zi);
end


aestus_plot_tide(U_tide, x, z, reg_tide, 'Uwholegrid_1cycle')
aestus_plot_tide(V_tide, x, z, reg_tide, 'Vwholegrid_1cycle')
aestus_plot_tide(W_tide, x, z, reg_tide, 'Wwholegrid_1cycle')
aestus_plot_tide(R_tide, x, z, reg_tide, 'Rwholegrid_1cycle')

I = find(x>39000 & x<44000);
J = find(z<140);

aestus_plot_tide(T_tide(J, I, :), x(I), z(J), reg_tide, 'Tslope')

aestus_plot_tide_rho(rho_tide(J, I, :), U_tide(J, I, :), W_tide(J, I, :), x(I), z(J), reg_tide, 'rhoslope')


toc
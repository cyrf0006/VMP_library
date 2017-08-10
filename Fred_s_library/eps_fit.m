clear
load '/home/cyrf0006/WINDEX/data_processing/BMix_study/partI.mat'

%data = load('/home/cyrf0006/WINDEX/data_processing/BMix_study/meanU.mat');
kappa = 0.41; 
%zmax = 50;
zmax = 10;
zmin = 1;



% $$$ nboot = 1000
% $$$ N = size(epsMat, 2)
% $$$ eps_boot_b = nan(length(hab),nboot);
% $$$ for iboot = 1:nboot
% $$$     r = rand(N,1);
% $$$     r = ceil(r*N/1);
% $$$     
% $$$     m = nanmean(log(epsMat(:,r)),2);
% $$$     s2 = nanvar(log(epsMat(:,r)),2);
% $$$     eps_boot_b(:,b) = exp(m+s2/2);
% $$$     
% $$$ % $$$     m = nanmean(log(S2Mat(:,r)),2);
% $$$ % $$$     s2 = nanvar(log(S2Mat(:,r)),2);
% $$$ % $$$     S2_boot_b(:,b) = exp(m+s2/2);
% $$$ end



global Z E 

Z = hab;
m = nanmean(log(epsMat),2);
s2 = nanvar(log(epsMat),2);
E = exp(m+s2/2);

I = find(Z<=zmax & Z>=zmin);
Z = Z(I);
E = E(I);

ustar_0 = 0.005;    

% 1st guess
E0 = ustar_0.^3/kappa./Z;


figure(1)
clf
plot(E, Z)
set(gca, 'ydir','reverse')
hold on
%plot(U0, Z, 'g--')


% Minimization
[asol vals] = fminsearch ( 'fitFun5', ustar_0);

% best fit
E1 = asol.^3/kappa./Z;
plot(E1, Z, 'r--')
hold off
title(sprintf('%d', i))



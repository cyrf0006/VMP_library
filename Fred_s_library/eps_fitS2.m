clear
load '/home/cyrf0006/WINDEX/data_processing/BMix_study/partI.mat'

%data = load('/home/cyrf0006/WINDEX/data_processing/BMix_study/meanU.mat');
kappa = 0.41; 
%zmax = 50;
zmax = 50;
zmin = 4;



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

m = nanmean(log(S2Mat),2);
s2 = nanvar(log(S2Mat),2);
S = exp(m+s2/2);

m = nanmean(log(NMat.^2),2);
s2 = nanvar(log(NMat.^2),2);
N = exp(m+s2/2);

I = find(Z<=zmax & Z>=zmin);
S = S(I);
E = E(I);
N = N(I);
Zz = Z(I);

Ri = N./S;

% $$$ ustar_0 = 0.005;    
% $$$ 
% $$$ % 1st guess
% $$$ E0 = ustar_0.^3/kappa./Z;

% 1st guess
%p = polyfit(S, E,1);
p = polyfit(1./Ri, E,1);

a0 = [p(1)/1.1 1e-8];
%E0 = a0(1).*S+a0(2);
E0 = a0(1)./Ri+a0(2);


figure(1)
clf
plot(E, Zz)
set(gca, 'ydir','normal')
hold on
plot(E0, Zz, 'g--')

Z = Ri;
% Minimization
%[asol vals] = fminsearch ( 'fitFun6', [a0 b0]);
[asol vals] = fminsearch ( 'fitFun7', a0);

% best fit
%E1 = asol(1).*S + asol(2);
E1 = asol(1)./Ri + asol(2);
plot(E1, Zz, 'r--')
plot(.0054^3/.4./Zz, Zz, '--k')
hold off
title(sprintf('%d', i))
set(gca, 'xscale', 'log')



% Minimization from mooring shear
load S2_fromtideshear.mat
m = nanmean(log(S2_bin),2);
s2 = nanvar(log(S2_bin),2);
S2_mooring = exp(m+s2/2);

S2itp = interp1(zhab_S2, S2_mooring, hab);
I = find(hab<=zmax & hab>=zmin);
S2itp = S2itp(I);

I = find(~isnan(S2itp)==1);
Z = N(I)./S2itp(I);
Zz = Zz(I);
E = E(I);

p = polyfit(1./Z, E,1); 
a0 = [p(1) -2.5e-8];

[asol vals] = fminsearch ( 'fitFun7', a0);
E2 = asol(1)./Z + asol(2);
hold on
plot(E2, Zz, 'm--')
hold off
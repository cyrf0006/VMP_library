clear
load('~/WINDEX/data_processing/ctdN080_sept2011_out.mat');
%load('~/WINDEX/data_processing/ctdRiki_sept2011_out.mat');


zmax = 300;
zmin = 1;
zfit = 82;

global Z T 

Z = Pbin;
T = meanT;
I = find(Z<=zmax & Z>=zmin);
Z = Z(I);
T = T(I);

TT = T;
ZZ = Z;
I = find(Z>82);
TT(I) = [];
ZZ(I) = [];
TT = [TT; 4.25; 5];
ZZ = [ZZ; 200; 300];

T = interp1(ZZ,TT,Z);


% $$$ % articial T
% $$$ I = find(Z<=zfit);
% $$$ newT = nan(length(Z),1);
% $$$ newT(I) = T(I);
% $$$ slope = (5-T(I(end)))/(Z(end)-zfit);
% $$$ I = find(Z>zfit);
% $$$ newT(I) = T(zfit)+slope*(Z(I)-zfit);
% $$$ T = newT;

% 1st guess
a0 = [-0.000000000078309; 0.000000076845478; -0.000028887612060; ...
      0.005077796184885; -0.381721371883246; 10.368071802373141];
T0 = a0(1)*Z.^5 + a0(2)*Z.^4 + a0(3)*Z.^3 + a0(4)*Z.^2 + a0(5)*Z + a0(6);

% $$$ % 1st guess
% $$$ a0 = [0.005077796184885; -0.381721371883246; 10.368071802373141];
% $$$ T0 = a0(1)*Z.^2 + a0(2)*Z + a0(3);
% $$$ 

figure(1)
clf
plot(T, Z)
set(gca, 'ydir','reverse')
hold on
plot(T0, Z, 'g--')


% Minimization
[asol vals] = fminsearch ( 'fitFun8', a0);

% best fit
T1 = asol(1)*Z.^5 + asol(2)*Z.^4 + asol(3)*Z.^3 + asol(4)*Z.^2 + asol(5)*Z + asol(6);
%T1 = asol(1)*Z.^2 + asol(2)*Z + asol(3);

plot(T1, Z, 'r--')
hold off
title(sprintf('%d', i))



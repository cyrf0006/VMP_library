function rho_fit(Sfile, Tfile)

% ex: rho_fit('~/PhD/CTD_IML4/BIN/renamed/S_bootclim_09.dat', '~/PhD/CTD_IML4/BIN/renamed/T_bootclim_09.dat')
%Sal = load('~/PhD/CTD_IML4/BIN/renamed/S_bootclim_09.dat');
Sal = load(Sfile);
Tem = load(Tfile);

global Z S

Z = Sal(:,1);
S = Sal(:,2);
T = Tem(:,2);

% For sept2011 profiles
load('/home/cyrf0006/WINDEX/data_processing/ctdRiki_sept2011_out.mat')
%load('/home/cyrf0006/WINDEX/data_processing/ctdN080_sept2011_out.mat')
S = meanS;
T = meanT;
Z = Pbin;

S = sw_dens(S, T, Z); %S is Rho

I = find(isnan(S)==1);
S(I) = [];
Z(I) = [];


a0 = [1000 .01 10 .01 25];

%S1 = 35.*exp(-9./(Z+40));
%S2 = a0(1) + a0(2).*Z + a0(3).*tanh((Z-a0(4))./a0(5));
%S1 = 1050.*exp(-2./(Z+80));

figure(1)
clf
plot(S, Z)
set(gca, 'ydir','reverse')
hold on
%plot(S2, Z, 'k--')

format long
% 1st fit
a0 = [1000 .01 10 .01 25];
asol = fminsearch ( 'fitFun1', a0)

S1 = asol(1) + asol(2).*Z + asol(3).*tanh((Z-asol(4))./asol(5));
plot(S1, Z, 'r--')



% 2nd fit
a0 = [1050 -2 80];
%a0 = [1027.5 -.0924 15.4]; % for sept2011 N080
asol = fminsearch ( 'fitFun2', a0)

S2 = asol(1).*exp(asol(2)./(Z+asol(3)));
plot(S2, Z, 'k--')

legend('obs', 'fit tanh', 'fit exp')

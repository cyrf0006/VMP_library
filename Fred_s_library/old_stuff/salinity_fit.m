function salinity_fit(infile)

%Sal = load('~/PhD/CTD_IML4/BIN/renamed/S_bootclim_09.dat');
Sal = load(infile);
    
global Z S

Z = Sal(:,1);
S = Sal(:,2);

a0 = [25 .01 10 .01 25];

%S1 = 35.*exp(-9./(Z+40));
%S2 = a0(1) + a0(2).*Z + a0(3).*tanh((Z-a0(4))./a0(5));



figure(1)
clf
plot(S, Z)
set(gca, 'ydir','reverse')
hold on
%plot(S2, Z, 'k--')


% 1st fit
a0 = [25 .01 10 .01 25];
asol = fminsearch ( 'fitFun1', a0)

S1 = asol(1) + asol(2).*Z + asol(3).*tanh((Z-asol(4))./asol(5));
plot(S1, Z, 'r--')



% 2nd fit
a0 = [35 -9 40];
asol = fminsearch ( 'fitFun2', a0)

S2 = asol(1).*exp(asol(2)./(Z+asol(3)));
plot(S2, Z, 'k--')

legend('obs', 'fit tanh', 'fit exp')

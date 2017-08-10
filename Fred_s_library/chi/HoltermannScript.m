% Run in /home/cyrf0006/research/NIOZ/IOW/CHI
clear all

load DAT_050p

t1 = uMp.T1_hres;
p = uMp.P_hres;
n = uMnc.time; % in seconds
W = uMp.w;
P = uMp.P;
freqHz = 1/(n(2)-n(1));
% (convention: capital letter is finescale small letter is high_res)


% Viscosity
%DENS = sw_dens(uMp.S, uMp.T1, uMp.P);   <-----------PETER YOU HAVE SALINITY?
%[nu, mu] = viscosity(SBS, uMp.T1, DENS);
nu = 1e-6; 

% round frequencies
fs = round(uMp.fs_fast);
FS = round(uMp.fs_slow);


% Limit the profile:
figure(10)
clf
plot(W, P)
hold on

% Limit on depth 
I  = find(P<=180 & P>=40);
P = P(I);
W = W(I);
I  = find(p<=180 & p>=40);
p = p(I);
t1 = t1(I);

% OR limit on velocity (you chose)


plot(W, P, 'r')
hold off
ylabel('P (dbar)')
xlabel('W (m/s)')
set(gca, 'ydir', 'reverse')

dt1dz = diff(t1)./diff(p');
dt1dz_p = p(1:end-1) + diff(p);

[chi1, eps1, p_chi1, ratio] = chi_holtermann(dt1dz,P, dt1dz_p, W, nu, 2);



figure(11)
clf
subplot(131)
semilogx(chi1, p_chi1)
xlabel('\chi')

subplot(132)
semilogx(eps1, p_chi1)
xlabel('\epsilon')

subplot(133)
semilogx(ratio, p_chi1)
xlabel('ratio')

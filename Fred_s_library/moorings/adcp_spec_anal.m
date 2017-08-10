% Files needed (M_N080_vel, etc. are created in ADCP2mat.m)

clear

theta = 33.5; %deg.


% ----------------------------- Spectral Analysis ---------------------------- %



% Comment one or other
%  --------- for M_N080 --------------- %
load M_N080_vel
% reduce time range:
east_vel(:,1:1000) = [];
east_vel(:,end-1000:end) = [];
north_vel(:,1:1000) = [];
north_vel(:,end-1000:end) = [];

freq = 1/3; % Hz
A = east_vel(40,:); % at certain depth
B = north_vel(40,:);
%A = nanmean(east_vel, 1);
%B = nanmean(north_vel, 1);
A2 = nanmean(east_vel, 1); % barotropic
B2 = nanmean(north_vel, 1);

[valong vcross] = rotate_vecd(A, B, theta);
[valong2 vcross2] = rotate_vecd(A2, B2, theta);
% ----------------------------------- %

%  --------- Temp for M_N080 --------------- %
load Tfield_M_N080.mat
dt = (time_vec_N080(2)-time_vec_N080(1))*86400;
freq_T = 1/dt; % Hz
T = Tgrid(11,:); % 10m hab
T2 = nanmean(Tgrid,1); % barotropic

% ----------------------------------- %
% $$$ %  --------- for M_RIKI --------------- %
% $$$ load M_RIKI_vel
% $$$ freq = 1/10; % Hz
% $$$ A = east_vel(35,:); % at certain depth (baroclinic); 38=pycno
% $$$ B = north_vel(35,:);
% $$$ A2 = nanmean(east_vel, 1); % barotropic
% $$$ B2 = nanmean(north_vel, 1);
% $$$ [valong vcross] = rotate_vecd(A, B, theta);
% $$$ [valong2 vcross2] = rotate_vecd(A2, B2, theta);
% $$$ % ----------------------------------- %

% $$$ % --------- for IML4 ADCP ------------ %
% $$$ load IML4_ADCP2009
% $$$ freq = 1/(30*60);
% $$$ A = ADCP.east_vel(14,:);
% $$$ B = ADCP.north_vel(14,:);
% $$$ time_adcp = ADCP.mtime;
% $$$ 
% $$$ A2 = nanmean(ADCP.east_vel, 1);
% $$$ B2 = nanmean(ADCP.north_vel, 1);
% $$$ [valong vcross] = rotate_vecd(A, B, theta);
% $$$ [valong2 vcross2] = rotate_vecd(A2, B2, theta);
% $$$ % ----------------------------------- %

% -------- Old VERSION was using this --------- %
% Baroclinic 
S = A.^2+B.^2;
I = find(isnan(S)==1);
S(I)=0;       
time_adcp(I)=0;
I = find(isnan(valong)==1);
valong(I)=0;
I = find(isnan(vcross)==1);
vcross(I)=0;

% Barotropic
S2 = A2.^2+B2.^2;
I = find(isnan(S2)==1);
S2(I)=0;       
time_adcp2(I)=0;
I = find(isnan(valong2)==1);
valong2(I)=0;
I = find(isnan(vcross2)==1);
vcross2(I)=0;


% mean removed
% (I checked that detrend and mean removed ~ the same , but lower
% Energy for mean removed)
% $$$ [ps, f] = pwelch(sqrt(S)-nanmean(sqrt(S)), [], [], [], freq*86400);
% $$$ [psa, fa] = pwelch(valong-nanmean(valong), [], [], [], freq*86400);   
% $$$ [psc, fc] = pwelch(vcross-nanmean(vcross), [], [], [], freq*86400);   

% $$$ [ps2, f2] = pwelch(sqrt(S2)-nanmean(sqrt(S2)), [], [], [], freq*86400);   
% $$$ [psa2, fa2] = pwelch(valong2-nanmean(valong2), [], [], [], freq*86400);   
% $$$ [psc2, fc2] = pwelch(vcross2-nanmean(vcross2), [], [], [], freq*86400);   
% $$$ 
% $$$ [ps, f] = pwelch(sqrt(S)-nanmean(sqrt(S)), length(S), [], [], freq);   
% $$$ [psa, fa] = pwelch(valong-nanmean(valong), length(valong), [], [], freq);   
% $$$ [psc, fc] = pwelch(vcross-nanmean(vcross), length(vcross), [], [], freq);   
% $$$ 
% $$$ [ps2, f2] = pwelch(sqrt(S2)-nanmean(sqrt(S2)), length(S2), [], [], freq);   
% $$$ [psa2, fa2] = pwelch(valong2-nanmean(valong2), length(valong2), [], [], freq);   
% $$$ [psc2, fc2] = pwelch(vcross2-nanmean(vcross2), length(vcross2), [], [], freq);   
% $$$ % ------------------------------------------------------- %


%  ---------------- new version (2013-02-17) ------------------- %
Sd = detrend(S);
VAd = detrend(valong);
VCd = detrend(vcross);

nx = max(size(Sd));
na = 10;
w = hanning(floor(nx/na));

%[Pxx,f]=pwelch(Sd,w,0,[],freq*86400);
[ps, f] = pwelch(Sd, w, 0, [], freq*86400); 
[psVA, fVA] = pwelch(VAd, w, 0, [], freq*86400); 
[psVC, fVC] = pwelch(VCd, w, 0, [], freq*86400); 

Td = detrend(T);
Td2 = detrend(T2);

nx = max(size(T));
na = 10;
w = hanning(floor(nx/na));
[pst, ft] = pwelch(Td, w, 0, [], freq_T*86400); 
[pst2, ft2] = pwelch(Td2, w, 0, [], freq_T*86400); 
% -------------------------------------------------------------- %



% $$$ % ------- Along_shore -------- %
% $$$ figure(1)
% $$$ clf
% $$$ %loglog(1./(fa2.*86400), psa2, 'k')  
% $$$ loglog(fa2, psa2, 'k')  
% $$$ title('U spectrum')  
% $$$ ylabel('PSD (m s^{-1} cpd^{-1})')
% $$$ xlabel('\sigma (d^{-1})')
% $$$ %xlim([1e-2 1e1])
% $$$ %ylim([1e-3 1e5])
% $$$ xlim([1e-1 2e1])
% $$$ ylim([1e-6 1e0])
% $$$ 
% $$$ hold on
% $$$ %loglog(1./(fa.*86400), psa, 'r')  
% $$$ legend('baroclinic (58m)', 'barotropic (45-145m)')
% $$$ %legend('baroclinic (58m)', 'barotropic (6-98m)')
% $$$ plot_harmonics
% $$$ hold off
% $$$ % ----------------------------- %
% $$$ 
% $$$ % ------- Cross_shore -------- %
% $$$ figure(2)
% $$$ %loglog(1./(fc.*86400), psc)  
% $$$ loglog(f, ps)  
% $$$ title('V spectrum')  
% $$$ ylabel('PSD (m s^{-1} cpd^{-1})')
% $$$ xlabel('period (day)')
% $$$ %xlim([1e-2 1e1])
% $$$ %ylim([1e-3 1e5])
% $$$ xlim([1e-1 1e2])
% $$$ ylim([1e-6 1e0])
% $$$ 
% $$$ hold on
% $$$ plot_harmonics
% $$$ hold off
% $$$ % ----------------------------- %
% $$$ 
% $$$ % ------- Mean Shear -------- %
% $$$ figure(3)
% $$$ %loglog(1./(f2.*86400), ps2, 'k')  
% $$$ loglog(f2, ps2, 'k')  
% $$$ title('E_k spectrum')  
% $$$ ylabel('PSD (m^2 s^{-2} cpd^{-1})')
% $$$ xlabel('\sigma (d^{-1})')
% $$$ %xlim([1e-2 1e1])
% $$$ %ylim([1e-3 1e5])
% $$$ xlim([1e-1 2e1])
% $$$ ylim([1e-6 1e0])
% $$$ 
% $$$ hold on
% $$$ %loglog(1./(f.*86400), ps, 'r')  
% $$$ legend('barotropic (55-145m)', 'baroclinic (58m)')
% $$$ %legend('barotropic (6-98m)', 'baroclinic (58m)')
% $$$ plot_harmonics
% $$$ hold off
% $$$ % ----------------------------- %

% ------- Temperature -------- %
I = find(f < 100);
J = find(ft < 100);

figure(4)
clf
%loglog(1./(fc.*86400), psc)  
loglog(ft(J), pst(J), 'm', 'linewidth', 2)  
hold on
loglog(ft2(J), pst2(J), 'k', 'linewidth', 2)  

title('T spectrum')  
ylabel('PSD (m s^{-1} cpd^{-1})')
xlabel('f (day^{-1})')
xlim([1e-1 1e2])
ylim([1e-8 1e0])


loglog(f(I), ps(I), 'r')
loglog(fVA(I), psVA(I), 'g')
loglog(fVC(I), psVC(I), 'b')

plot_harmonics2
hold off
legend('T', 'T_{10m}', 'U^2+V^2', 'U', 'V')

% ----------------------------- %

% ------- Temperature -------- %
I = find(f < 100);
J = find(ft < 100);

figure(4)
clf
%loglog(1./(fc.*86400), psc)  
loglog(ft(J), pst(J), 'm', 'linewidth', 2)  
hold on
loglog(ft2(J), pst2(J), 'k', 'linewidth', 2)  

title('T spectrum')  
ylabel('PSD')
xlabel('f (day^{-1})')
xlim([1e-1 1e2])
ylim([1e-8 1e0])


loglog(f(I), ps(I), 'r')
loglog(fVA(I), psVA(I), 'g')
loglog(fVC(I), psVC(I), 'b')

plot_harmonics2
hold off
legend('T', 'T_{10m}', 'U^2+V^2', 'U', 'V')
% ----------------------------- %


keyboard

% PSD Contour plot 

for i = 1:50
    A = east_vel(i,:); % at certain depth
    B = north_vel(i,:);
    S = A.^2+B.^2;
    I = find(isnan(S)==1);
    S(I)=0;     
    [ps, f] = pwelch(sqrt(S), [], [], [], freq);       
    mat_ps(i, :) = ps'; 
end


imagesc(f, 1:50, log10(mat_ps))
colorbar
caxis([-1 ,1])   
set(gca, 'xscale', 'log')






% ---------------------------------------------------------------------------- %

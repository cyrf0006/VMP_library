% Artificail_K %
%fits a profile on epsilon...



% multiplying factor
factor = 1;

% % some constants
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.035e3;%kg/m^3
GAMMA = 0.2; %mixing efficiency
freeze_pt = -1.8; 
cp = 3.99; %Kj/Kg/K
CIL_def = 1;

eps_filename = 'EPSILON_mean';
max_value = 4e-7;
max_depth=170;

Tprof = load('T_climato_04.dat');
Sprof = load('S_climato_04.dat');

% -- initial profiles -- %
P = Tprof(:,1);
T = Tprof(:,2);
S = Sprof(:,2);
DENS = sort(sw_dens(S, T, P));

% $$$ % -- Mean epsilon / 10 pts window runmean -- %
% $$$ load eps_ave;
% $$$ P_eps = P_eps';
% $$$ eps_ave10 = eps_ave10';
% $$$ EPS = interp1(P_eps, eps_ave10 , P(1:length(P)-1));

% -- Mean epsilon / 90% -- %
load(eps_filename);
P_eps = P_eps(1:350)';
eps_ave = eps_mean5_95p100(1:350)';
EPS = interp1(P_eps, eps_ave , P(1:length(P)-1));


if ~isempty(find(isnan(EPS)==1)) %remove NaN
    I = find(~isnan(EPS)==1);
    EPS(1:I(1))=EPS(I(1));
    EPS(I(end):end)=EPS(I(end));
end


% -- Brunt-Vaisala -- %
dz = P(2)-P(1);
drho = diff(DENS);
N2 = g/rho_0*drho/dz;


% -- K computation if epsilon method -- %
K = GAMMA.*EPS./N2; %K from epsilon


% Approx epsilon and K profile    
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 16 10])

subplot(1,3,1)
%semilogx(EPS, P(1:end-1))
semilogx(EPS(1:max_depth), P(1:max_depth), 'k')
set(gca, 'ydir', 'reverse')
ylabel('depth(m)')
title({'\epsilon (W/Kg)'})
set(gca, 'xtick', [1e-9 1e-8 1e-7 1e-6 1e-5])
%set(gca, 'xticklabel', [{'' '10^{-8}' '' '10^{-6}' ''}])
set(gca, 'xgrid', 'on', 'ygrid', 'on', 'xminorgrid', 'off')

e_app = (exp(-P./90)./300000)./P; 
%e_app = fit_eps(P_eps, eps_ave, P(1:end), max_value);
I = find(e_app>4e-7);
e_app(I) = 4e-7;

% multiplication by factor
e_app = e_app*factor;

hold on
semilogx(e_app(1:max_depth), P(1:max_depth), 'k--');

smooth_e = runmean(EPS, 5);
semilogx(smooth_e(1:max_depth), P(1:max_depth), 'r--');

hold off

subplot(1,3,2)
%semilogx(N2, P(1:end-1))
semilogx(N2(1:max_depth), P(1:max_depth), 'k')
set(gca, 'ydir', 'reverse')
set(gca, 'xtick', [1e-5 1e-4 1e-3 1e-2], 'xlim', [1e-5 1e-2])
%set(gca, 'xticklabel', [10^{-4} 10^{-3}])
set(gca, 'yticklabel', [], 'xgrid', 'on', 'ygrid', 'on', 'xminorgrid', 'off')
title({'N^2 (s^{-2})'})

subplot(1,3,3)
%semilogx(K, P(1:end-1))
semilogx(K(1:max_depth), P(1:max_depth), 'k')
set(gca, 'ydir', 'reverse')
set(gca, 'xtick', [1e-7 1e-6 1e-5 1e-4], 'xlim', [1e-7 1e-4])
%set(gca, 'xticklabel', [10^{-6} 10^{-5}])
set(gca, 'yticklabel', [], 'xgrid', 'on', 'ygrid', 'on', 'xminorgrid', 'off')

k_app = (exp(-P./500)./1000)./P; 
I = find(k_app>5e-5);
k_app(I) = 5e-5;

% $$$ hold on
% $$$ semilogx(k_app(1:max_depth), P(1:max_depth), 'k--');
% $$$ hold off
%title('Diffusivity')

title('K_{\rho} (m^2/s)')
%legend('obs', 'approx', 'location', 'northwest')

%k_app(end)

print('-deps2', 'approx_e_subplot.eps')
dlmwrite('approx_e.dat', [P e_app],'delimiter',' ','precision',6)
%dlmwrite('approx_e.dat', [P smooth_e],'delimiter',' ','precision',6)

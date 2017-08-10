
fname_in =  'CFL_03aDEC2007_000.mat';
fname_out = 'CFL_03aDEC2007_000_redux.mat';

year = 2007; month = 12; day = 03;

% Threshold for falling speed
wmin=0.4; 
  
% Mixing efficiency
Gamma = 0.2;
g = 9.81;

% Parameters for despiking
thresh_spike = 5; smooth_spike = 0.5; N_spike = 5;

% Parameters for calculating dissipation
nu0 = 1.d-6;
dp = 5.0;     % Bin interval
pmin = 8.0;   % Minimum depth
Options.CleanSpec = 0;
Options.Nfft = 512; %256; %512;
Options.AveOverlap = dp;
Options.Kprb = 48;
Options.Fs = 512;
Options.Kmin = 1/dp;      
Options.Kmax = 200;      


load(fname_in);
  
% Despike the shear signals
[sh1 spike_index1] = despike(sh1, thresh_spike, smooth_spike, fs, N_spike);
[sh2 spike_index2] = despike(sh2, thresh_spike, smooth_spike, fs, N_spike);

 ax = ang2acc(pitch); ay = ang2acc(roll);
 J = find(pitch > 10 | roll > 10);
 sh1(J)= NaN; sh2(J) = NaN;
    
% Calculate viscosity
[NU,MU] = viscosity(SBS,SBT,sigmat(SBT,SBS)+1000);
nu = interp1(P,NU,p);
inan = find(isnan(nu)==1);
nu(inan) = nu0;
    
% ** CALCULATE DISSIPATION RATES **
Options.AveLen = pmin:dp:max(p);
[E1, Spec1] = epsprofile_db(sh1, p, w, nu, Options, [ax, ay, az]);
[E2, Spec2] = epsprofile_db(sh2, p, w, nu, Options, [ax, ay, az]);

    
I = find(E1.W >= wmin);
eps1 = E1.e(I);
eps2 = E2.e(I);
eps = (eps1 + eps2)/2;

p_eps = E1.z(I); 
visc = 0.5*(E1.nu(I)+E1.nu(I));
    
mtime_eps  = interp1(P,(MTIME-mean(MTIME))/std(MTIME),p_eps,'linear');
mtime_eps  = mtime_eps*std(MTIME)+mean(MTIME);
    
  
save(fname_out,'mtime_eps','p_eps','eps','eps1','eps2','visc');
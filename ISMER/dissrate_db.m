function [Results,Spectra,Options] = dissrate(dudz,W,nu,Options,accel)
% DISSRATE Dissipation rate of turbulent kinetic energy.
%
% Calculates the dissipation rate of TKE for the shear signal dudz by iterative
% integration of the dissipation shear spectrum. The algorithm described in more
% detail below. W and nu are the sinking velocity of the instrument and the 
% kinematic viscosity, respectively. They can be specified as scalars or vectors.
% In the latter case, they will be reduced to their mean values during the 
% computation. If Options.CleanSpec==1 and accelerometer data are included
% (vector or matrix, one signal per column), vibrational contamination is
% removed from the shear spectra.  The results of the computation are
% returned in the structures E  and S (see output parameter section below).
%
% Usage: 
%   [E,S] = dissrate(dudz,W,nu)
%   [E,S] = dissrate(dudz,W,nu,Options)
%   [E,S] = dissrate(dudz,W,nu,Options,accel)
% The default behaviour of the function can be changed by specifying one or more
% fields in the Options structure. Note that the fields are case sensitive. The 
% option fields and their default values are:
%
% Field:           Description:                                  
% Options.Nfft     Length of the FFT sections                    
% Options.Overlap  Amount of overlap of FFT sections             
% Options.Kmin     Mininum integration wavenumber                
% Options.Kmax     Maximum integration wavenumber                
% Options.Kprb     Response limit of the shear probe             
% Options.H0       Transfer function of the differentiator board 
% Options.F0       Frequency bins of H0                          
% Options.Fs       Sampling frequency                            
% Options.CleanSpec Flag for removal of vibrational contamination
% Options.KCINITIAL Initial value used for the upper limit of integration
%
% Options.Kprb is the assumed averaging wavelength of the probe. If this
% field is set to 0 no spectral correction for shear probe averaging effects will 
% be performed. Options.H0 contains the transfer function of the shear probe 
% differentiator (see XFRCORR for details). If Options.H0 is empty no correction
% for the differentiator characteristic is performed.
%
% [E,S,Options] = dissrate(...)
% Returns the complete set of options used in the calculation. This is useful 
% meta data. 
%
% Input parameters:
%  dudz     shear signal                [1/s]
%  nu       viscosity                   [m^2/s]
%  W        mean fall speed             [m/s]
%  Options  Option parameters (structure)
%  accel    accelerometer signals       [m/s^2]
%
% Output parameters:
%  E.e      Dissipation rate                         [W/kg]
%  E.ks     Kolmogorov wavenumber                    [cpm]
%  E.kc     Cutoff wavenumber of the integration     [cpm]
%  E.nu     Viscosity used in the calculation        [m^2/s]
%  E.W      Sinking velocity used in the calculation [m/s]
%
%  S.P      Disspation spectrum                      [s^(-2)/cpm]
%  S.k      Wavenumber vector for P                  [cpm]
%  S.Pn     Nasmyth spectrum for E.e                 [s^(-2)/cpm]
%
%  Options  Computation options
%
% Algorithm: 
%   1) Compute the power spectrum P of the shear signal using 
%      Welch's averaged periodogram method (nfft-point spectra).
%      Optionally, remove the portion of the signal that is coherent with
%      the accelerometer signal (using clean_shear_spec.m)
%   2) Apply spectral correction for differentiator response (optional).
%   3) Apply spectral correction for shear probe averaging (optional).
%   4) Set the intial upper integration limit to kc=15 cpm.
%   5) Compute the variance by integrating P between kmin and kc. 
%   6) Compute dissipation rate and Kolmogorov wavenumber ks based on 
%      the variance computed in step 5.
%   6) Increase kc by kinc.
%   7) Re-iterate steps 5 to 7 until a stop criterion occurs.
%  10) Adjust the disspation rate by extrapolating the measured spectrum 
%      above kc and below kmin using the Nasmyth spectrum scaled by the 
%      last computed epsilon value.
%
%  Stop criteria: a) kc is within 2 cpm of the last computed ks.
%                 b) kc is greater than kmax.
%                 c) There are more than 50  iterations.
% The wavenumber increment kinc is determined based on the first estimate of the 
% dissipation rate (step 5). If the initial disspation rate estimate is 
%         e < 1e-9 kinc = 2 cpm
% 1e-9 <= e < 1e-8 kinc = 3 cpm 
% 1e-8 <= e < 1e-6 kinc = 5 cpm
% 1e-6 <= e        kinc = 10 cpm.
%
% Sytnax summary:
%   [E,S] = dissrate(dudz,W,nu)
%   [E,S] = dissrate(dudz,W,nu,Options)
%   [E,S] = dissrate(dudz,W,nu,Options,accel)
%   [E,S,Options] = dissrate(...)

% (C) 2002 Rockland Oceanographic Services Inc.
% Author: Fabian Wolk
% Revision: 2002/08/01
% Isabelle Gaboury, 17 Mar. 2005 - added options to remove vibrational
%   contamination
% IG, 20 Apr. 2005 - modified definition of Kmax (now use a constant value)

% Error checking:
%
warning('on')
if nargin < 4 | isempty(Options)
   Options.Fs = parameter_list('default','Fs');
   Options.Nfft = 2*Options.Fs;
   Options.Overlap = 0;
   Options.Kmin = 0;
   Options.Kmax = 50;
   Options.Kprb= 50;
   Options.H0 = [];
   Options.F0 = [];
   Options.CleanSpec = 0;
   Options.KCINITIAL = 15;
else
   Options = checkOptions(Options); % local function
end
if nargin<5
    accel = [];
    if isfield(Options,'CleanSpec') & Options.CleanSpec==1
        disp('Unable to remove vibrational contamination');
        Options.CleanSpec=0;
    end
end

% Reduce sinking velocity and viscosity are vectors
% to their mean values:
nu = mean(nu);
W = abs(mean(W));

% Some hard coded parameters:
%
MAXITER   = 1000;        % maximum number of iterations
DELTA     = 2;           % maximum difference in wavenumber between iterations
inv2pi    = 1/2/pi;      % a constant
nu7p5     = 7.5*nu;      % a constant
SPORDER   = 2;           % the order of the response correction filter for the shear probe

% Compute the spectrum and turn it into a wavenumber spectrum:
%
Fs = Options.Fs;
if Options.CleanSpec==1
    [P,junk1,junk2,junk3,f] = clean_shear_spec(accel-repmat(mean(accel),size(accel,1),1),dudz,Options);
    P = squeeze(P(1,1,:));
    if max(abs(imag(P)./real(P)))>1e-5    % Cleaning may leave some imaginary values; just ignore if they're small enough
        error('Computed spectrum is complex. Check input data consistency.');
    else
        P = real(P);
    end
else
    [P,f] = psd(dudz,Options.Nfft,Fs,Options.Nfft,Options.Overlap,'linear');
    P = (2*P/Fs); % Scale the spectrum so that the integral P(f)*df is equal to the variance
end
if any(isnan(P))
   warning('(DISSRATE.M) Computed spectrum contains NaN''s. Check input data consistency.');
elseif all(P==0)
   error('Computed spectrum is zero. Check input data consistency.');
elseif ~isreal(P)
   error('Computed spectrum is complex. Check input data consistency.');
elseif any(isinf(P))
   warning('(DISSRATE.M) Computed spectrum contains infinite values. Check input data consistency');
end

% Apply spectral correction for analog circuit.
% This must be done in frequency space.
if ~isempty(Options.H0)
   P = xfrcorr(P,f,Options.H0,Options.F0);
end

% Now turn the spectrum into a wavenumber spectrum:
%
P = P*W; 
k = f/W;

% Apply the spectral correction for the shear probe
% wavenumber cutoff.
if Options.Kprb > 0
   P = respcorr(P,k,Options.Kprb,SPORDER);
end

%% Original method provided by Rockland

% Zero-th iteration. Here we find a first approximation to
% the dissipation rate and the Kolmogorov wavenumber. This
% gives us an idea of "where we are on the scale".
%kc = Options.KCINITIAL; 
%kmin = Options.Kmin;
%e = nu7p5*integral(P,k,kmin,kc);
%ks = inv2pi*(e/nu^3)^0.25;
%
%%% Initialize the iterative integration:
%%%
%kmax = Options.Kmax;
%
%%% A function that finds kmax by fitting a polynomial.
%%% kmax = findkmax(P,k)
%if kc >= kmax
%   kc = kmax;
%end
%
%%% Determine step size:
%%%
%if e < 1e-9
%   kinc = 2;
%elseif e < 1e-8
%   kinc = 3;
%elseif e < 1e-6
%   kinc = 5;
%else
%   kinc = 10;
%end
%kinc=1.0;
%
%%% Main iteration loop:
%%%
%for n = 1:MAXITER
%   if (abs(ks-kc) < DELTA | kc >= kmax), % stop criteria
%      break;
%   end
%   kc = kc+kinc;
%   e = nu7p5*integral(P,k,kmin,kc);
%   ks = inv2pi*(e/nu^3)^0.25;
%   %fprintf('Iteration %d: e=%0.2e ks=%0.2f kc=%6.2f abs(ks-kc)=%7.3f\n', n, e, ks, kc, abs(ks-kc))   
%end
%%%fprintf('\n');
%% New Method by D. Bourgault
kmin = Options.Kmin;

% Compute e10 for the determination of kc = f(e10)
e10 = nu7p5*integral(P,k,kmin,10);

% Conmpute a safe kmax
kmax = 180 - 165.*exp(-e10/(5e-7));

e = nu7p5*integral(P,k,kmin,kmax);
ks = inv2pi*(e/nu^3)^0.25;

%%
% Fit the observations to the Nasmyth spectrum
options=optimset('TolX',1d-12,'TolFun',1.d-12);
if ks < kmax
  ii = find(k <= ks & k >= kmin);
else
  ii = find(k <= kmax & k >= kmin);
end
en = fminsearch(@spectral_fit,e,options,nu,k(ii),P(ii));
Pn = PSD_Nasmyth(en,nu,k);
ksN = inv2pi*(en/nu^3)^0.25;

% Calculate the error between theory and measurements
%misfit=mean(((log(Pn(ii)+1.d-12)-log(P(ii)+1.d-12))./log(Pn(ii)+1.d-12)).^2);
%misfit=sqrt(mean(((log(Pn(ii)+1.d-12)-log(P(ii)+1.d-12))./log(Pn(ii)+1.d-12)).^2));
%misfit = mad(P(ii)./Pn(ii));
misfit=sum((log(Pn(ii)+1.d-12)-log(P(ii)+1.d-12)).^2);

%% Extrapolation using Nasmyth
% Extrapolate the spectrum using the Nasmyth fitted function if needed
% Bourgault's method
if kmax < ksN
  %e = e + nu7p5*integral(Pn,k,kc,ks);
  e = e + nu7p5*integral(Pn,k,kmax,ksN);
end
if kmin > 0 
   e = e + nu7p5*integral(Pn,k,0,kmin);
end
ks = inv2pi*(e/nu^3)^0.25;
kc = kmax;

% Extrapolate the spectrum using the Nasmyth fitted function
% Original method provided by Rockland
% Extrapolate the spectrum up to ks:
%
% Pn = nasmyth(e,nu,k); 
% if kc < ks
%    e = e+nu7p5*integral(Pn,k,kc,ks);
%    ks = inv2pi*(e/nu^3)^0.25;
%    Pn = nasmyth(e,nu,k); 
% end
% 
% % Extrapolate the spectrum below kmin:
% %
% if kmin>0 
%    e = e+nu7p5*integral(Pn,k,1,kmin);
%    ks = inv2pi*(e/nu^3)^0.25;
%    Pn = nasmyth(e,nu,k); 
% end
%%
if e > 0.1
%if e > 0.01
    %Results.e = e.*NaN;
    %Results.e = NaN;
    e=NaN;
end

visual_inspection = 0;
misfit_max = 1.2;
if visual_inspection==1
%if misfit >= misfit_max | e > 1.d-3
%if e > 1.d-3 & e < 0.1
  misfit;
  loglog(k,P);
  hold on;
  loglog(k,Pn,'m');
  xlabel('{\it k} (cpm)');
  ylabel('{\it k} (s^{-2} cpm^{-1})');
  loglog([ks ks],[min(P) max(P)],'k');
  loglog([kc kc],[min(P) max(P)],'r');
  %loglog([kmax kmax],[min(P) max(P)],'r');
  
  %legend('Spectrum','Nasmyth','Kolmo','k-cutoff');
  title(['Kolmogorov: ',num2str(ks),'\epsilon = ',num2str(e),' \epsilon_{fit} = ',num2str(en),' mfit =', num2str(misfit)]);

  answer = input('Is this ok? (y/n): ','s');
  %if answer == 'y'
  %  load cutoff.mat;
  %  iii=length(kc0);
  %  [kc0(iii+1) dummy] = ginput(1)
  %  dudz2(iii+1)=mean(dudz.^2)
  %  maxdudz(iii+1)=max(abs(dudz));
  %  WW(iii+1)=W;
  %  E10(iii+1)=e10;
  %  fprintf('%d %d',kc0(iii+1), dudz2(iii+1));
  %  fprintf('\n');
  %  save cutoff.mat kc0 dudz2 maxdudz WW E10
  %end
  clf;

  if answer == 'y'
    Results.e = e;
  elseif answer == 'n'
    Results.e = e.*NaN;
  end
else
  Results.e = e;
end

Results.ks = ks;
Results.kc = kc;
Results.nu = nu;
Results.W = W;

Spectra.k = k;
Spectra.P = P;
Spectra.Pn = Pn;

return;
%%% end of main function 'dissrate'

function Options = checkOptions(Options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function for dissrate. Checks the Options structure for consistency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(Options)
   error('Options parameter is not a valid structure.');
end

% Make sure that all options fields exist. The user may 
% only pass a subset of all possible options.
if ~isfield(Options,'Fs')                  Options.Fs = parameter_list('default','Fs'); end;
if ~isfield(Options,'Nfft')              Options.Nfft = 2*Options.Fs; end;
if ~isfield(Options,'Overlap')        Options.Overlap = 0; end;
if ~isfield(Options,'Kmin')               Options.Kmin = 0; end;
if ~isfield(Options,'Kmax')              Options.Kmax = 50; end;
if ~isfield(Options,'Kprb')              Options.Kprb = 50; end;
if ~isfield(Options,'H0')                  Options.H0 = []; end;
if ~isfield(Options,'F0')                  Options.F0 = []; end;
if ~isfield(Options,'CleanSpec')    Options.CleanSpec = 0; end
if ~isfield(Options,'KCINITIAL')    Options.KCINITIAL = 15; end

% Now all options fields are defined. Check them for consistency.
if Options.Overlap>=Options.Nfft
   warning('(DISSRATE.M) Overlap greater or equal to NFFT. Setting Overlap to 0.');
   Options.Overlap = 0;
end
   

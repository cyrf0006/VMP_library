function [E,S,Options]  = epsprofile(dudz,z,W,nu,Options,accel)
% EPSPROFILE Compute a profile of dissipation rates.
% 
% Computes a depth profile of the dissipation rate by calling the
% function DISSRATE for consecutive sections of the shear signal
% dudz. The depth information is taken from the input vector z. W and nu
% are the sinking velocity and the viscosity, respectively. If the shear
% signal is to be decontaminated, a vector or matrix of accelerometer
% signals (one per column) may be included as a sixth input.  The results
% are returned in the two structures E and S (see 'Output parameters').  If
% a horizontal profiler was used, du/dx, distance, and instrument velocity
% may be passed instead of the first three input parameters.
% 
% Usage:
%   [E,S] = epsprofile(dudz,z,W,nu)
%   [E,S] = epsprofile(dudz,z,W,nu,Options)
%   [E,S] = epsprofile(dudz,z,W,nu,Options,accel)
% The length of each section is controlled by the Options structure. The
% fields .AveLen and .AveOverlap give the number of points and the
% amount of overlap of consecutive sections. To control the actual compu-
% tation of the dissipation rate, any valid field from the DISSRATE 
% Options structure can be specified (e.g. Options.Nfft).
%
% If Options.AveLen is a vector the dissipation rate is computed for the 
% depth bins centered on the elements in .AveLen. The width of the bins
% is controlled by Options.AveOverlap.  The bin-widths can be specified
% for each depth individually (if Options.AveOverlap is a vector) or 
% uniformly (if .AveOverlap is a scalar). For example, 
% Options.AveLen = 10:5:100 and Options.AveOverlap = 5 computes the 
% dissipation for 5-m bins at depth 10,15,20, etc.
% 
% [E,S,Options_out] = epsprofile(dudz,z,W,nu,Options)
% Returns the full set of options used in the computation.
% 
% Input parameters:
%    dudz                shear      [1/s]
%    z                   depth      [m or dbar]
%    W                   speed      [m/s]
%    nu                  viscosity  [m^2/s]
%    Options.AveLen      length of each section [points, default 1024]
%    Options.AveOverlap  overlap of consecutive sections [points, default 0]
%    Options.CleanSpec   whether or not to decontaminate shear spectra [1/0, default is 1]
%    accel               accelerometer signals, one per column [m/s^2] 
% 
% Output parameters:
%   E.e      Dissipation rate                         [W/kg]
%   E.z      Mean depth of the section                [m or dbar]
%   E.ks     Kolmogorov wavenumber                    [cpm]
%   E.kc     Cutoff wavenumber of the integration     [cpm]
%   E.nu     Viscosity used in the calculation        [m^2/s]
%   E.W      Sinking velocity used in the calculation [m/s]
%   
%   S.P      Disspation spectrum                      [s^(-2)/cpm]
%   S.k      Wavenumber vector for P                  [cpm]
%   S.Pn     Nasmyth spectrum for E.e                 [s^(-2)/cpm]
%   
%   Note: The fields in the E structure are vectors, which can be used in the
%   following way: e.g., plot(log10(E.e),E.z). 
% 
%   S is an Nx1 structure, where each field contains a vector
%   quantity. Element n in S contains the spectra for the n-th element
%   in E.e. For example, S(7).Pn is the Nasmyth spectrum for E.e(7).
% 
% Syntax summary:
%  [E,S] = epsprofile(dudz,z,W,nu)
%  [E,S] = epsprofile(dudz,z,W,nu,Options)
%  [E,S] = epsprofile(dudz,z,W,nu,Options,accel)
%  [E,S,Options_out] = epsprofile(...)
    
% (C) 2002 Rockland Oceanographic Services Inc.
% Author: Fabian Wolk
% Revision: 2002/08/01
% Isabelle Gaboury, 17 Mar. 2005 - added options to clean shear spectra

% Error checking:
%
error(nargchk(4,6,nargin));
warning('on') % make sure to display warning messages

% Set default options if necessary 
%
DEFAULTAVELEN = 1024;
%DEFAULTAVELEN = 512;

if nargin<5
    disp('Using default options to calculate dissipation rates');
    %DEFAULTAVELEN = 1024;
    Options.AveLen = DEFAULTAVELEN;
    Options.AveOverlap = 0;
    Options.Nfft = DEFAULTAVELEN;
    Options.CleanSpec = 0;
else
    if ~isfield(Options,'AveLen') Options.AveLen = DEFAULTAVELEN; end
    if ~isfield(Options,'AveOverlap') Options.AveOverlap = 0; end
    if ~isfield(Options,'Nfft') Options.Nfft = DEFAULTAVELEN; end
    if ~isfield(Options,'CleanSpec'), Options.CleanSpec = 0; end
end
if nargin<6
    accel = [];
    if ~isfield(Options,'CleanSpec') & Options.CleanSpec==1
        disp('No accelerometer data provided, unable to clean shear spectra');
        Options.CleanSpec=0;
    end
end 
if ~all([size(dudz)==size(z)] == [size(W)==size(nu)])
   error('Dimensions of all input vectors must agree.');
end

if length(Options.AveLen) == 1
   %
   % Branch for equidistant sections of length .AveLen 
   %
   
   % Compute the number of sections:
   N = length(dudz);
   nSec = fix((N-Options.AveOverlap)/(Options.AveLen-Options.AveOverlap));
   
   % Intialize variables:
   zz = zeros(nSec,1);
   E = struct('e',zz,'ks',zz,'kc',zz,'nu',zz,'W',zz,'z',zz);
   S = repmat(struct('k',zeros(Options.Nfft/2+1,1),...
      'P',zeros(Options.Nfft/2+1,1),...
      'Pn',zeros(Options.Nfft/2+1,1)), nSec,1);
   
   % Main loop:
   index = (1:Options.AveLen)';
   for n = 1:nSec
      fprintf('nSec = %3d: index = [%6d,%6d] \n',n,min(index),max(index));
      inan=find(isnan(dudz(index))==1); [jjj,kkk]=size(inan); clear inan
      if jjj ==0 
        [ee,ss,Options] = dissrate_db(dudz(index),W(index),nu(index),Options,accel(index,:));
        E.e(n) = ee.e;
        E.ks(n) = ee.ks;
        E.kc(n) = ee.kc;
        E.nu(n) = ee.nu;
        E.W(n) = ee.W;
        E.z(n) = mean(z(index));
        S(n) = ss;
      else
        E.e(n) = NaN;
        E.ks(n) = NaN;
        E.kc(n) = NaN;
        E.nu(n) = NaN;
        E.W(n) = ee.W;
        E.z(n) = mean(z(index));
        S(n) = NaN;
      end
      index = index+Options.AveLen-Options.AveOverlap;
   end
else
   %
   % Branch for computation at given depth bins:
   %
   
   N        = length(dudz);
   Bins     = Options.AveLen(:);
   BinWidth = Options.AveOverlap(:);
   nSec     = length(Bins); 
   if length(BinWidth) == 1
      BinWidth = BinWidth*ones(nSec,1);
   end
   if nSec > N
      error('There are more depth bins than in the shear signal');
   end
   if any(~BinWidth)
      warning('One or more bins have width zero. Check Options.AveOverlap parameter.')
   end
      
   % Intialize variables:
   zz = zeros(nSec,1);
   E = struct('e',zz,'ks',zz,'kc',zz,'nu',zz,'W',zz,'z',zz);
   S = repmat(struct('k',zeros(Options.Nfft/2+1,1),...
      'P',zeros(Options.Nfft/2+1,1),...
      'Pn',zeros(Options.Nfft/2+1,1)), nSec,1);
   E.z = Bins;
   
   % Make the depth vectors independend of the direction of
   % the coordinate system.
   depth = abs(z); 
   Bins = abs(Bins);
   BinWidth = abs(BinWidth)/2;
   
   % Main loop:
   for n = 1:nSec
      index = find( (Bins(n)-BinWidth(n))<=depth & depth<=(Bins(n)+BinWidth(n)) ); 
        
      if isempty(index)
         % The bin is not complete
         E.e(n) = nan;
         E.ks(n) = nan;
         E.kc(n) = nan;
         E.nu(n) = nan;
         E.W(n) = nan;
      else
         if length(index)<Options.Nfft
            warning('Bin width is smaller than FFT length.')
         end
         %fprintf('nSec = %3d: index = [%6d,%6d] ',n,min(index),max(index));
         %fprintf('depth = [%6.2f %6.2f] \n',z(index([1 end])));
         inan=find(isnan(dudz(index))==1); [jjj,kkk]=size(inan); clear inan
         if jjj == 0 
           [ee,ss,Options] = dissrate_db(dudz(index),W(index),nu(index),Options,accel(index,:));
           E.e(n) = ee.e;
           E.ks(n) = ee.ks;
           E.kc(n) = ee.kc;
           E.nu(n) = ee.nu;
           E.W(n) = ee.W;
           E.z(n) = mean(z(index));
           S(n) = ss;
         else
           E.e(n) = NaN;
           E.ks(n) = NaN;
           E.kc(n) = NaN;
           E.nu(n) = NaN;
           E.W(n) = NaN;
           E.z(n) = mean(z(index));
           %S(n) = NaN;
         end
      end   
   end % for loop
   
end % case of given depth bins

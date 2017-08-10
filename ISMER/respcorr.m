function Pc = respcorr(P,k,varargin)
% RESPCORR Spectral correction for shear probe response.
%
% Adjusts the spectrum P defined at wavenumbers (or frequencies)
% k by the probe response function H(k) = 1/[1+(k/kc)^p], where 
% kc = 50 is the 3-db point of the assumed filter transfer function of 
% the probe, and p = 2 is the filter order. For arrays, the function
% works along the columns of P and k. 
%
% Usage:
% Pc = respcorr(P,k)
% Pc = respcorr(P,k,kc)
% Pc = respcorr(P,k,kc,p) can be used to specify different values 
% for kc and p.
%
% Input parameters:
%   P   spectrum to be adjusted [m/s^2/cpm]
%   k   wavenumber bins of P [cpm]
%   kc  cutoff wavenumber [cpm, default 50]
%   p   filter order      [default 2]
% 
% Ouput parameters:
%   Pc  the adjusted spectrum


% Ref: Oakey, N., 1982, J. Phys. Ocean., 12, 256-271.

% (C) 2002 Rockland Oceanographic Services Inc.
% Author: Fabian Wolk
% Revision: 2002/08/01

error(nargchk(1,4,nargin));
% -FC modifs. for Octave
if nargin < 4 || ~exist(p) 
    p = 2; 
end
if (nargin < 3 || ~exist(kc))
    kc = 50; 
end

if ~all(size(P) == size(k)),
   error('Dimensions of P and k must agree.')
end

[r,c] = size(P);
Pc = zeros(r,c);
for n = 1:c
   Pc(:,n) = P(:,n).*(1+(k(:,n)/kc).^p);
end

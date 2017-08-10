function y=integral(P,k,kmin,kmax);
% INTEGRAL Integrate spectrum (or function) in specified wavenumber range.
%
% Integrates the spectrum P in the frequency (or wavenumber) space k
% between kmin and kmax. Uses trapezoidal integration. 
%
% Usage: y = integral(P,k,k1,k2);
%
% Input parameters:
%   P     Spectrum or signal to be integrated
%   k     Abscissa vector for which P is defined, i.e., P=P(k)
%   k1    Lower integration limit
%   k2    Upper integration limit
%
% Ouput parameters:
%   y     Value of the integral

% (C) 2002 Rockland Oceanographic Services Inc.
% Author: Fabian Wolk
% Revision: 2002/08/01

error(nargchk(2,4,nargin));
if nargin == 2
   y = trapz(k,P);
   return
end
if nargin < 4
   kmax = max(k);
end
if nargin < 3
   kmin = 0;
end
if ~(size(P) == size(k))
   error('Dimensions P and k must agree.');
end

Nk = length(k);
NP = length(P);

% the integration range
dk = kmax-kmin;
if dk<0,
   error('k2 must be strictly greater than k1')
end

% integration range must straddle at least one point k vector
if (kmin > k(Nk) | kmax < k(1)),
   warning('Integration limits are out of range');
   y = 0;
   return;
end

if Nk == 1; % there is only one point in the vector and it is within dk
   y = P(1)*dk;
else
   index = find(k>=kmin & k<=kmax);
   n = length(index);
   if n==0, % dk is between two points ==> interpolate
      x1= find(k<kmin); % find the point just below kmin
      x1 = x1(length(x1));
      x2 = find(k>kmax); % find the point just above kmax
      x2 = x2(1);
      y = 0.5*(P(x1)+P(x2)) * dk;
   elseif n==1 % dk straddles one point (which could be the first or the last one)
      index = find(k<=kmax); % the point just below kmax lies in the interval dk
      index = index(length(index));
      y = P(index)*(kmax-kmin);
   else   % dk includes a range of points
      y = trapz(k(index), P(index));
   end
end

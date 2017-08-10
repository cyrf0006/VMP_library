function r = fc_chi2rnd(v,varargin);
%FC_CHI2RND Random arrays from chi-square distribution.
%   R = FC_CHI2RND(V) returns an array of random numbers chosen from the
%   chi-square distribution with V degrees of freedom.  The size of R is
%   the size of V.
%
%   R = FC_CHI2RND(V,M,N,...) or R = FC_CHI2RND(V,[M,N,...]) returns an
%   M-by-N-by-... array. 
%
%   See also CHI2CDF, CHI2INV, FC_CHI2PDF, CHI2STAT, NCX2RND, RANDOM.

%   FC_CHI2RND generates values using the definition of the chi-square
%   distribution, as a special case of the gamma distribution.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.

% No copyright, FD modif of matworks fn
%   $Revision: 2.12.4.3 $  $Date: 2004/01/24 09:33:12 $

if nargin < 1
    error('stats:fc_chi2rnd:TooFewInputs','Requires at least one input argument.');
end

[err, sizeOut] = fc_statsizechk(1,v,varargin{:});
if err > 0
    error('stats:fc_chi2rnd:InputSizeMismatch','Size information is inconsistent.');
end

% Generate gamma random values, and scale them.
r = 2.*randg(v./2, sizeOut);

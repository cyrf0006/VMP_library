function y = fc_chi2pdf(x,v)
%FC_CHI2PDF Chi-square probability density function (pdf).
%   Y = FC_CHI2PDF(X,V) returns the chi-square pdf with V degrees  
%   of freedom at the values in X. The chi-square pdf with V 
%   degrees of freedom, is the gamma pdf with parameters V/2 and 2.
%
%   The size of Y is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also CHI2CDF, CHI2INV, FC_CHI2RND, CHI2STAT, PDF.

%   Notice that we do not check if the degree of freedom parameter is integer
%   or not. In most cases, it should be an integer. Numerically, non-integer
%   values still gives a numerical answer, thus, we keep them.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986, pages 402-403.
%
% No copyright, FD modif of matworks fn
%   $Revision: 2.12.2.2 $  $Date: 2004/07/05 17:02:23 $

if nargin < 2, 
    error('stats:fc_chi2pdf:TooFewInputs','Requires two input arguments.'); 
end

[errorcode x v] = distchck(2,x,v);

if errorcode > 0
    error('stats:fc_chi2pdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y=zeros(size(x));

y = fc_gampdf(x,v/2,2);

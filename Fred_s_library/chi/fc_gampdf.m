function y = fc_gampdf(x,a,b)
%FC_GAMPDF	Gamma probability density function.
%	Y = FC_GAMPDF(X,A,B) returns the gamma probability density function 
%	with parameters A and B, at the values in X.
%
%	The size of Y is the common size of the input arguments. A scalar input  
%	functions as a constant matrix of the same size as the other inputs.	 
%
%	Some references refer to the gamma distribution with a single
%	parameter. This corresponds to the default of B = 1.

%	References:
%	   [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%	   Springer-Verlag, 1986, pages 401-402.

%       No copyright, FD modif of matworks fn
%	$Revision: 1.2 $  $Date: 1993/09/10 22:00:05 $

if nargin < 3, 
    b = 1; 
end

if nargin < 2, 
    error('Requires at least two input arguments'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

% Initialize Y to zero.
y = zeros(size(x));

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    y(k1) = NaN * ones(size(k1));
end

k=find(x > 0 & ~(a <= 0 | b <= 0));
if any(k)
    y(k) = (a(k) - 1) .* log(x(k)) - (x(k) ./ b(k)) - gammaln(a(k)) - a(k) .* log(b(k));
    y(k) = exp(y(k));
end
k1 = find(x == 0 & a < 1);
if any(k1)
  y(k1) = Inf*ones(size(k1));
end

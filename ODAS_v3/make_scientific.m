%% make_scientific
% Convert number into a string in scientific notation
%%
% <latex>\index{Type B!make\_scientific}</latex>
%
%%% Syntax
%   scString = make_scientific( N, Digits )
%
% * [N] Number to convert to scientific notation
% * [Digits] Number of significant digits required in the output
% * []
% * [scString] String representation of 'N' in scientific notation.  When 
%         'N' is a vector of length greater than 1, 'scString' is a cell 
%         array of strings.
%
%%% Description
% Converts the number 'N' into a string representation in scientific 
% notation. The resulting string is optimized for use in Matlab plotting.
% When a matrix of values is inputed, each value is converted and a cell
% array of strings is returned.
%
%%% Examples
%
%   >> plotLabel = make_scientific( pi, 5 )
%
% Returns Pi with 5 significant digits: ``3.1416 \times 10^{0}''

% Version History
%
% * 2012-09-03 (WID) documentation changed for Matlab publishing
% * 2012-09-03 (WID) added support for matrix input / string array output
% * 2012-11-05 (WID) documentation update

function scString = make_scientific(N, Digits)

if (~isnumeric(N) || ~isnumeric(Digits))
    error ('Input must be a number')
end

scString = {};

for Num = N,
    n = find(Num < 0); % find negative values
    Num(n) = -Num(n); % temporarily make negative values positive
    M = log10(Num); % take base 10 logarithm
    exponent = floor(M);
    mantissa = M - floor(M);
    mantissa = 10^(mantissa);
    mantissa(n) = -mantissa(n); % return negative values to les than zero 
    sc_string = [num2str(mantissa, Digits) ' \times 10^{' num2str(exponent) '}'];
    scString{end + 1} = sc_string;
end

if length(scString) == 1, scString = scString{1}; end

end
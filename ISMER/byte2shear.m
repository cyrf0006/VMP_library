function shear = byte2shear(sh,S,G,w);

% sh:  2-byte number.
% S:   The shear probe sensitivity.
% G:   Gain of the shear circuit relative to an ideal differentiator.
% shear: The shear in s^{-1}.
%
% See Application note AN-005 from Rockland.
%

if (length(w)>1 & length(sh) ~= length(w))
  error('Input vectors must be the same length');
end

shear = sh*(5/2^16)./(2*sqrt(2)*G*S*w.^2);

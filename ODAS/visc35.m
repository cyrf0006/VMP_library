%% visc35
% Approximation for the kinematic viscosity of seawater for S = 35.0
%%
% <latex>\index{Type A!visc35}</latex>
%
%%% Syntax
%
%   v = visc35( t )
%
% * [t] temperature in degrees Celsius 
% * []
% * [v] viscosity in m^2/s 
%
%%% Description
%
% Return an approximation of the kinematic viscosity, based on temperature (in
% degrees C). The viscosity is derived from a 3-rd order polynomial fit of nu
% against T for salinity 35. The error of the approximation is less than 1% for
% (30 <= S <= 40) and (0 <= T <= 35) at atmospheric pressure.
%
% (see also viscosity)

% *Version History*
%
% * 1999-06-01 (FW)
% * 2002-10-09 (FW) revised
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-24 (WID) documentation update for publishing

function v = visc35(t)

pol=[-1.131311019739306e-011
      1.199552027472192e-009
     -5.864346822839289e-008
      1.828297985908266e-006];
v = polyval(pol,t);    
     

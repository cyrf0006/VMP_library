%% visc00
% Approximation for the kinematic viscosity of freshwater for _S_ = 0.0
%%
% <latex>\index{Type A!visc00}</latex>
%
%%% Syntax
%
%   v = visc00( t )
%
% * [t] temperature in degrees Celsius 
% * []
% * [v] viscosity in m^2/s 
%
%%% Description
%
% Returns an approximation of the kinematic viscosity, based on temperature (in
% degrees C). The viscosity is derived from a 3-rd order polynomial fit of nu
% against T for salinity 0. The error of the approximation is less than 1% for
% (0 <= T <= 20) at atmospheric pressure.
%
% (see also viscosity)

% Version History:
%
% * 1999-06-01 (FW)
% * 2002-10-09 (FW) revised
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-24 (WID) documentation update for publishing

function v = visc00(t)

pol=[-1.8377764060e-011
      1.4447664328e-009
     -6.0996917629e-008
      1.7903678481e-006];
v = polyval(pol,t);    
     

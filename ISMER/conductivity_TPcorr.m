function sbc_corr=conductivity_TPcorr(sbc,sbt,pres,CPcor,CTcor);

% CONDUCTIVITY_TPCORR - correct conductivity data for thermal expansion/compressibility
% sbc_corr = conductivity_TPcorr(sbc,sbt,pres,CPcor,CTcor);
%
% Correct conductivity data for the effects of thermal expansion and
% compressibility.  
% Inputs: conductivity, temperature, and pressure (all must be in physical
%   units: mS/cm, deg. C, dBar); coefficients for compressibility and
%   thermal expansion (optional, default values included in this M-file)
% Output: corrected conductivity
%
% Isabelle Gaboury, 14 Feb. 2005. From OTL sb_c.m routine.

% Default coefficients for compressibility and thermal expanion (from
% Sea-Bird calibration sheet, May 2004)
if nargin<5
    [CTcor,CPcor] = parameter_list('odas','CTcor','CPcor');
%     CPcor = -9.5700e-8;
%     CTcor = 3.2500e-6;
end

% Correct the conductivity data.  Note that this is typically a fairly
% small correction.
sbc_corr = sbc./(1+CTcor*sbt+CPcor*pres);
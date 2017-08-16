%% cal_FP07_in_situ
% Calibrate thermistor probe using in situ data
%%
% <latex>\index{Functions!cal\_FP07\_in\_situ}</latex>
%
%%% Syntax
%   [T_0, beta, Lag] = cal_FP07_in_situ( file_name, T_ref, T, SN, ... )
%
% * [file_name] Name of mat-file containing data used to
%       calibrate a thermistor.
% * [T_ref] Name of vector within the mat-file that contains the
%       reference temperature in degrees celsius. Usually from a SBE4F
%       thermometer or a JAC_CT.
% * [T] Name of thermistor to calibrate, typically 'T1' or 'T2'.
% * [SN] Serial number of thermistor.
% * [...] Optional elements describing a profile.  Can be provided in a
%       structure or as key / value pairs.  See below for more details.
% * []
% * [T_0] Value of parameter T0, used in the Steinhart-Hart equation. When
%       called with no input parameters, a structure containing default 
%       input parameter values is returned for reference.
% * [beta] beta coefficients, in ascending order, of the fit to the SH 
%       equation. i.e. beta_1, beta_2, beta_3
% * [Lag] Delay in seconds between the thermistor and the reference
%       thermometer. Typically a negative value because the reference
%       sensor is usually behind the thermistor being calibrated.
%
%%% Description
%
% Function to calibrate a FP-07 thermistor probe using in-situ data.  Reliable
% temperature data must be within the specified mat-file. They will be the 
% refercence temperature -- usually a Sea-Bird SBE3F or a JAC-CT. 
%
% This function makes 5 figures. Figure 1 shows you the portion of the data
% file that will be used for the calibration. It then detrends the
% thermistor data and scales it so that it is approximately aligned with
% the reference thermometer (figure 2). It then plots the cross-correlation
% coefficient between the thermistor and the reference thermometer and
% estimates the lag between these two signals (Figure 3). Next it plots the
% the natural logarithm of the resistance ratio against the inverse of the
% absolute temperature to show the quality of the regression (Figure 4).
% Finally, it plots a depth profile of the lagged reference temperature,
% the newly calibrated thermistor temperature, and their difference.
%
% Optional input parameters include;
%
% * [profile_num] Integer identifying the profile number to be
%       analyzed. Default = 1. Some files may contain multiple profiles.
% * [direction] String identifying the direction of profiling.  Either
%      'up' or 'down' - default = 'down'.
% * [min_duration] Minimum time, in seconds, for a profile to be deemed a
%       profile.  Direction of travel must be monotonic, rate of change of 
%       pressure must be above a minimum level, and pressure must 
%       exceed a minimum value. Default = 20 s.
% * [P_min] Minimum pressure for a profile. Default = 1 dBar
% * [min_speed] Minimum speed of a profile, actually dP/dt.
%       Default = 0.4 dBar/s. Use a smaller value, ~0.1 dBar/s, for
%       gliders.
% * [order] Fit order to the Steinhart-Hart equation. Value can be 1,
%       2, or 3. Default = 2.
%
%___________________________
%
% Version History
%
% * 2013-12-05 (RGL) original version.
% * 2013-12-06 (RGL) added varargin to define a profile with default values.
% * 2015-04-10 (WID) revised documentation for publishing.
% * 2015-04-27 (RGL) modified to allow the specification of the fit order.
% * 2015-07-27 (WID) use texstr in place of fix_underscore.
% * 2015-07-29 (WID) return default values when called with no input
%                       parameters.
% * 2015-10-27 (RGL) Changed description section.

function [T_0,beta,Lag] = cal_FP07_in_situ(file_name,T_ref,T,SN,varargin)

%%%%
% Default values for optional fields
default_direction     = 'down'; % downwards profiling is the default
default_P_min         = 1; % in dBar
default_min_speed     = 0.4; % in dbar/s
default_profile_num   = 1; % process the first profile
default_min_duration  = 20; % minimum duration [s] for a segment to be considered a profile
default_order         = 2; % The order of the fit to the SS equation
if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    T_0 = result;
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_string      = @(x) ischar(x);

addRequired(  p, 'file_name',    val_string);
addRequired(  p, 'T_ref',        val_string);
addRequired(  p, 'T',            val_string);
addRequired(  p, 'SN',           val_string);
addParamValue(p, 'profile_num',  default_profile_num,    val_numeric);
addParamValue(p, 'P_min',        default_P_min,          val_numeric);
addParamValue(p, 'min_speed',    default_min_speed,      val_numeric);
addParamValue(p, 'direction',    default_direction,      val_string);
addParamValue(p, 'min_duration', default_min_duration,   val_numeric);
addParamValue(p, 'order',        default_order,          val_numeric);

% Parse the arguments.
parse(p, file_name, T_ref, T, SN, varargin{:});

% Perform last stages of input validation.

P_min         = p.Results.P_min;
min_speed     = p.Results.min_speed;
profile_num   = p.Results.profile_num;
direction     = p.Results.direction;
min_duration  = p.Results.min_duration;
order         = p.Results.order;

% end of input argument checking.
%%%%
% Get started by making sense of the inputs

T_string = [T '_d' T]; % should be T1_dT1 or T2_dT2
T_without_pre_emphasis_string = T; % need this to get coefficients from the configuration file string
T_ref_string = T_ref;

load(file_name, 'setupfilestr', T, T_string, T_ref,'t_fast', 't_slow', ...
    'P_slow', 'P_fast', 'W_slow', 'fs_fast', 'fs_slow')

ratio = round(fs_fast / fs_slow); % sampling ratio
eval(['T = ' T_string ';']) % T is the thermistor signal with pre-emphasis
eval (['T_ref = ' T_ref ';']) % T_ref is the reference thermometer, usually SBT
eval(['T_without_pre_emphasis = ' T_without_pre_emphasis_string ';']) % usually T1 or T2

T = deconvolve(T_string,[], T, fs_fast, setupfilestr,6);

T = reshape(T, ratio, []); % down size to match T_ref
T = mean(T)';


title_string{1} = ['Thermistor in situ calibration , SN-' SN];

profile = get_profile(P_slow, W_slow, P_min, min_speed, direction, min_duration, fs_slow);

profile_start = profile(1,profile_num);
profile_end   = profile(2,profile_num);
m = (profile_start:profile_end)';

%%%%
% Provide some plots to varify the operations
% first test
figure(1)
plot(t_slow, T, t_slow(m), T(m),'r');grid on 
legend(...
     T_without_pre_emphasis_string, ...
    [T_without_pre_emphasis_string ' Profile'], 1)
title (title_string)
xlabel('\it t \rm [s]')
ylabel('[counts]')

% Visualize the lag and correlation between the thermistor and the
%   reference thermometer.
figure(2)
junk_T      = detrend(T(m));
junk_T_ref  = detrend(T_ref(m));
range_T_ref = max(junk_T_ref) - min(junk_T_ref);
range_T     = max(junk_T) - min(junk_T);
junk_T      = junk_T * range_T_ref / range_T; % T should now span the same range as T_ref.

title_string{2} = [...
    'Detrended ' texstr(T_ref_string) ' & scaled ' texstr(T_without_pre_emphasis_string)];

plot(t_slow(m), [junk_T  junk_T_ref]);grid on
title (title_string)
xlabel('\it t \rm [s]')
ylabel('[ ^{\circ}C ]')
legend(texstr(T_without_pre_emphasis_string), T_ref_string)

figure(3) % plot the cross-correlation coefficient for T and T_ref
title_string{2} = [...
    'Cross-correlation coefficinet of ' T_without_pre_emphasis_string ...
    ' and ' T_ref_string];
max_lag = round(5*fs_slow); % my estimate of the max lag required to find the actual lag.
[bb, aa] = butter(2,4/(fs_slow/2)); % 4 Hz smoother to suppress high-frequency noise

[correlation, lags] = xcorr(...
    filter(bb,aa,detrend(diff(T(m)))),...
    filter(bb,aa,detrend(diff(T_ref(m)))),max_lag,'coeff');
[max_corr, m_lag] = max(abs(correlation));
junk_m = m_lag; % needed for figure
m_lag = m_lag - max_lag - 1;

plot(lags/fs_slow, correlation, m_lag/fs_slow, correlation(junk_m), 'r*');grid on
xlabel('Lag [ s ]')
legend_string_1 = ['max X_{corr} = ' num2str(max_corr,2)];
legend_string_2 = ['@ \tau = ' num2str(m_lag/fs_slow,2) ' s'];
legend(legend_string_1, legend_string_2)
title(title_string)

%%%%
% Now the hard work begins
% First align the T and T_ref signals using m_lag. m_lag is expected to be
%   negative.
T_ref2 = T_ref(m); % Only the profile and use a copy
T2     = T(m);

if m_lag >0, m_lag = 0; end

T_ref2 = T_ref2(1-m_lag:end);
T2     = T2(1:end+m_lag);
T_ref_regress = T_ref2 + 273.15; % in kelvin
T_ref_regress = 1 ./ T_ref_regress;

% Now gather information about the electronics for this thermistor.
my_object = setupstr( setupfilestr );
E_B = str2double(char(setupstr( my_object, T_without_pre_emphasis_string, 'E_B')));
a   = str2double(char(setupstr( my_object, T_without_pre_emphasis_string, 'a'  )));
b   = str2double(char(setupstr( my_object, T_without_pre_emphasis_string, 'b'  )));
G   = str2double(char(setupstr( my_object, T_without_pre_emphasis_string, 'G'  )));
adc_fs   = str2double(char(setupstr( my_object, T_without_pre_emphasis_string, 'adc_fs'  )));
adc_bits = str2double(char(setupstr( my_object, T_without_pre_emphasis_string, 'adc_bits'  )));
factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);

Z = factor*(T2 - a)/b;
RT_R0 = (1 - Z) ./ (1 + Z); % This is the resistance ratio for this thermistor.
RT_R0 = log(RT_R0);
% Next confirm that the thermistor follows the Steinhart-Stein equation. The
%   plot should be a nearly straight line.

figure(4)
plot(T_ref_regress, RT_R0, '.');grid on
xlabel ('T^{-1} [k^{-1}]')
ylabel ('log_e (R_T / R_0)')

beta = zeros(1,order);
% Generate the coefficients for this thermistor.
p = polyfit(RT_R0, T_ref_regress, order);
p = 1 ./ p;
p = fliplr(p); % place in ascending order
T_0    = p(1);
for index = 2:order+1
    beta(index-1) = p(index);
end
Lag    = m_lag / fs_slow; % in seconds and should be negative.

title_string{2} = ...
    ['T_0 = ' num2str(T_0) ', \beta = ' num2str(beta) ];
title(title_string)
x_limits = get(gca,'xlim');
x_text = x_limits(1) + (x_limits(2) - x_limits(1))/20;
y_limits = get(gca,'ylim');
y_text = y_limits(2) - (y_limits(2) - y_limits(1))/20;
my_text = [...
    '$$\frac{1}{T} = \frac{1}{T_0} + ' ...
    '\frac{1}{\beta_1} \log_{\ e}\left(\frac{R_T}{R_0}\right) + ' ...
    '\frac{1}{\beta_2} \log^{\ 2}_{\ e} \left(\frac{R_T}{R_0}\right)$$'];

text(x_text, y_text, my_text, 'interpreter','latex')

%%%%
% Now make a profile to show how well the calibration worked
figure(5)
Z = factor*(T - a)/b;
RT_R0 = (1 - Z) ./ (1 + Z); % This is the resistance ratio for this thermistor.
RT_R0 = log(RT_R0);
pp = fliplr(1./p);

T_calibrated = polyval(pp, RT_R0);
%T_calibrated = 1 ./ T_0 + RT_R0/beta_1 + RT_R0.^2/beta_2;
T_calibrated = 1 ./ T_calibrated;
T_calibrated = T_calibrated - 273.15;

plot(...
    T_ref(m-m_lag), P_slow(m), ...
    T_calibrated(m), P_slow(m), ...
    T_calibrated(m) - T_ref(m-m_lag), P_slow(m));grid on
set(gca,'ydir','rev')
xlabel('[ ^{\circ}C]')
ylabel('\itP \rm[dBar]')
title(title_string{1})
legend(...
    T_ref_string, ...
    T_without_pre_emphasis_string, ...
   [T_without_pre_emphasis_string  ' - ' T_ref_string], 4)

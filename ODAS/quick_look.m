%% quick_look
% Visualize contents of a RSI data file, compute spectra for a selected
% range, and return a profile of the rate of dissipation of kinetic energy.
%%
% <latex>\index{Functions!quick\_look}</latex>
%
%%% Syntax
%   diss = quick_look (fname, P_start, P_end, ql_info, ...)
%
% * [fname] Name of the binary data file to process (extension optional).
% * [P_start] Starting point, in pressure, of the segment to be used for
%           displaying a spectrum of all shear, scalar-gradient, and
%           acceleration signals. Can be empty, if P_end is also empty, to
%           suppress the display of the spectra. 
% * [P_end] End point, in pressure, for the segment used to show spectra. Can
%           be empty, if P_start is also empty. 
% * [ql_info] Structure containing configuration parameters. A template 
%           is generated and returned when quick_look is called with no
%           input parameters. The parameters are described below.
% * [...] Optional configuration parameters to supplement, or override,
%           those values included within ql_info.
% * []
% * [diss] Dissipation estimates for the specified profile. A default
%           ql_info structure is returned when quick_look is called with no
%           input parameters.
%
%%% Description
% This function generates a variety of figures that visualize the data in
% the file $\texttt{fname}$. The function works in three stages. In the
% first stage, data is converted into physical units using the
% $\texttt{odas\_p2mat}$ function, if a mat-file does not already exist. In
% the second stage, the data are plotted in several figures. In the third
% stage, the function computes a profile of the rate of dissipation of
% kinetic energy from the shear probe data.
%
% Converting data into physical units is performed by the
% $\texttt{odas\_p2mat}$ function. $\texttt{quick\_look}$ calls this
% function if a mat-file with the same name as the data file does not
% exist. Or, if the parameters for the $\texttt{odas\_p2mat}$ function
% (which may be included in the call to the $\texttt{quick\_look}$
% function) are not the same as those that were used to create the
% mat-file. You will be prompted for permission to overwrite the mat-file.
% Therefore, input parameters used with $\texttt{quick\_look}$ can include
% the parameters required to convert your data into physical units, but
% these parameters are only required if the mat-file does not exist. 
%
% If $\texttt{quick\_look}$ is called without input arguments, it returns
% the default parameters used by $\texttt{quick\_look}$ and by
% $\texttt{odas\_p2mat}$, within a single structure, that you can customize
% to your particular processing requirements. For example,
%  
%    >> ql_info = quick_look
% 
%    ql_info = 
%                   HP_cut: 0.4000
%                   LP_cut: 30
%                     YD_0: 0
%                despike_A: [8 0.5000 0.0400]
%                despike_C: [10 1 0.0400]
%               despike_sh: [8 0.5000 0.0400]
%              diss_length: 8
%                     f_AA: 98
%                  f_limit: Inf
%               fft_length: 2
%                fit_2_isr: 1.5000e-05
%                fit_order: 3
%             make_figures: 1
%                  op_area: 'open_ocean'
%                  overlap: 4
%            profile_min_P: 1
%            profile_min_W: 0.2000
%     profile_min_duration: 20
%              profile_num: 1
%          MF_extra_points: 0
%                     MF_k: 4
%                 MF_k_mag: 1.7000
%                   MF_len: 256
%                MF_st_dev: 'st_dev'
%             MF_threshold: []
%                      aoa: []
%           constant_speed: []
%            constant_temp: []
%               hotel_file: []
%             speed_cutout: 0.0500
%                speed_tau: []
%              time_offset: 0
%                  vehicle: ''
% 
% The configuration parameters (fields) within the structure
% $\texttt{ql\_info}$ that control the behaviour of $\texttt{quick\_look}$,
% are listed below. They are grouped for clarity, but are all part of the
% single structure $\texttt{ql\_info}$, including those required only by
% $\texttt{odas\_p2mat}$. 
%
%%% Parameters that specify a profile 
% A single data file may contain multiple profiles. For example, a verticle
% profiler that was raised and lowered multiple times while recording data
% into a single file. $\texttt{quick\_look}$ detects profiles using 5
% specifications. These are:
%
% * [profile_num] Index to the requested profile from the set of detected 
%       profiles.  The first profile, 1, is the default value.
% * [profile_min_P] The minimum pressure of a profile, [dbar]. Default = 1.
% * [profile_min_W] The minimum vertical speed of profiling, [dbar/s].
%      Default = 0.2 
% * [profile_min_duration] The minimum duration in which the minimum
%      pressure and speed must be satisfied [s]. Default = 20.
%
% The $\texttt{direction}$ of profiling is the implict 5th specification
% and is determined by the $\texttt{vehicle}$ used for profiling. See the
% file $\texttt{default\_vehicle\_attributes.ini}$ for the direction of each
% recognized $\texttt{vehicle}$.
%
%%% Parameters that control the calculation of dissipation rate 
% The rate of dissipation of turbulent kinetic energy, $\epsilon$, is
% estimated using the function $\texttt{get\_diss\_odas}$ and follows the
% method described in RSI Technical Note 028. The controlling parameters
% are:
%
% * [diss_length] The time span, in seconds, over which to make each
%      estimate of the rate of dissipation. Default = 8. 
% * [overlap] The overlap, in seconds of each dissipation estimate. Default
%      = 4.
% * [fft_length] The length, in seconds, of the fft-segments that will be
%      ensemble-averaged into a spectrum. Default = 2.
% * [fit_2_isr] The rate of dissipation, in W/kg, above which the function will
%      switch from the method of spectral-integration to the method of
%      fitting to the inertial subrange, to estimate the rate of
%      dissipation. Default = 1.5e-5.
% * [fit_order] The order of the polynomial to fit to the shear spectrum,
%      in log-log space, to estimate the wavenumber at which the spectrum
%      has a minimum. This is one of several constraints on the upper limit
%      of spectral integration. Default = 3.
% * [ f_limit] The upper frequency limit, in Hz, to be used for estimating
%      the rate of dissipation. Default = inf (no unconditional limit).
% * [f_AA] The cut-off frequency, in Hz, of the anti-aliasing filter in
%      your instrument. Default = 98. This value is instrument dependent
%      but is almost always 98, unless you have an instrument that has been
%      customized to samples at rates faster than 512 per second.
%
%%% Parameters that control the processing of microstructure signals 
% The parameters that control the processing of the microstructure signals
% pertain to the high-pass filtering of the shear-probe signals, and the
% despiking of the shear-probe, acceleration and micro-conductivity signals.
% 
% * [HP_cut] The cut-off frequency, in Hz, of the high-pass filter applied
%      to the shear-probe signals. Default = 0.4. 
% * [despike_sh] The triplet of parameters for the despike
%      function applied to the shear-probe signals. The first value is the
%      threhold. The second value is the cut-off frequency, in Hz, of the
%      low-pass smoothng filter. The third value is the duration, in
%      seconds, of data to remove around a spike. Default = [8 0.5 0.04].
%      See the function despike for more information on the parameters. You
%      can suppress the despike function by specifying an infinite
%      threshold, for example [inf 0.5 0.07].
% * [despike_A] The triplet of parameters for the despike
%      function applied to the accelerometer signals. For data collected
%      with a glider, it may be necessary to supress despiking. The
%      intermittent vibrations from battery movement and fin actuators
%      creates short duration vibrations that are easily confused with
%      spikes, but such data is needed for coherent-noise removal. Default
%      = [8 0.5 0.04].
% * [despike_C] The triplet of parameters for the despike function applied
%      to the micro-conductivity signals. Default = [10 1.0 0.04]. 
%
%%% Parameters for other purposes 
% 
% * [LP_cut] The cut-off frequency, in Hz, of the low-pass filter applied
%      to the microstructure profile signals, for graphical display only.
%      It does not affect the estimation of the rate of dissipation.
%      Default = 30.
% * [YD_0] The year-day subtracted from the time axis of figures. It is
%      currently not used. Default = 0.
% * [op_area] The operational area of your instrument. Recognized values are
%     'open_ocean' and 'tidal_ch'. It controls the
%     scale on certain figures. Default = 'open\_ocean'.
% * [make_figures] The parameter that determins if figures are generated.
%     make_figures = false suppresses the generation of figures to speed up
%     the data processing. The default is make_figures = true.
%
%%% Parameters for the odas_p2mat function 
% The parameters starting with $\texttt{MF\_extra\_points}$ are used by the
% $\texttt{odas\_p2mat}$ function, to convert your data into physical
% units. They are described in the section for that function. 
%
%
%%% The diss structure output
% The output structure from $\texttt{quick\_look}$ depends slightly on the
% channels in your instrument. The typical fields are shown below, in groups.
% The first group is associated with the calculation of the profile of the
% rate of dissipation, $\epsilon$, and is described in the section for
% $\texttt{get\_diss\_odas}$.
%
%                         [e]  [2x217 double]
%                     [K_max]  [2x217 double]
%                   [warning]  [2x217 double]
%                    [method]  [2x217 double]
%              [Nasmyth_spec]  [513x2x217 double]
%                        [sh]  [4-D double]
%                  [sh_clean]  [4-D double]
%                        [AA]  [4-D double]
%                        [UA]  [4-D double]
%                         [F]  [513x217 double]
%                         [K]  [513x217 double]
%                 [Data_fast]  [18x217 double]
%                 [Data_slow]  [19x217 double]
%                     [speed]  [217x1 double]
%                        [nu]  [217x1 double]
%                         [P]  [217x1 double]
%                         [T]  [217x1 double]
%                         [t]  [217x1 double]
%                       [AOA]  []
%                      [f_AA]  88.2000
%                   [f_limit]  Inf
%                 [fit_order]  3
%               [diss_length]  4096
%                   [overlap]  2048
%                [fft_length]  1024
%
% The next group is associated with the despiking of the shear-probe,
% micro-conductivity and acceleration signals.
%
%                  [spikes_A]  {[4x1 double], []}
%              [pass_count_A]  [2x1 double]
%                [fraction_A]  [2x1 double]
%                 [spikes_sh]  {[247x1 double], [269x1 double]}
%             [pass_count_sh]  [2x1 double]
%               [fraction_sh]  [2x1 double]
%                  [spikes_C]  {}
%              [pass_count_C]  [0x1 double]
%                [fraction_C]  [0x1 double]
%
% The indices to the spikes located in the signals from the accelerometers,
% are given in $\texttt{spikes\_A}$. The number of passes of the despike
% function used to remove the spikes is in $\texttt{pass\_count\_A}$. The
% fraction of data removed by the despike function is in
% $\texttt{fraction\_A}$. Similarly for the shear probe and the
% micro-conductivity signals.
%
% The next group is associated with the scalar signals. 
%
%            [scalar_spectra]  [1x1 struct]
%        [scalar_vector_list]  {'gradT1'  'gradT2'}
%               [scalar_info]  [1x1 struct]
% 
% The structure
% $\texttt{scalar\_spectra}$ (a structure within a structure) is described in
% the section for the function $\texttt{get\_scalar\_spectra\_odas}$. The
% processing parameters are in $\texttt{scalar\_info}$. The names of the
% scalar vectors are in $\texttt{scalar\_vector\_list}$.
%
% The remaining fields are: 
%
%                   [fs_fast]  511.9454
%                   [fs_slow]  63.9932
%     [profiles_in_this_file]  1
%              [speed_source]  'Rate of change of pressure'
%                 [fast_list]  {1x18 cell}
%                 [slow_list]  {1x19 cell}
%                   [ql_info]  [1x1 struct]
%                [ql_info_in]  []
%
% $\texttt{fs\_fast}$ and $\texttt{fs\_slow}$ are the actual fast and slow
% sampling rates of the data. The number of profiles detected in this data
% file is given in $\texttt{profiles\_in\_this\_file}$. The source of the
% speed of profiling is identified in $\texttt{speed\_source}$.
%
% All column vectors that have a length matching the length of the time
% vector $\texttt{t\_fast}$ are combined into a single matrix and passed to
% the dissipation function for averaging over the interval of each
% dissipation estimate. Similarly for all column vectors that match
% $\texttt{t\_slow}$. Each row of $\texttt{Data\_fast}$ and
% $\texttt{Data\_slow}$ hold the values from a single vector. The names of
% the signals are identified in the cell arrays $\texttt{fast\_list}$ and
% $\texttt{slow\_list}$. In this example there are 19 fast vectors and 19
% slow vectors. Use the Matlab $\texttt{find}$ and $\texttt{strcmpi}$
% functions to identify the row of a particular signal.
%
% The structure that was input to this function is saved in field
% $\texttt{ql\_info\_in}$ and is empty in this example because the call
% used default values. The structure that was actually used to process
% the file, in this case a structure of default values, is given in the
% field $\texttt{ql\_info}$.

% * 2015-10-31 RGL Update documentation
% * 2015-11-04 RGL Changed position of legend on Fig 1 for gliders.
%      Supressed spectragram outputs. Improved handling of pressure ranges
%      that are empty or do not exist. Changed CTD figure.
% * 2015-11-09 RGL Amended parameters associated with despike.
% * 2015-11-12 RGL Changed scaling on piezo-accelerometers for case of
%      Sea-glider. Changed profile_P_min to profile_min_P for consistency
%      with other profile parameters. Change profile_min_speed to
%      profile_min_W to make it more explictly the vertical component of
%      speed. Signicantly changed the way we handle accelerometers by
%      distinquishing between piezo and linear accelerometers, in a maner
%      that is backwards compatible with older data files that use
%      type=accel instead of type=piezo. Removed the ability to make pdf
%      fileas of the figures. Added option to supress making figures for
%      greater speed of execution.
% * 2015-11-18 RGL Force drawnow after every figure.
% * 2015-11-20 Changed recognisition of piezo-accelerometers to accomodate
%      changes to setupstr. Removed BP shear from figure for tidal_channel
%      data.
% * 2015-12-07 Check for existance of P_limit before using it. Required by
%      XMPs.
% * 2015-12-23 Fixed frequency spectra plot. Indices were erroneous so the 
%      wrong data range was being plotted.  Some minor cleaning of the
%      code.
% * 2016-01-21 Removed third linear accelerometer plot from Figure 1.
% * 2016-04-19 RGL, added factor of 10 scaling of mico-conductivity
%           gradient for case of "tidal_ch". Corrected legend bug --
%           numbers are no longer supported for location of legend. Fixed
%           bug with despiking the piezo-accelerometers. Now the statistics
%           are only for the profile and not the entire file.
% * 2016-04-25 RGL, added factor of 10 scaling of thermistor temperture
%           gradient for case of "tidal_ch". 

function result = quick_look( fname, P_start, P_end, varargin )


default_op_area       = 'open_ocean'; % Adjust some figure scales
default_make_figures  = true; % render figures for data visulization

% Default values used to determin profiles
default_profile_num          = 1; % process the first profile
default_profile_min_P        = 1; % in dBar
default_profile_min_W        = 0.2;
default_profile_min_duration = 20; % minimum duration [s] for a segment to be considered a profile

default_fft_length  = 2;% in units of seconds
default_diss_length = 4* default_fft_length;
default_overlap     = round(default_diss_length/2); % 50 percent over_lap
default_HP_cut      = 0.4; % high-pass shear probes at 0.5 Hz
default_LP_cut      = 30; % low-pass profile at 30 Hz, FOR PROFILE DISPLAY PURPOSES ONLY!
default_fit_2_isr   = 1.5e-5; % W/kg
default_fit_order   = 3;
default_f_AA        = 98; % anti-aliasing filter
default_f_limit     = inf;
default_YD_0        = 0;

% The despiking parameters are [thresh, smooth, and length (in seconds)
default_despike_sh  = [ 8  0.5 0.04]; % for shear probes
default_despike_A   = [ 8  0.5 0.04]; % for piezo-accelerometers
default_despike_C   = [10  1.0 0.04]; % for micro-C

if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    % Add the arguments from odas_p2mat.  They will not be required but can
    % be used as a reference when calling quick_look directly.
    p2mat = odas_p2mat();
    for name = fieldnames(p2mat)'
        result.(name{1}) = p2mat.(name{1});
    end
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_speed       = @(x) isnumeric(x) && isscalar(x)   || isempty(x);
val_vector      = @(x) (isnumeric(x) && isvector(x)) || isempty(x);
val_string      = @(x) ischar(x);
val_logical     = @(x) islogical(x);

addRequired(  p, 'fname',   val_string);
addRequired(  p, 'P_start', val_speed);
addRequired(  p, 'P_end',   val_speed);
addParamValue(p, 'fft_length',           default_fft_length,           val_numeric);
addParamValue(p, 'diss_length',          default_diss_length,          val_numeric);
addParamValue(p, 'overlap',              default_overlap,              val_numeric);
addParamValue(p, 'fit_order',            default_fit_order,            val_numeric);
addParamValue(p, 'HP_cut',               default_HP_cut,               val_numeric);
addParamValue(p, 'LP_cut',               default_LP_cut,               val_numeric);
addParamValue(p, 'profile_num',          default_profile_num,          val_numeric);
addParamValue(p, 'profile_min_P',        default_profile_min_P,        val_numeric);
addParamValue(p, 'profile_min_W',        default_profile_min_W,        val_numeric);
addParamValue(p, 'make_figures',         default_make_figures,         val_logical);
addParamValue(p, 'profile_min_duration', default_profile_min_duration, val_numeric);
addParamValue(p, 'fit_2_isr',            default_fit_2_isr,            val_numeric);
addParamValue(p, 'f_limit',              default_f_limit,              val_numeric);
addParamValue(p, 'f_AA',                 default_f_AA,                 val_numeric);
addParamValue(p, 'despike_sh',           default_despike_sh,           val_vector);
addParamValue(p, 'despike_A',            default_despike_A,            val_vector);
addParamValue(p, 'despike_C'  ,          default_despike_C,            val_vector);
addParamValue(p, 'YD_0',                 default_YD_0,                 val_numeric);
addParamValue(p, 'op_area',              default_op_area,              val_string);

ql_info_in = join_arguments(varargin); % save the input for record keeping

% Parse the arguments.
parse(p, fname, P_start, P_end, varargin{:});

% Perform last stages of input validation.
if p.Results.diss_length < 2*p.Results.fft_length,
  error('Invalid size for diss_length - must be greater than 2 * fft_length.');
end
if p.Results.P_end < p.Results.P_start,
  error('Starting pressure must be less than end pressure.');
end

names = fieldnames(p.Results);
for name = names'
    eval(['ql_info.' char(name) ' = p.Results.' char(name) ';']);
end
names = fieldnames(p.Unmatched);
for name = names'
    eval(['ql_info.' char(name) ' = p.Unmatched.' char(name) ';']);
end

fft_length           = p.Results.fft_length;
diss_length          = p.Results.diss_length;
overlap              = p.Results.overlap;
profile_min_P        = p.Results.profile_min_P;
profile_min_W        = p.Results.profile_min_W;
profile_num          = p.Results.profile_num;
HP_cut               = p.Results.HP_cut;
LP_cut               = p.Results.LP_cut;
profile_min_duration = p.Results.profile_min_duration;
fit_2_isr            = p.Results.fit_2_isr;
fit_order            = p.Results.fit_order;
f_AA                 = p.Results.f_AA;
f_limit              = p.Results.f_limit;
despike_sh           = p.Results.despike_sh;
despike_A            = p.Results.despike_A;
despike_C            = p.Results.despike_C;
YD_0                 = p.Results.YD_0;
op_area              = p.Results.op_area;
make_figures         = p.Results.make_figures;

YD_0 = floor(YD_0); % Force it to be w whole number

% Now save these parameter values so that they can be placed into the
% dissipation structure that is returned by this function, quick_look

%ql_info = p.Results;

d = odas_p2mat(fname, ql_info);
[P,N,E] = fileparts(d.fullPath);
File_Name = [N E];
if ~exist([P filesep N '.mat'], 'file')
    disp(['Saving into MAT-file: ' P filesep N '.mat']);
    save([P filesep N '.mat'], '-struct', 'd', '-v6');
else
    disp(['Loading from MAT-file: ' P filesep N '.mat']);
end

for field = fieldnames(d)'
    eval([char(field) ' = d.' char(field) ';']);
end

% The new odas_p2mat saves input parameters within a "params" structure. If
% a variable of the same name already exists, keep it.  This ensures that
% parameters passed into quick_look are used in place of those parameters
% used when generating the MAT file.
names = fieldnames(params);
for name = names'
    if ~exist(name{1}, 'var')
        eval([name{1} ' = params.' name{1} ';']);
    end
end

model       = char(setupstr(cfgobj,'instrument_info', 'model'));
serial_num  = char(setupstr(cfgobj,'instrument_info', 'serial_num'));
profile_dir = char(setupstr(cfgobj,'instrument_info', 'profile_dir'));
if isempty(profile_dir), profile_dir = vehicle_info.profile_dir; end

fft_num = round(fft_length*fs_fast);% fft length in units of samples

% Start of data processing and visulization
% ____________________________________
% extract a profile out of this file
profile = [1 ; length(t_slow)]; % start using entire file
% horizontal profilers get only a single profile, also if P_slow and
% W_slow do not exist
if exist('P_slow','var') && exist('W_slow','var')
        
    if strcmpi(profile_dir, 'up') || strcmpi(profile_dir, 'down')
        % get profile based on direction, duration, etc.
        profile = get_profile(P_slow, W_slow, profile_min_P, ...
                            profile_min_W, profile_dir, ...
                            profile_min_duration, fs_slow);
                            
    elseif strcmpi(profile_dir,'glide')
        profile_down = get_profile(P_slow, W_slow, profile_min_P, ...
                            profile_min_W, 'down', ...
                            profile_min_duration, fs_slow);
        profile_up   = get_profile(P_slow, W_slow, profile_min_P, ...
                            profile_min_W, 'up',   ...
                            profile_min_duration, fs_slow);
        % sort columns in ascending order
        profile = sort([profile_down profile_up],2);
    end
end

profiles_in_this_file = size(profile,2);
if profile_num > profiles_in_this_file
    warning(['There are only ' num2str(profiles_in_this_file) ' profiles in this file'])
    result = [];
    return
end

start_index_slow = profile(1, profile_num);
  end_index_slow = profile(2, profile_num);


start_index_fast = 1 + round((fs_fast/fs_slow)*(start_index_slow - 1));
  end_index_fast =     round((fs_fast/fs_slow)*(  end_index_slow    ));
n = (start_index_fast : end_index_fast)';
m = (start_index_slow : end_index_slow)';

% calculate a better estimate of the fall rate, than W_slow because W_slow
% is heavily biased by the zero speed after the profile. It was calculated
% using a 4th order zero-phase filter that "sees" the stop after the
% profile.

fc = 0.5;
[b,a]=butter(4,fc/(fs_slow/2));
dP_dt = gradient(P_slow(m),1/fs_slow);
dP_dt_LP = filtfilt(b,a,dP_dt);

% Here I trim the vectors to this profile to make the code look simpler.
profile_speed = speed_fast(n);

% The "hdr" variable contains the header of the starting profile - used to
% extract the time.
% NOTE: Time formatted in ISO 8601 format.
hd = header(ceil((start_index_slow - 1) / fs_slow + 1),:);
line1 = sprintf('%s;  %d-%02d-%02dT%02d:%02d:%02dZ', ...
                File_Name, hd(4), hd(5), hd(6), hd(7), hd(8), hd(9));
line2 = sprintf('Profile = %02d', profile_num);
title_string = texstr({line1, line2});

tidal_channel = strcmpi(op_area,'tidal_ch');

% Make figures
%
% define some colors
fuj  = [0.75 0    0.5]; % a shade of fujia
gold = [0.75 0.65 0];
% Pitch, Roll Fall-Rate
fig_num = 0;
LP_ADIS = 2; % low-pass filter applied to ADIS Inclinometer

% The names Incl_X and Incl_Y are standard for inclinometers
[b,a] =   butter(1,LP_ADIS/(fs_slow/2));% low-pass to show only the gravity signal
if make_figures && exist('Incl_X','var') && exist('Incl_Y','var') % then we have inclinometers
    fig_num = fig_num + 1;
    figure (fig_num); clf

    X_data = Incl_X;
    % check if inclinometer channels are reveresed
    if abs(mean(Incl_X(m))) > 50
        X_data = Incl_Y;
    end

    X_data = filtfilt(b,a,X_data);

    fs_ADIS = 482; % Internal sampling rate of the ADIS Inclinometer
    tau_N = 256; % The nuber of samples in the running average internal to the ADIS
    ADIS_delay = (tau_N/2) / fs_ADIS; % time delay of output from ADIS
    ADIS_delay = round(ADIS_delay * fs_slow); % time delay in samples

    mm = m;% new index for inclinometer
    if m(end) + ADIS_delay <= length(P_slow) % Check if we run past the end
        mm = m + ADIS_delay; % shift ADIS data forward in time to compensate its delay
    end
    MM = 2;
    if ~isempty(constant_speed),  MM = MM + 1;end
    if strcmpi(profile_dir,'glide'), MM = MM + 2; end

    subplot(1,MM,1)
    h = plot(X_data(mm), P_slow(m)); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('\theta_X [  ^{\circ} ]')
    title(title_string)
    P_limits = get(gca,'ylim');

    subplot(1,MM,2)
    h = plot([W_slow(m) dP_dt_LP], P_slow(m) ); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(h(2), 'linewidth', 2,   'color','r')
    set(gca,'ylim',P_limits)
    set(gca, 'ydir', 'rev')
%    set(gca,'yticklabel',[])
%    ylabel('\it P \rm [dBar]')
    xlabel('d\itP\rm / d\itt \rm[ dBar s^{-1}]')
    legend('W_1','W_2', 'location', 'southeast')

if ~isempty(constant_speed) && (strcmp(profile_dir,'up') || strcmp(profile_dir,'down'))
    subplot(1,MM,3)
    h = plot(constant_speed-dP_dt_LP, P_slow(m) ); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca,'ylim',P_limits)
    set(gca, 'ydir', 'rev')
    legend('W_2 - S_0','location','northwest')
    xlabel('w\prime [ m s^{-1}]')
end

if  strcmpi(profile_dir,'glide')
    subplot(1,MM,3)
    h = plot(speed_slow(m), P_slow(m) ); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca,'ylim',P_limits)
    set(gca, 'ydir', 'rev')
%    set(gca,'yticklabel',[])
%    ylabel('\it P \rm [dBar]')
    xlabel('speed [ m s^{-1}]')

    subplot(1,MM,4)
    h = plot([Incl_X(m) Incl_Y(m)], P_slow(m) ); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(h(2), 'linewidth', 2,   'color','r')
    set(gca,'ylim',P_limits)
    set(gca, 'ydir', 'rev')
    legend('\theta_X', '\theta_Y', 'location','northeast')
    xlabel('[  ^{\circ} ]')
end
drawnow

end

% Identify if we have the true linear (DC) response accelerometers, or if
% we have piezo-accelerometers. A lot of data has been collected using type
% = accel for piezo acceleromters, with coefficients of 0 and 1. On
% 2015-11-12 we introduced the type=piezo to make it easier to handle both
% types of accelerometers, especially with respect to scaling the spectra of
% acceleration. This used to be a nightmare of if-statements. Here we
% figure out what type of accelerometers we have (including the possibility
% of both types, such as with the Nemo system), and place the data into
% matrices with appropriate names. In addition, we place AA_piezo into AA
% for coherent-noise removal. If the are no piezo acceleromters, then we
% place the linear accelerometers, AA_linear, into the matrix AA.

AA_linear       = [];
AA_names_linear = {};
AA_piezo        = [];
AA_names_piezo  = {};
for ch = setupstr(cfgobj, '', 'type', 'accel|piezo')
    name = setupstr(cfgobj, ch{1}, 'name');
    type = setupstr(cfgobj, ch{1}, 'type');
    
    if strcmp(type, 'accel')
        AA_names_linear{end+1} = name{1};
        AA_linear = [AA_linear eval(name{1})];
    else
        AA_names_piezo{end+1} = name{1};
        AA_piezo = [AA_piezo eval(name{1})];
    end
end

piezo = ~isempty(AA_piezo);

AA = [];
if piezo
    AA = AA_piezo;
    AA_piezo  = AA_piezo (n,:);
elseif ~isempty(AA_linear)
    AA = AA_linear;
    AA_linear = AA_linear(n,:);
end
%-----------------------------
if ~isempty(AA_linear) && make_figures
% Then we have old-fashioned linear DC response accelerometers. Plot the
% pitch and roll
fig_num = fig_num + 1;
    figure (fig_num); clf
    f_low_pass = 1; % in Hz
    [b,a] = butter(2, f_low_pass /(fs_fast/2));

    X_data = filtfilt(b,a,AA_linear(n,1:min(2,end))); % because we want rotation around y-axis
    X_data = asind(X_data / 9.81); % in degrees of angle

    index_2_Ax = find(strcmpi('Ax', AA_names_linear));
    if isempty (index_2_Ax), error('Cannot find channel with name = Ax'), end
    X_data (:,index_2_Ax) = -X_data(:,index_2_Ax);
    
    subplot(1,2,1)
    h = plot(X_data, P_fast(n)); grid on
    if ~exist('P_limits','var')
        P_limits = get(gca,'ylim');
    end
    set(gca,'ylim',P_limits)
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('[  ^{\circ} ]')
    legend('\theta_Y', '\theta_X',4)
    title(title_string)

    subplot(1,2,2)
    h = plot(W_slow(m), P_slow(m) ); grid on
    set(gca,'ylim',P_limits)
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    set(gca,'yticklabel',[])
    xlabel('d\it P \rm/ d\itt \rm[  dBar s^{-1}]')
end

drawnow

% for an XMP, only
if exist('Pitch','var') && make_figures
    fig_num = fig_num + 1;
    figure (fig_num); clf

    X_data = Pitch;
    X_data = filtfilt(b,a,X_data);

    subplot(1,2,1)
    h = plot(X_data(m), P_slow(m)); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('Pitch [  ^{\circ} ]')
    title(title_string)

    subplot(1,2,2)
    h = plot(W_slow(m), P_slow(m) ); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    set(gca,'yticklabel',[])
    xlabel('\it W \rm[  m s^{-1}]')
end

drawnow

%
%-----------------------------------------------------------------------
% Now we plot the profile of acceleration.
%
% Piezo accelerometers before and after despiking.
if ~isempty(AA_piezo)
    spikes_A     = {}; % pre-assign
    pass_count_A = zeros(size(AA_piezo,2),1);
    fraction_A   = zeros(size(AA_piezo,2),1);
    
    for index = 1:2 % first round no despiking, second round with despiking
        [b,a] =   butter(1,1/(fs_fast/2));% low-pass to show only the gravity signal
        piezo_accel_num = size(AA_piezo,2);

        if index == 1
            for junk = 1:piezo_accel_num
                legend_string{junk}                   =  AA_names_piezo{junk};
                legend_string{piezo_accel_num + junk} = [AA_names_piezo{junk} ' LP'];
            end
        end

        if index == 2 % despike the piezo-accelerometer signals
            if despike_A(1) ~= inf
                for probe = 1:piezo_accel_num
                    [AA_piezo(:,probe), spikes_A{probe}, ...
                        pass_count_A(probe), fraction_A(probe)]  = ...
                        despike(AA_piezo(:,probe),  despike_A(1), ...
                        despike_A(2), fs_fast, round(despike_A(3)*fs_fast));
                end
            end
        end
        
        if make_figures
           fig_num = fig_num + 1;
           figure (fig_num); clf

           plot([AA_piezo filtfilt(b, a, AA_piezo)], P_fast(n));grid on
           if exist('P_limits', 'var')
               set(gca,'ylim',P_limits)
           end
           set(gca, 'ydir', 'rev')
           xlabel('[counts]')
           ylabel(' \itP\rm [dbar]')
           legend(legend_string,'location','NorthEastOutside');
           
           new_title_string = title_string;
           if ~isempty(spikes_A)
               despike_string = [];
               for accel_index = 1:length(spikes_A)
                   despike_string = [despike_string num2str(length(spikes_A{accel_index})) ' '];
               end
               new_title_string{end+1} = ...
                   ['with despiking = ' despike_string];
           end
           %                ['with despiking = ' num2str([length(spikes_A{1}) ...
           %                length(spikes_A{2})])]; ...
           legend(legend_string,'location','NorthEastOutside');
           title(new_title_string)
           if ~exist('P_limits', 'var')
               P_limits = get(gca,'ylim');
           end
        end
    end
end

if ~isempty(AA_linear) && make_figures
    % Then we have a 3-axis linear DC response accelerometer, such as used
    % in older instruments.
    fig_num = fig_num + 1;
    figure (fig_num); clf
    
    plot(AA_linear(n,:), P_fast(n));grid on
    set(gca,'ylim',P_limits)
    set(gca, 'ydir', 'rev')
    xlabel('[m s^{-1}]')
    ylabel(' \itP\rm [dbar]')
    legend(AA_names_linear,'location','NorthEastOutside');

    title(title_string)
end
drawnow

%
% High accuracy CTD, if available
if make_figures
we_have_T = false;
we_have_C = false;
we_have_S = false;
T_CTD = [];
C_CTD = [];
S_CTD = [];
sbt_offset = 0.35;
JAC_offset = 0.14;
legend_string = [];

if exist('SBT1','var')
    we_have_T = true;
    T_CTD = SBT1;
    offset = sbt_offset;
    legend_string{end+1} = 'SBT1';
elseif exist('SBT','var')
    we_have_T = true;
    T_CTD = SBT;
    offset = sbt_offset;
    legend_string{end+1} = 'SBT';
elseif exist('sbt','var')
    we_have_T = true;
    T_CTD = sbt;
    offset = sbt_offset;
    legend_string{end+1} = 'sbt';
elseif exist('JAC_T','var')
    we_have_T = true;
    T_CTD = JAC_T;
    offset = JAC_offset;
    legend_string{end+1} = 'JAC\_T';
end
if exist('SBC1','var')
    we_have_C = true;
    C_CTD = SBC1;
    offset = sbt_offset;
    legend_string{end+1} = 'SBC1';
elseif exist('SBC','var')
    we_have_C = true;
    C_CTD = SBC;
    offset = sbt_offset;
    legend_string{end+1} = 'SBC';
elseif exist('sbc','var')
    we_have_C = true;
    C_CTD = sbc;
    offset = sbt_offset;
    legend_string{end+1} = 'sbc';
elseif exist('JAC_C','var')
    we_have_C = true;
    C_CTD = JAC_C;
    offset = JAC_offset;
    legend_string{end+1} = 'JAC\_C';
end

windows = 0;
if we_have_T
    windows = windows + 1;
    T_CTD = T_CTD(m);
end
if we_have_C
    windows = windows + 1;
    C_CTD = C_CTD(m);
end
if we_have_T && we_have_C
    we_have_S = true;
    windows = windows + 1;
    S_CTD = salinity(P_slow(m),T_CTD,C_CTD);
    legend_string{end+1} = 'S';
end

if windows > 0
    fig_num = fig_num +1;
    figure (fig_num); clf
    plot_data = [T_CTD C_CTD S_CTD];
    for index = 1:windows
        subplot(1,windows,index)
        plot(plot_data(:,index), P_slow(m) + offset);grid on
        if index == 1, ylabel('\it P \rm[ dbar ]'), end
        legend(legend_string{index},'location','northwest')
        set(gca, 'ydir', 'rev')
    end
    title(title_string)

end
end
drawnow
% Thermistors
if ...
        make_figures           && ...
       (exist('T1_fast','var') || ...
        exist('T2_fast','var') || ...
        exist('T_fast','var'))
    fig_num = fig_num +1;
    figure (fig_num); clf
    legend_string = [];
    plot_data     = [];
    P_data_fast   = [];
    T_count       = 0;
    if exist('T1_fast','var')
        T_count = T_count + 1;
        legend_string{T_count} = '\itT\rm_1';
        plot_data = [plot_data T1_fast(n)];
        P_data_fast =[P_data_fast P_fast(n)];
    end
    if exist('T2_fast','var')
        T_count = T_count + 1;
        plot_data = [plot_data T2_fast(n)];
        legend_string{T_count} = '\itT\rm_2';
        P_data_fast =[P_data_fast P_fast(n)];
    end
    if exist('T_fast','var')
        T_count = T_count + 1;
        plot_data = [plot_data T_fast(n)];
        P_data_fast =[P_data_fast P_fast(n)];
        legend_string{T_count} = '\itT';
    end
    P_data_slow    = [];
    slow_data = [];
    if exist('SBT1','var') || exist('SBT','var') || exist('sbt','var')
        offset = 0.34; % typical offset of SBT4, in metres
        T_count = T_count + 1;
        P_data_slow = [P_data_slow P_slow(m) - offset];
        if exist('SBT1','var'), slow_data = [slow_data SBT1(m)];end
        if exist('SBT', 'var'), slow_data = [slow_data SBT(m)];end
        if exist('sbt', 'var'), slow_data = [slow_data sbt(m)];end
        legend_string{T_count} = '\itT_{SB}';
    end

    if exist('JAC_T','var')
        offset = 0.14; % in metres
        T_count = T_count + 1;
        P_data_slow = [P_data_slow P_slow(m) - offset];
        slow_data = [slow_data JAC_T(m)];
        legend_string{T_count} = '\itT_{JAC}';
    end
    h = plot(plot_data, P_data_fast, slow_data, P_data_slow);grid on
    set(gca,'ylim',P_limits)
    legend(legend_string,  'location','NorthEastOutside')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('[  ^{\circ}C ]')
    title(title_string)
end
drawnow
%
%----------------------------------------------------------
%micro-conductivity sensors
if make_figures  && (exist('C1_fast','var') || exist('C2_fast','var'))
    fig_num = fig_num +1;
    figure (fig_num); clf
    legend_string = [];
    plot_data     = [];
    C_count       = 0;
    if exist('C1_fast','var')
        C_count = C_count + 1;
        legend_string{C_count} = '\itC\rm_1';
        plot_data = [plot_data C1_fast(n)];
    end

    if exist('C2_fast','var')
        C_count = C_count + 1;
        plot_data = [plot_data C2_fast(n)];
        legend_string{C_count} = '\itT\rm_2';
    end
    P_data    = [];
    slow_data = [];

    if exist('JAC_C','var')
        C_count = C_count + 1;
        P_data = P_slow(m) - offset;
        slow_data = JAC_C(m);
        legend_string{C_count} = '\itC_{JAC}';
    end


    if exist('SBC1','var') || exist('SBC','var') || exist('sbc','var')
        offset = 0.34; % typical offset of SBT4, in metres
        C_count = C_count + 1;
        P_data = P_slow(m) - offset;
        if exist('SBC1','var'), slow_data = SBC1(m);end
        if exist('SBC', 'var'), slow_data = SBC (m);end
        if exist('sbc', 'var'), slow_data = sbc (m);end
        legend_string{C_count} = '\itC_{SB}';
    end

    h = plot(plot_data, P_fast(n), slow_data, P_data);grid on
    set(gca,'ylim',P_limits)
    legend(legend_string,  'location','NorthEastOutside')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('[ mS cm^{-1} ]')
    title(title_string)
end
drawnow
%
%----------------------------------------------------------
% Fluorometer and Backscatter sensor. Must be in fast channels.
% Suna must be a slow channel.
if ...
        make_figures               && ...
       (exist ('Fluo','var')       ||...
        exist('Turbidity','var')   ||...
        exist('Chlorophyll','var') ||...
        exist('BS','var')          || ...
        exist('Suna','var'))
    fig_num = fig_num +1;
    figure (fig_num); clf
    legend_string = [];
    Y_fast        = [];
    Y_slow        = [];
    if exist('Fluo','var')
        legend_string{end+1} = 'Fluoro';
        Y_fast = [Y_fast Fluo(n)];
    end
    if exist('BS','var')
        Y_fast = [Y_fast BS(n)];
        legend_string{end+1} = 'BS';
    end
    if exist('Turbidity','var')
        Y_fast = [Y_fast Turbidity(n)];
        legend_string{end+1} = 'Turbidity';
    end
    if exist('Chlorophyll','var')
        Y_fast = [Y_fast Chlorophyll(n)];
        legend_string{end+1} = 'Chlorophyll';
    end
    if exist('Suna','var'),
        Y_slow = [Y_slow Suna(m)];
        legend_string{end+1} = 'Suna';
    end
    P_junk_fast = [];
    P_junk_slow = [];
    if ~isempty(Y_fast), P_junk_fast = P_fast(n); end
    if ~isempty(Y_slow), P_junk_slow = P_slow(m); end

    h = plot(Y_fast, P_junk_fast, Y_slow, P_junk_slow);grid on
    set(gca,'ylim',P_limits)
    legend(legend_string,  'location','NorthEastOutside')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('[ ppb ]  [ FTU ]')
    title(title_string)
end
drawnow

% Battery Voltage
if exist('V_Bat','var') && make_figures
    fig_num = fig_num + 1;
    figure(fig_num); clf
    [b,a]=butter(1,1/(fs_slow/2));
    h = plot([V_Bat(m) filtfilt(b,a,V_Bat(m))], P_slow(m));grid on
    set(gca,'ylim',P_limits)
    set(h(2),'linewidth',3,'color','m')
    set(gca, 'ydir', 'rev')
    ylabel('\it P \rm [dBar]')
    xlabel('\itV_{Bat}\rm [volts]')
    title(title_string)
    legend('V_{Bat}','B_{Bat}-LP','location','eastoutside')
end
drawnow

%
%--------------------------------------------------------------------------
% Temperature gradient, microconductivity gradient, and shear

shear_wish_list = {'sh1', 'sh2', 'sh3', 'sh4'};
shear_list  = {};
shear_count = 0;
SH          = [];
plot_gradient_list = {};
legend_list = {};
for index = 1:length(shear_wish_list)
    if exist(shear_wish_list{index},'var')
        shear_count = shear_count + 1;
        shear_list{end+1} = shear_wish_list{index};
        plot_gradient_list{end+1} = [shear_wish_list{index} '_HP'];
        legend_list{end+1} = [shear_wish_list{index} ' HP'];
        junk = eval(shear_wish_list{index});
        SH = [SH junk(n)];
    end
end

% Now a loop for the shear band-pass data SH_BP, yet to be created
if ~tidal_channel
    for index = 1:length(shear_wish_list)
        if exist(shear_wish_list{index},'var')
            plot_gradient_list{end+1} = [shear_wish_list{index} '_BP'];
            legend_list{end+1} = [shear_wish_list{index} ' BP'];
        end
    end
end

T_scale = 1;
T_scale_string = '';
if tidal_channel
    T_scale = 0.1;
    T_scale_string = '0.1\times';
end

gradT_wish_list = {'gradT', 'gradT1', 'gradT2'};
gradT_list = {};
T_count = 0;
GRADT = [];
for index = 1:length(gradT_wish_list)
    if exist(gradT_wish_list{index},'var')
        T_count = T_count + 1;
        gradT_list{end+1} = gradT_wish_list{index};
        plot_gradient_list{end+1} = [gradT_wish_list{index} '_LP'];
        legend_list{end+1} = [T_scale_string '\partialT_' num2str(T_count) '/\partialz LP'];
        junk = eval(gradT_wish_list{index});
        GRADT = [GRADT T_scale*junk(n)];
    end
end

C_scale = 1;
C_scale_string = '';
if tidal_channel
    C_scale = 0.1;
    C_scale_string = '0.1\times';
end
    
gradC_wish_list = {'gradC', 'gradC1', 'gradC2'};
gradC_list = {};
C_count = 0;
GRADC = [];
for index = 1:length(gradC_wish_list)
    if exist(gradC_wish_list{index},'var')
        C_count = C_count + 1;
        gradC_list{end+1} = gradC_wish_list{index};
        plot_gradient_list{end+1} = gradC_wish_list{index};
        legend_list{end+1} = [C_scale_string '\partialC_' num2str(C_count) '/\partialz'];
        junk = eval(gradC_wish_list{index});
        GRADC = [GRADC C_scale*junk(n)];
    end
end

num_of_C = C_count; % the number of micro-conductivity signals

C_count = 0;
for index = 1:length(gradC_list)
        C_count = C_count + 1;
        plot_gradient_list{end+1} = [gradC_wish_list{index} '_LP'];
        legend_list{end+1} = [C_scale_string '\partialC_' num2str(C_count) '/\partialz LP'];
end

% Despike the shear signals
spikes_sh     = {}; % they must exist
pass_count_sh = zeros(shear_count,1);
fraction_sh   = zeros(shear_count,1);

if despike_sh(1) ~= inf
    for index = 1:shear_count
        [SH(:,index), spikes_sh{index}, pass_count_sh(index), fraction_sh(index) ] =  ...
            despike(SH(:,index), despike_sh(1), despike_sh(2), fs_fast, round(despike_sh(3)*fs_fast));
    end
end


% Despike the micro-conductivity signals.
spikes_C = {}; % they must exist
pass_count_C = zeros(num_of_C,1);
fraction_C   = zeros(num_of_C,1);

if despike_C(1) ~= inf
    for index = 1:num_of_C
        [GRADC(:,index), spikes_C{index}, pass_count_C(index), fraction_C(index)] =  ...
            despike(GRADC(:,index), despike_C(1), despike_C(2), fs_fast, round(despike_C(3)*fs_fast));
    end
end

% Apply gentle high-pass filter
[bh,ah] = butter(1, HP_cut/(fs_fast/2), 'high'); % HP at ~1 cpm
[bl,al] = butter(4, LP_cut/(fs_fast/2));% for profile plotting only

% High-pass the shear data with sensible initial conditions.
SH_HP = filtfilt(bh, ah, SH);
% Low-pass the shear data with sensible initial conditions.
SH_BP = filtfilt(bl, al, SH_HP);

% Low-pass filter the scalar gradients.
GRADT_LP = [];
GRADC_LP = [];
if T_count>0, GRADT_LP = filtfilt(bl, al, GRADT);end
if C_count>0, GRADC_LP = filtfilt(bl, al, GRADC);end

if make_figures
fig_num = fig_num + 1;
figure(fig_num); clf

%plot_gradient_list = [];

plot_scalar_wish_list = {'T_fast', 'T1_fast', 'T2_fast', 'C1_fast', 'C2_fast'};
legend_scalar_wish_list = {'T', 'T_1','T_2', 'C_1','C_2' };
plot_scalar_list = {};

x_gradient_step   = 5;
if tidal_channel, x_gradient_step = 20; end
left_edge = -x_gradient_step;

x_gradient_offset = left_edge + x_gradient_step;

LL = length(plot_gradient_list);
%index_gradient; % the number of gradient channels

for k = 1:length(plot_scalar_wish_list)
    if exist(plot_scalar_wish_list{k},'var')
        plot_scalar_list{end+1} = plot_scalar_wish_list{k};
        legend_list{end+1} = legend_scalar_wish_list{k};
    end
end
num_of_scalars = length(plot_scalar_list);

% assume all are fast channels with same sampling rate
X_scalar_data = zeros(length(n),num_of_scalars);
x_scalar_offset = 0;

if tidal_channel
    X_gradient_data = [SH_HP       [GRADT_LP GRADC GRADC_LP]];
else
    X_gradient_data = [SH_HP SH_BP [GRADT_LP GRADC GRADC_LP]];
end
for k = 1:LL
%    eval(['X_gradient_data(:,k) = '  plot_gradient_list{k} '(n);'])
   X_gradient_data(:,k) = X_gradient_data(:,k) + x_gradient_offset;
   x_gradient_offset = x_gradient_offset + x_gradient_step;
end

for k = 1: num_of_scalars
   eval(['X_scalar_data(:,k) = '  plot_scalar_list{k} '(n);'])
   X_scalar_data(:,k) = X_scalar_data(:,k) + x_scalar_offset;
   x_scalar_offset = x_scalar_offset + 0;
end

X_data = [];
if ~isempty(X_gradient_data)
    X_data = X_gradient_data;
end
if ~isempty(X_scalar_data)
    X_data = [X_data X_scalar_data];
end

plot(X_data, P_fast(n) ); grid on
set(gca,'ylim',P_limits)
set(gca,'ydir','rev')
x_lim = get(gca, 'xlim');
x_lim = [left_edge x_lim(2)];
if tidal_channel, x_lim = [-20 x_gradient_offset]; end
set(gca,'xlim',x_lim)
legend(legend_list, 'location','eastoutside')
ylabel('\itP\rm [dBar]')

if tidal_channel
    junk = ['HP = ' num2str(HP_cut) ' Hz'];
else
    junk = ['BP = ' num2str(HP_cut) ' \cdot\cdot\cdot ' num2str(LP_cut) ' Hz'];
end
title([title_string junk])
drawnow

%
%----------------------------------------------------------------
% Frequency spectra
fig_num = fig_num + 1;
figure(fig_num); clf

if isempty(P_start) || isempty(P_end)
    range = [];
else
    range = find((P_fast(n) > P_start) & (P_fast(n) <= P_end)) + (n(1)-1);
    if isempty(range)
        warning(['The pressure range of ' num2str(P_start) ' to ' ...
            num2str(P_end) ' was not found in this profile'])
    end
end

if isempty(range)
    fig_num = fig_num + 1;
    figure(fig_num); clf
else

P_start = round(P_fast(range(1))); % So the labels give the actual range of data used
P_end   = round(P_fast(range(end)));

Data_fast   = [];
Data_slow = [];

info.fft_length    = fft_num;
info.diss_length   = length(range);
info.overlap       = info.diss_length/2;
info.fs_fast       = fs_fast;
info.fs_slow       = fs_slow;
info.speed         =       speed_fast(range);
info.T             = temperature_fast(range);
info.t             =           t_fast(range);
info.P             =           P_fast(range);
info.Data_fast     = Data_fast;
info.Data_slow     = Data_slow;
info.fit_order     = fit_order;
info.f_AA          = f_AA;
info.f_limit       = f_limit;
info.fit_2_isr     = fit_2_isr;

% range-(n(1)-1) ==== data points with respect to the start of profile.
diss = get_diss_odas(SH_HP(range-(n(1)-1),:), AA(range,:), info);

mean_speed = diss.speed;
F = diss.F;
K = diss.K;

% get_diss_odas returns wavenumber spectra. Convert to frequency spectra by
% dividing by the mean speed. The spectra are corrected for the wavenumber
% response of the shear probe.
for index = 1:shear_count
    r = num2str(index); % put it into a string
    junk = ['P_' shear_list{index}]; % The name of the spectrum
    eval([junk '       = diss.sh      (:,' r ',' r ') / mean_speed;' ])
    eval([junk '_clean = diss.sh_clean(:,' r ',' r ') / mean_speed;' ])
end

if piezo 
    % We have piezo-accelerometers so scale-down the acceleration spectra
    % to bring them in line with the shear spectra. 
    diss.AA = 1e-6*diss.AA;
end

accel_count = size(diss.AA,2);
if accel_count > 0, P_Ax = diss.AA(:,1,1);end
if accel_count > 1, P_Ay = diss.AA(:,2,2);end
if accel_count > 2, P_Az = diss.AA(:,3,3);end

% produce spectra of scalar variable, if any.
scalar_vector_list = [];
scalar_vectors     = [];
diff_gain          = [];

%                   'grad' prepended to each name
scalar_wish_list = {'T', 'T1', 'T2', 'C', 'C1', 'C2'};

for index = 1:length(scalar_wish_list)
    scalar_name = ['grad' scalar_wish_list{index}];
    vector_name = [scalar_wish_list{index} '_d' scalar_wish_list{index}];
    if exist(scalar_name,'var')
        scalar_vector_list{end+1} = scalar_name;
        scalar_vectors = [scalar_vectors eval(scalar_name)];
        tmp = setupstr(cfgobj, vector_name, 'diff_gain');  % Find diff_gain
        diff_gain = [diff_gain str2double(tmp{1})];% Assume it is numeric
    end
end

scalar_info.fft_length      = fft_num;
scalar_info.spec_length     = length(range);
scalar_info.overlap         = info.overlap;
scalar_info.fs              = fs_fast;
scalar_info.gradient_method = 'first_difference';
scalar_info.diff_gain       = diff_gain;
scalar_info.f_AA            = 98;

if ~isempty(scalar_vector_list)
    scalar_spectra = get_scalar_spectra_odas(scalar_vectors(range,:), ...
        info.P, info.t, info.speed, scalar_info);

    for index = 1:length(scalar_vector_list)
        % extract the scalar spectra
        name = ['P_' scalar_vector_list{index}];
        spectra = scalar_spectra.scalar_spec(:,index) / scalar_spectra.speed;
        eval([name ' = spectra;']);
    end
end

epsilon = 1e-10*[1 10 100 1e3 1e4 1e5];
if tidal_channel, epsilon = 10*epsilon; end
nu = visc35(mean(temperature_fast(range)));
[Pn, kn]=nasmyth(epsilon, nu);
fn = kn*mean_speed; Pn = Pn/mean_speed;

plot_wish_list = {...
    'P_Ax', 'P_Ay', 'P_Az', ...
    'P_sh1', 'P_sh2', 'P_sh3', 'Psh4', ...
    'P_gradT', 'P_gradT1', 'P_gradT2', ...
    'P_gradC1', 'P_gradC2'};
legend_wish_list = {...
    'A_x', 'A_y', 'A_z', ...
    'sh_1', 'sh_2', 'sh_3', 'sh_4', ...
    '\partialT/\partialz', '\partialT_1/\partialz', '\partialT_2/\partialz', ...
    '\partialC_1/\partialz', '\partialC_2/\partialz'};

index = 0;
plot_list = [];
legend_list = [];
for k = 1:length(plot_wish_list)
    if exist(plot_wish_list{k},'var')
        index = index + 1;
        plot_list{index}   = plot_wish_list{k};
        legend_list{index} = legend_wish_list{k};
    end
end
legend_list{index+1} = 'Nasmyth';

num_vectors = index;
Y_data = zeros(size(F,1),num_vectors);
for k = 1: num_vectors
   eval(['Y_data(:,k) = '  plot_list{k} ';'])
end

h = loglog(F, Y_data, fn, Pn, 'k');grid on
hh=legend(legend_list,'location','NorthEastOutside');

for index = 1: size(Y_data,2)
    set(h(index), 'linewidth',1.5)
end

set(h(1),'linewidth',3,'color','b'); % Ax
set(h(2),'linewidth',3,'color',[0 0.5 0]); % Ay
if exist('Az','var'), set(h(3),'linewidth',3 ,'color','r');end % Az

f_upper = max(F);

y_limit = [1e-9 10]; % limits for y-axis
x_limit = [1./fft_length f_upper];% limits for x-axis
set(gca, 'ylim', y_limit, 'xlim',x_limit)

for index = 1: length (epsilon)
    f_location = find(fn(:,index) > x_limit(1));
    if isempty(f_location)
        f_location = 1;
    else
        f_location = f_location(1);
    end
    h=text(fn(f_location,index), Pn(f_location,index), ...
        ['10^{' num2str(log10(epsilon(index)),2) '}'], ...
         'fontweight','bold');
end
ylabel('[Variance Hz^{-1}]')
xlabel('\it f \rm [Hz]')
new_title_string = title_string;
new_title_string{size(new_title_string,2)+1} = ...
    [ num2str(P_start) ' < P < ' num2str(P_end) ' m, f_{HP} = ' num2str(HP_cut) ' Hz'];
title(new_title_string)
drawnow

%
%_____________________________
% Wavenumber Spectra
fig_num = fig_num + 1;
figure(fig_num); clf

% Extract wavenumber spectra of shear
for index = 1:shear_count
    r = num2str(index); % put it into a string
    eval(['P_sh' r       ' = diss.sh(:,'       r ',' r ');'])
    eval(['P_sh' r '_clean = diss.sh_clean(:,' r ',' r ');'])

end

F_0 = 25*sqrt(mean_speed); % frequency response of thermistor after Vachon and Lueck

% compensate frequency spectra of thermistors
thermistor_list = {'P_gradT', 'P_gradT1', 'P_gradT2'};
% These are still frequency spectra

for index = 1:length(thermistor_list)
    if exist(thermistor_list{index},'var')
       eval([thermistor_list{index} ' = ' ...
             thermistor_list{index} '.* (1 + (F/F_0).^2);'])
    end
end

% Convert all scalar_spectra to wavenumber spectra
for index = 1:length(scalar_vector_list)
    junk = ['P_' scalar_vector_list{index}];
    eval([junk ' = ' junk '.* mean_speed;'])
end

e1_string = []; % They must exist
e2_string = [];
e3_string = [];
e4_string = [];

% find index to the limit of integration.
K_max = diss.K_max;
nu = diss.nu;
e = diss.e;
for index = 1:shear_count
    r = num2str(index); % Put it into a string
    eval(['e' r ' = e(' r ');'])
    eval(['K_max_index_' r ' = find (K == K_max(' r '));'])
    eval(['phi' r   ' = diss.Nasmyth_spec(:,' r ');'])
    e_value = eval(['e' r ';']);
    junk = ['\epsilon=' make_scientific(e_value,2) 'W kg^{-1}'];
    eval(['e' r '_string =  junk ;'])

end

plot_wish_list = {
    % Plot settings for all channels.
    % Name          Legend Title                        Colour  Size
    {'P_sh1_clean', '\partial u_1/\partial z_{clean}',  'b',    3   }
    {'P_sh1',       '\partial u_1/\partial z',          'b',    1.5 }
    {'phi1',        e1_string,                          'k',    1.5 }
    {'P_sh2_clean', '\partial u_2/\partial z_{clean}',  'r',    3   }
    {'P_sh2',       '\partial u_2/\partial z',          'r',    1.5 }
    {'phi2',        e2_string,                          'k',    1.5 }
    {'P_sh3_clean', '\partial u_3/\partial z_{clean}', gold,    3   }
    {'P_sh3',       '\partial u_3/\partial z',         gold,    1.5 }
    {'phi3',        e3_string,                          'k',    1.5 }
    {'P_sh4_clean', '\partial u_4/\partial z_{clean}',  fuj,    3   }
    {'P_sh4',       '\partial u_4/\partial z',          fuj,    1.5 }
    {'phi4',        e4_string,                          'k',    1.5 }
    {'P_gradT',     '\partial T/\partial z ',           'g',    2   }
    {'P_gradT1',    '\partial T_1/\partial z ',         'g',    2   }
    {'P_gradT2',    '\partial T_2/\partial z ',         'm',    2   }
    {'P_gradC1',    '\partial C_1/\partial z ',         'c',    2   }
    {'P_gradC2',    '\partial C_2/\partial z ',         'c',    2   }};

point_wish_list = {
    % X pos         Y pos                           Label                                                 Colour  Size
    {'K_max(1)',    'P_sh1_clean(K_max_index_1)',   '[''k_{max} u_1='' num2str(round(K_max(1))) ''cpm'']',  'b',    18  }
    {'K_max(2)',    'P_sh2_clean(K_max_index_2)',   '[''k_{max} u_2='' num2str(round(K_max(2))) ''cpm'']',  'r',    18  }
    {'K_max(3)',    'P_sh3_clean(K_max_index_3)',   '[''k_{max} u_3='' num2str(round(K_max(3))) ''cpm'']', gold,    18  }
    {'K_max(4)',    'P_sh4_clean(K_max_index_4)',   '[''k_{max} u_4='' num2str(round(K_max(4))) ''cpm'']',  fuj,    18  }};

% Remove all channels that are not available from the wish list.
plot_list = {};
for ch = plot_wish_list'
    if exist(ch{1}{1}, 'var') plot_list{end+1,1} = ch{1}; end
end

% Remove all points that are not available from the wish list.  Evaluate
% and replace the values.
point_list = {};
for pt = point_wish_list'
    try
        point_list{end+1,1} = {eval(pt{1}{1}), eval(pt{1}{2}), eval(pt{1}{3}), pt{1}{4}, pt{1}{5}};
    catch continue; end
end

% Extract the required data then plot on a loglog plot.
Y_data = [];
for ch = plot_list', Y_data(:,end+1) = eval( ch{1}{1} ); end
log_plot = loglog(K, Y_data);

% Set the requested line width and colour settings.
for ch = 1:size(plot_list,1)
    set(log_plot(ch),'linewidth',plot_list{ch}{4},'color',plot_list{ch}{3});
end

% Plot the requested points with correct colours.
hold on;
for pt = point_list'
    loglog( pt{1}{1}, pt{1}{2},         ...
        'Marker',           '^',        ...
        'MarkerSize',       pt{1}{5},   ...
        'MarkerFaceColor',  pt{1}{4},   ...
        'MarkerEdgeColor',  'w');
end
hold off
grid on;

set(gca, 'ylim',[1e-8 1e0], 'xlim', [K(2) K(end)])

if tidal_channel, set(gca,'ylim',[1e-7 1e1]), end

xlabel('\itk \rm [cpm]')
ylabel('\Phi (\itk\rm)  [s^{-2} cpm^{-1}]')

legend_list = {};
for plt = plot_list',  legend_list{end+1} = plt{1}{2}; end
for plt = point_list', legend_list{end+1} = plt{1}{3}; end
legend(legend_list,'location','NorthEastOutside');

title_string_2_old = title_string{2};
title_string{2} = [title_string{2} ...
    ',   Method = ' num2str(diss.method') ];

new_title_string = [title_string , ...
        sprintf('P = %d - %d m, W = %0.2f m s^{-1}, f_{HP} = %0.2f Hz', ...
        P_start, P_end, mean_speed, HP_cut)];

title(new_title_string)
title_string{2} = title_string_2_old;

end
end
drawnow

%
% Now we calculate a dissipation profile

Data_fast = []; % make empty so that it does not get counted in the search for vectors.
Data_slow = [];

% Develop the lists of fast and slow vectors in this data file. We will
% only recognize column vectors of the right length.
junk = [];
S = whos; % The structure containing all of the variables:
% name       -- variable name
% size       -- variable size
% other fields can be ignored.

junk = []; for k = 1:size(S,1), junk = [junk; S(k).size]; end
index_fast = find((junk(:,1) == size(t_fast,1)) & junk(:,2) == 1); % index_fast points to the fast channels.
index_slow = find((junk(:,1) == size(t_slow,1)) & junk(:,2) == 1); % index_slow points to the slow channels.
% Move vectors into a matrix to pass to the get_diss_odas function, where
%  they will be averaged over the range of each dissipation estimate. All
%  slow channels go into matrix Data_slow, and the fast channels go into
%  Data_fast. The names of the vectors go ito the cell arrays fast_list and
%  slow_list.

fast_list = [];
Data_fast = zeros(size(t_fast,1), size(index_fast,1)); % prefill
for k = 1:size(index_fast,1),
    fast_list{k} = S(index_fast(k)).name;
    Data_fast(:,k) = eval( S(index_fast(k)).name );
end

slow_list = [];
Data_slow = zeros(size(t_slow,1), size(index_slow,1));
for k = 1:size(index_slow,1)
    slow_list{k} = S(index_slow(k)).name;
    Data_slow(:,k) = eval( S(index_slow(k)).name );
end

info.fft_length    = fft_num;
info.diss_length   = round(diss_length*fs_fast);
info.overlap       = round(overlap    *fs_fast);
info.fs_fast       = fs_fast;
info.fs_slow       = fs_slow;
info.speed         = profile_speed;
info.T             = temperature_fast(n);
info.t             = t_fast(n);
info.P             = P_fast(n);
info.fit_2_isr     = fit_2_isr;
info.f_AA          = f_AA;
info.f_limit       = f_limit;
info.Data_fast     = Data_fast(n,:);
info.Data_slow     = Data_slow(m,:);
info.fast_list     = fast_list;
info.slow_list     = slow_list;

diss = get_diss_odas(SH_HP, AA(n,:), info);
% Flag piezo acceleromter
diss.piezo         = piezo;

diss.fast_list     = fast_list; % lock into returned structure
diss.slow_list     = slow_list;

if exist('spikes_A')
    diss.spikes_A      =     spikes_A; % The spikes found in piezo-accelerometer signals
    diss.pass_count_A  = pass_count_A; % The number of passes of the despike function
    diss.fraction_A    =   fraction_A; % The fraction of accelerometer data removed by the despike function
end
    
diss.spikes_sh     =     spikes_sh; % spikes found in the shear signals
diss.pass_count_sh = pass_count_sh; 
diss.fraction_sh   =   fraction_sh; 

diss.spikes_C      =     spikes_C; % The spikes found in micro_C signals
diss.pass_count_C  = pass_count_C; 
diss.fraction_C    =   fraction_C;  

diss.profiles_in_this_file = profiles_in_this_file;
diss.speed_source = speed_source;

% produce spectra of scalar variable, if any.
scalar_info.spec_length   = info.diss_length;
scalar_info.overlap       = info.overlap;


% produce spectra of scalar variable, if any.
count              = 0; % initialize
scalar_vector_list = [];
scalar_vectors     = [];
diff_gain          = [];

obj = setupstr(setupfilestr); % get differentiator gains
scalar_wish_list = {'gradT', 'gradT1', 'gradT2', 'gradC', 'gradC1', 'gradC2'};
name_wish_list   = {'T_dT' , 'T1_dT1', 'T2_dT2', 'C_dC' , 'C1_dC1', 'C2_dC2'};
we_have_scalars = false;

for index = 1:length(scalar_wish_list)
    if exist(scalar_wish_list{index},'var')
        count = count + 1;
        scalar_vector_list{count} = scalar_wish_list{index};
        eval(['scalar_vectors = ' ...
            '[ scalar_vectors ' scalar_wish_list{index} '(n)];'])
        tmp = setupstr(obj, name_wish_list{index}, 'diff_gain');  % Find diff_gain
        diff_gain = [diff_gain str2double(tmp{1})];% Assume it is numeric

    end
end
num_of_scalar_spectra = count;
if count >0, we_have_scalars = true; end
% the number of scalar vectors for which we have a spectrum

scalar_info.fft_length    = fft_num;
scalar_info.spec_length   = info.diss_length;
scalar_info.overlap       = info.overlap;
scalar_info.fs            = fs_fast;
scalar_info.gradient_method = 'first_difference';
scalar_info.diff_gain     = diff_gain;
scalar_info.f_AA          = 98;

scalar_spectra = []; % in case there are no scalars for spectra


if we_have_scalars
    scalar_spectra = get_scalar_spectra_odas(scalar_vectors, info.P, ...
        info.t, info.speed, scalar_info);
end

diss.scalar_spectra     = scalar_spectra; % lock into returned structure
diss.scalar_vector_list = scalar_vector_list;
diss.scalar_info        = scalar_info;
% Set the return value.
diss.ql_info    = ql_info;
diss.ql_info_in = ql_info_in;
result = diss;

%
% Now extract dissipation rates and spectra for plotting
% Skip this section if you do not have at least 1 shear probes

if shear_count < 1 ||  ~make_figures, return; end
fig_num = fig_num + 1;
figure(fig_num); clf

e1 = diss.e(1,:)';
if shear_count > 1
    e2 = diss.e(2,:)';
    ratio_limit = 2.5;
    e_ratio = e1 ./ e2; % dissipation ratio
    bad_ratio = find(e_ratio > ratio_limit);
    if ~isempty(bad_ratio),e1(bad_ratio) = nan(size(bad_ratio));end% e1 is too large
    bad_ratio = find(e_ratio < 1/ratio_limit);
    if ~isempty(bad_ratio),e2(bad_ratio) = nan(size(bad_ratio));end% e2 is too large
end
P = diss.P;
method = diss.method;
%assuming sh1 exists
% [junk, J1_0] = find(method(1,:) == 0);
% [junk, J1_1] = find(method(1,:) == 1);
% if size(method,1) > 1
%     [junk, J2_0] = find(method(2,:) == 0);
%     [junk, J2_1] = find(method(2,:) == 1);
% end

    h=semilogx(diss.e', P, '-o','linewidth', 2); hold on
%     semilogx(e(1, J1_0)', P(J1_0), 'ob', 'markersize', 10, 'markerfacecolor', 'b');
%     semilogx(e(1, J1_1)', P(J1_1), '^b', 'markersize', 10, 'markerfacecolor', 'b');
%     semilogx(e(2, J2_0)', P(J2_0), 'or', 'markersize', 10, 'markerfacecolor', 'r');
%     semilogx(e(2, J2_1)', P(J2_1), '^r', 'markersize', 10, 'markerfacecolor', 'r');
hold off
grid on
set(h(1),'color','b', 'markersize', 10, 'markerfacecolor', 'b')
if shear_count>1,
    set(h(2),'color','r', 'markersize', 10, 'markerfacecolor', 'r')
end

set(gca,'ylim',P_limits)

x_limits = get(gca,'xlim');
if log10(x_limits(2) / x_limits(1)) < 1
    x_limits(1) = 10^(floor(log10(x_limits(1))));
    x_limits(2) = 10^( ceil(log10(x_limits(2))));
    set(gca,'xlim',x_limits)
end

set(gca,'ydir','rev')
ylabel ('\it P \rm [dBar]')
xlabel ('\epsilon [W kg^{-1}]')
%legend_string = {'\epsilon_1','\epsilon_2', ...
%    'Method=0','Method=1', 'Method=0','Method=1'};
legend_string = [];
if shear_count >= 1
    legend_string{1} = '\epsilon_1';
end
if shear_count >=2
    legend_string{2} = '\epsilon_2';
end

legend(legend_string,'location','NorthEastOutside');

title(title_string)
drawnow

return
%
% Now a spectrogram of shear
% I am assuming that sh1 and sh2 exist and do not handle the case of more
% shear probes
fig_num = fig_num + 1;
figure(fig_num); clf

sh1_spec = squeeze(diss.sh_clean(:,1,1,:));
if shear_count>1, sh2_spec = squeeze(diss.sh_clean(:,2,2,:)); end

subplot(1,2,1)
pcolor(diss.K', diss.P, log10(sh1_spec'));grid on
%set(gca,'XScale','log')
set(gca,'xlim',[1 150])
%set(gca,'zlim',[-6 -1])
caxis([-8 -1])

set(gca,'ydir','rev')
set(gca,'xdir','rev')
set(gca,'tickdir','out')
colorbar('location','eastoutside')
shading interp
ylabel('\it P \rm [dBar]')
xlabel('\it k \rm [cpm]')

hold on
h=plot(diss.K_max(1,:), diss.P);
set(h(1),'linewidth', 2, 'color','w')
if any(diss.method(1,:) == 1)
    index = find(diss.method(1,:) == 1);
    h=plot(diss.K_max(1,index), diss.P(index), '*');
    set(h(1),'linewidth', 2, 'color','w')
end
hold off
set(gca,'layer','top')

new_title = title_string;
new_title{end} = [new_title{end} '    \Phi(\itk\rm) ' shear_list{1}];
title(new_title)

set(gca,'ylim',P_limits)

if shear_count > 1
    subplot(1,2,2)
    pcolor(diss.K', diss.P, log10(sh2_spec'));grid on
    %set(gca,'XScale','log')
    set(gca,'xlim',[1 150])
    %set(gca,'zlim',[-6 -1])
    caxis([-8 -1])

    set(gca,'ydir','rev')

    colorbar('location','eastoutside')
    shading interp
    xlabel('\it k \rm [cpm]')

    hold on
    h=plot(diss.K_max(2,:), diss.P);
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(2,:) == 1)
        index = find(diss.method(2,:) == 1);
        h=plot(diss.K_max(2,index), diss.P(index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end
    hold off
    set(gca,'tickdir','out')
    set(gca,'layer','top')
    set(gca,'ylim',P_limits)
    title(['\Phi(\itk\rm) ' shear_list{2}])
end

% Now a spectragram of shear using frequency
fig_num = fig_num + 1;
figure(fig_num); clf

sh1_spec = sh1_spec ./ repmat(diss.speed', size(sh1_spec,1), 1); % scale to preserve variance
if shear_count>1
    sh2_spec = sh2_spec ./ repmat(diss.speed', size(sh1_spec,1), 1);
end

x_lim = [1 150]; if tidal_channel, x_lim = [1 250]; end

subplot(1,2,1)
pcolor(diss.F', diss.P, log10(sh1_spec'));grid on
%set(gca,'XScale','log')
set(gca,'xlim',x_lim)
%set(gca,'zlim',[-6 -1])
caxis([-8 -1])

set(gca,'ydir','rev')
set(gca,'xdir','rev')
set(gca,'tickdir','out')
colorbar('location','eastoutside')
shading interp
ylabel('\it P \rm [dBar]')
xlabel('\it f \rm [Hz]')

hold on
h=plot(diss.K_max(1,:)' .* diss.speed, diss.P);
set(h(1),'linewidth', 2, 'color','w')
if any(diss.method(1,:) == 1)
    index = find(diss.method(1,:) == 1);
    h=plot(diss.K_max(1,index)' .* diss.speed(index), diss.P(index), '*');
    set(h(1),'linewidth', 2, 'color','w')
end
hold off
set(gca,'layer','top')
set(gca,'ylim',P_limits)
new_title = title_string;
new_title{end} = [new_title{end} '    \Phi(\itf\rm) ' shear_list{1}];
title(new_title)

if shear_count >1
    subplot(1,2,2)
    pcolor(diss.F', diss.P, log10(sh2_spec'));grid on
    %set(gca,'XScale','log')
    set(gca,'xlim',x_lim)
    %set(gca,'zlim',[-6 -1])
    caxis([-8 -1])

    set(gca,'ydir','rev')
    set(gca,'tickdir','out')

    colorbar('location','eastoutside')
    shading interp
    xlabel('\it f \rm [Hz]')

    hold on
    h=plot(diss.K_max(2,:)' .* diss.speed, diss.P);
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(2,:) == 1)
        index = find(diss.method(2,:) == 1);
        h=plot(diss.K_max(2,index)' .* diss.speed(index), diss.P(index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end
    hold off
    set(gca,'layer','top')
    set(gca,'ylim',P_limits)
    title(['\Phi(\itk\rm) ' shear_list{2}])
end

% Now a spectragram of accelleration
fig_num = fig_num + 1;
figure(fig_num); clf

c_lim = zeros(  accel_count,2); % used to set color scale
h     = zeros(1,accel_count);

if accel_count <1, return, end
if accel_count >0, Ax_spec = squeeze(diss.AA(:,1,1,:));end
if accel_count >1, Ay_spec = squeeze(diss.AA(:,2,2,:));end
if accel_count >2, Az_spec = squeeze(diss.AA(:,3,3,:));end

h(1)=subplot(1,accel_count,1);
pcolor(diss.F', diss.P, log10(Ax_spec'));grid on
set(gca,'xlim',x_lim)
c_lim(1,:) = get(gca,'clim');

set(gca,'ydir','rev')
set(gca,'tickdir','out')
colorbar('location','eastoutside')
shading interp
ylabel('\it P \rm [dBar]')
xlabel('\it f \rm [Hz]')

set(gca,'layer','top')
set(gca,'ylim',P_limits)
new_title = title_string;
new_title{end} = [new_title{end} '    \Phi(\itf\rm) ' accel_list{1}];
title(new_title)

if accel_count > 1
    h(2)=subplot(1,accel_count,2);
    pcolor(diss.F', diss.P, log10(Ay_spec'));grid on
    %set(gca,'XScale','log')
    set(gca,'xlim',x_lim)
    c_lim(2,:) = get(gca,'clim');

%    caxis([0 6])

    set(gca,'ydir','rev')
    set(gca,'tickdir','out')

    colorbar('location','eastoutside')
    shading interp
    xlabel('\it f \rm [Hz]')

    set(gca,'layer','top')
    set(gca,'ylim',P_limits)
    title(['\Phi(\itf\rm) ' accel_list{2}])
end
if accel_count > 2
    h(3)=subplot(1,accel_count,3);
    pcolor(diss.F', diss.P, log10(Az_spec'));grid on
    %set(gca,'XScale','log')
    set(gca,'xlim',x_lim)
    c_lim(3,:) = get(gca,'clim');

%    caxis([0 6])

    set(gca,'ydir','rev')

    shading interp
    xlabel('\it f \rm [Hz]')
    colorbar('location','eastoutside')

    set(gca,'layer','top')
    set(gca,'ylim',P_limits)
    set(gca,'tickdir','out')
    title(['\Phi(\itf\rm) ' accel_list{3}])
end
c_lim = [min(c_lim(:,1)) , max(c_lim(:,2))];

for index = 1:accel_count
    set(h(index),'clim', c_lim)
end

end
% This is the end


function printfile(fig_num, name, profile, eps, render)
    if nargin < 5
        render = '-painters';
    end
    if eps
        if ~exist('export_fig','file')
            disp('You are missing the function "export_fig".');
            disp('It can be found in the Matlab file exchange.');
            disp('Please download this function and add to your path');
            disp('to enable exporting figures');
            disp('  ');
            return
        end
        handle = figure(fig_num);
        filename = sprintf('%s_P_%02d_Fig_%02d.pdf',name,profile,fig_num);
        export_fig(handle, filename, render, '-transparent');
        %saveas(fig_num, print_file)
        %fig2pdf(gcf,    print_file, [], '')
    end
end



% Assemble input arguments into a 1X1 stracture.  Contents of the input
% struture are joined with parameter name / value pairs.  The resulting
% structure can then be saved for future reference.
function args = join_arguments(input)
    if isempty(input)
        args = [];
        return
    end

    skip = 0;
    for k = 1:length(input)
        if skip, skip = 0; continue; end
    
        if isstruct(input{k})
            args = input{k};
            continue
        else
            if ~ischar(input{k})
                error('Invalid input argument');
            end
            args.(input{k}) = input{k+1};
            skip = 1;
        end
    end
end
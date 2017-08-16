%% quick_look_HMP - EXCLUDED
% Old, outdated version of quick_look for HMP instruments.
%%
% <latex>\index{Functions!quick\_look\_HMP}</latex>
%
%%% Syntax
% @defgroup typeB Type B functions
%
% function result = quick_look_HMP (fname, t_start, t_end, ql_info)
%
% ===========================================================================
% @ingroup typeB
% @file quick_look_HMP.m
% @brief Quickly evaluate data file recorded by your instrument
% ===========================================================================
%
% This is a function for the rapid generation of plots and figures
% from profiles recorded with a horizontal profiler.
% It is intended for the quick evaluation of an instrument during a sea trail.
% This function only works with Version 6 and later data files that have
% the setup-configuration file embedded within the data file, so that all
% variable can be converted into physical units.
% The first 3 input arguments are required. gl_info is a structure that
% may contain information to taylor the processing done by this function.
%
% @param fname name of the file to be processed (extension optional) [string]
% @param t_start is the starting time of the data section used for
%       spectral estimates [integer].
% @param t_end is the end of this segment. Only the spectral plots use this
%       value [integer]. t_end > t_start, only.
%
% A structure containing the default parameters expected within ql_info is
% returned when quick_look is called with no input parameters.  This
% structure can be examined, edited, and used as the ql_info input
% parameter.
%
% The optional parameters in the ql_info structure are:
%  @param min_speed - the minimum speed, in dBar per second, of a profile.
%       Default = 0.2.
%  @param profile_num - an integer identifying the profile number to be
%       analyzed. Default = 1. Some files may contain multiple profiles.
%  @param fft_length - length, in seconds, of the fft-segments used for
%       spectral estimation. Default = 2.
%  @param fit_order - order of the polynomial used to fit to the shear
%       spectra in order to determine the minimum in the shear spectrum.
%       This minimum is used to set the upper limit of spectral integration
%       to derive the variance of shear. A good choice is 6, but this
%       depends on the nature of your instrument. Default = 6.
%  @param diss_length - length, in seconds, of dissiaption estimates. It
%       must be >= 2*fft_length. Default = 4*default_fft_length.
%  @param overlap - the overlap, in seconds, used for dissipation estimates.
%       It should be diss_length/2. Default = default_diss_length / 2.
%  @param HP_cut - frequency, in Hz, of high-pass filter applied to the shear
%       data. Use a value equivalent to ~0.5 cpm. Default = 0.4.
%  @param LP_cut - frequency, in Hz, of low-pass filter applied to the
%       shear, thermistor, and micro-conductivity data. It is used only for
%       the plotting of depth profiles. It is not used for spectral
%       estimation. Default = 30.
%  @param eps - 0 if you do not want to produce graphics files of the figures.
%       Set it to 1, if you do want to save the graphics files. The
%       file-formats are 'fig', and 'pdf'. By uncommenting lines, you can
%       also generated 'png' and 'eps' files. Default = 0.
%  @param fit_2_isr - The dissipation rate above which the dissipation
%       algorithm will make its estimate by fitting to the inertial subrange
%       rather than by spectral integration. Default fit_2_isr = 1.5e-5 W/kg.
% @param despike_sh_1 a vector of three values [A B C] for despiking the
%       shear probe signals on the first pass. A is the threshold B is the
%       low-pass smoothing frequency and C is the length of data to be
%       replaced around a spike in units of seconds. Default value is
%       [7 0.5 0.07]. Set threshold to A=inf, to not despike the shear
%       data.
% @param despike_A_1 a vector of three values [A B C] for despiking the
%       piezo-accelerometer signals on the first pass. Default value is
%       [7 0.5 0.07]. Set threshold to A=inf, to not despike the
%       piezo-accelerometers.
% @param despike_C a vector of three values [A B C] for despiking the
%       micro-conductivity signals on the first pass. Default value is
%       [10 0.5 0.08].
% @parme f_limit. The absolute maximum frequency to use when estimating the
%       rate of dissipation of TKE with the shear probe signals. Use this
%       parameter if you have vibrational contamination that is not removed
%       by the Goddman coherent noise removal function
%       (clean_shear_spectrum). Default value is f_limit=inf.
% @param YD_0 the year day of the start of data collection. The time axis
%       of all plots will be moved by subtracting this day from the time
%       coordinate. This reduces the number of digits on the time axis and makes
%       it reading it easier. Default = 0.
% @param info.fit_2_isr. The value of dissipation rate at which the estimation
%       will be based on a fit to the inertial subrange rather than by spectral
%       integration. Default fit_2_isr = 1e-5 W/kg.
% @param constant_speed. Force the data conversion to use a constant speed
%       for converting shear probe data into physical units, and for creating
%       gradients. The velocity components provided by the Nortek CM might be
%       very noisy and not directly suitable for data conversion. Cannot be a
%       negative value and must be larger than 0.1 m/s.
%
% Output:
%   @param result - a structure containing spectral information and dissipation
%       rates, that were returned by the function get_diss. See get_diss for
%       details. This structure may be all that you need to see all spectra in
%       a profile of all shear probes, thermistors, conductivity probes and
%       accelerometers. It also contsins the rate of
%       dissipation, wavenumbers and the limits of integration for each
%       estimate of the rate of dissipation. Suitable for paper-ready plots.
%
%       Optionally, if quick_look is called with no arguments, result will
%       contain a copy of the default ql_info input structure.
%
% <b>Revision history:</b>
% - 2006-04-15 2012-06-26 (RGL) A great many changes over these years.
% - 2012-07-13 (RGL) major changes to input arguments and how they are
%      handled.
% - 2013-01-12 (RGL) modified quick_look_float to make it useful for all
%       types of vertical profilers, except XMP. Now called quick_look_new.
%       This name will be change to quick_look once I am happy with the
%       results.
% - 2013-04-09 (RGL) added spectogram, chnaged the wave figures are saved,
%       modified tolerence of gediss so that bandwidth of integration is
%       not overly reduced by low shear spectrum below 5 cpm.
%       Renamed this function to quick_look. The original quick_look
%       function will be named quick_look_old
% - 2013-04-24 (WID) return copy of default argument structure when no
%       inputs are provided.
% - 2013-05-07 (WID) input parser set to not result in an error when an
%       input parameter does not match.
% - 2013-07-16 (RGL) added spectrogram of Ax and Ay
% - 2013-07-27 (RGL) added rejection of epsilon if their ratio
%       exceeds 2.5 for the dissipation profile figure
% - 2013-11-08 (RGL) added spike removal from magnetometer data
% - 2014-03-13 (RGL) update to act more like quick_look. Added scalar
%       spectra to diss-structure. Added an upper limit to frequency to avoid the
%       bad vibrations above 200 Hz, at high (U>2 m/s) speeds.
% - 2014-03-19 (RGL) got spectragrams to work. Added feature to force using
%       a constant speed for conversion of data into physical units in case the
%       Vector CM data is very noisy.
% - 2014-04-02 (RGL) final tweeking of all sorts of little things.
% - 2014-04-03 (RGL) Corrected year-day value so that it is relative to the
%       year of the file rather than to 2013.
% - 2015-07-27 (WID) Use texstr in place of fix_underscore
%
% ==============================================
function result = quick_look_HMP( fname, t_start, t_end, varargin )

%
% Default values for optional fields
default_min_speed     = 0.2;
default_fft_length    = 2;% in units of seconds
default_diss_length   = 60;
default_overlap       = 0; % 0 percent over_lap
default_profile_num   = 1; % process the first profile
default_HP_cut        = 0.5; % high-pass shear probes at 0.5 Hz
default_LP_cut        = 100; % low-pass profile at 30 Hz, FOR PROFILE DISPLAY PURPOSES ONLY!
default_eps           = 0; % do not render any graphics files
default_fit_2_isr     = 1.5e-5; % W/kg
default_fit_order     = 3;
default_f_AA          = 392; % anti-aliasing filter
default_f_limit       = inf;
default_YD_0          = 0;
default_constant_speed= [];
% The despiking parameters are [thresh, smooth, and length (in seconds)
default_despike_sh_1  = [ 7  0.5 0.02]; % for first pass shear probes
default_despike_A_1   = [ 8  0.5 0.02]; % for first pass piezo-accelerometers
default_despike_C     = [10  0.5 0.08]; % for first and only pass of micro-C

if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    return
end

p = inputParser;
p.KeepUnmatched = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_speed       = @(x) isnumeric(x) && isscalar(x)   || isempty(x);
val_vector      = @(x) (isnumeric(x) && isvector(x)) || isempty(x);
val_string      = @(x) ischar(x);

addRequired(  p, 'fname',        val_string);
addRequired(  p, 't_start',      val_numeric);
addRequired(  p, 't_end',        val_numeric);
addParamValue(p, 'fft_length',   default_fft_length,     val_numeric);
addParamValue(p, 'diss_length',  default_diss_length,    val_numeric);
addParamValue(p, 'overlap',      default_overlap,        val_numeric);
addParamValue(p, 'fit_order',    default_fit_order,      val_numeric);
addParamValue(p, 'profile_num',  default_profile_num,    val_numeric);
addParamValue(p, 'min_speed',    default_min_speed,      val_numeric);
addParamValue(p, 'HP_cut',       default_HP_cut,         val_numeric);
addParamValue(p, 'LP_cut',       default_LP_cut,         val_numeric);
addParamValue(p, 'eps',          default_eps,            val_numeric);
addParamValue(p, 'f_AA',         default_f_AA,           val_numeric);
addParamValue(p, 'fit_2_isr',    default_fit_2_isr,      val_numeric);
addParamValue(p, 'f_limit',      default_f_limit,        val_numeric);
addParamValue(p, 'despike_sh_1', default_despike_sh_1,   val_vector);
addParamValue(p, 'despike_A_1',  default_despike_A_1,    val_vector);
addParamValue(p, 'despike_C'  ,  default_despike_C,      val_vector);
addParamValue(p, 'YD_0',         default_YD_0,           val_numeric);
addParamValue(p, 'constant_speed', default_constant_speed, val_speed);

% Parse the arguments.
parse(p, fname, t_start, t_end, varargin{:});

% Perform last stages of input validation.
if p.Results.diss_length < 2*p.Results.fft_length,
  error('Invalid size for diss_length - must be greater than 2 * fft_length.');
end
if p.Results.t_end < p.Results.t_start,
  error('Starting pressure must be less than end pressure.');
end

fft_length    = p.Results.fft_length;
diss_length   = p.Results.diss_length;
overlap       = p.Results.overlap;
min_speed     = p.Results.min_speed;
profile_num   = p.Results.profile_num;
HP_cut        = p.Results.HP_cut;
LP_cut        = p.Results.LP_cut;
eps           = p.Results.eps;
fit_order     = p.Results.fit_order;
f_AA          = p.Results.f_AA;
f_limit       = p.Results.f_limit;
fit_2_isr     = p.Results.fit_2_isr;
despike_sh_1  = p.Results.despike_sh_1;
despike_A_1   = p.Results.despike_A_1;
despike_C     = p.Results.despike_C;
YD_0          = p.Results.YD_0;
constant_speed = p.Results.constant_speed;

YD_0 = floor(YD_0); % Force it to be w whole number

% Now save these parameter values so that they can be placed into the
% dissipation structure that is returned by this function, quick_look
ql_info = p.Results;

% end of input argument checking.
%
% ____________________________________________________________
% File opening etc.
% Check file name, check for a .mat file, and open the .p file if no *.mat
% file
% _________________________________________________________________

we_have_old_mat_file = 0;% if mat-file already exists, assume it is in physical units.

% Look for the .p and .mat files, deal with problems
error_str = ['File "' fname '" not found.  Only data files (.p) are valid.'];
[fpath,N,E,file] = file_with_ext( fname, {'.p', '.P', ''}, error_str );

if exist([fpath filesep N '.mat'],'file') % Then the mat-file already exist and we will
    % use it and go straight to the plotting.
    disp(['Loading file  = ' fpath filesep N '.mat'])
    load([fpath filesep N '.mat'])
    we_have_old_mat_file = 1;
else
    % check for bad buffers
    fix_manually = [];
    bad_records = check_bad_buffers(file);
    if ~isempty(bad_records)
        disp(['The file "' file '" has bad buffers'])
        disp('Making backup copy')
        copyfile(file, [fpath filesep N '_original' E]);
        disp('Calling function patch_odas')
        [bad_records, fix_manually] = patch_odas(file);
        if exist(fix_manually, 'var') && ~isempty(fix_manually)
            disp('Some records cannot be patched')
            disp('Calling function fix_bad_buffers')
            fix_bad_buffers(file);
        end
    end

    variable_list = read_odas(file); % convert to a mat-file
    if isempty(variable_list)
        error('Could not convert data into a mat-file');
    end
    disp(['Loading file  = ' fpath filesep N '.mat'])
    load([fpath filesep N '.mat']); % use the one we just converted
end

%
%__________________________________________________________________
%
% This is where we get some information about the data sampling

if (~we_have_old_mat_file), %We do not have a file that is already in
    % physical units. So, we will create this mat-file.

    number_of_records      = size(header,1); % number of records in this file.
    length_of_record       = (header(1,19) - header(1,18))/2; % number of samples in a record
    number_of_fast_columns = header(1,29);
    number_of_slow_columns = header(1,30);
    number_of_columns      = number_of_fast_columns + number_of_slow_columns;
    number_of_rows         = header(1,31);
    length_of_fast_vector  = (length_of_record / number_of_columns)*number_of_records; % for standard practice
    length_of_slow_vector  = length_of_fast_vector / number_of_rows; % for standard practice
    % double sampled channels and ground "fill" channels may have different
    % lengths

    t_slow = (0:length_of_slow_vector - 1)'/fs_slow; % make time vectors for plotting
    t_fast = (0:length_of_fast_vector - 1)'/fs_fast;

    Year  = header(1,4);
    Month = header(1,5);
    Day   = header(1,6);
    Hour  = header(1,7);
    Minute= header(1,8);
    Second= header(1,9);

    % Deconvolve the channels with pre-emphasis
    deconvolve_list = {'T1_dT1', 'T2_dT2', 'C1_dC1', 'C2_dC2'};
    new_name_list   = {'T1_hres', 'T2_hres', 'C1_hres', 'C2_hres'};
    for k = 1:length(deconvolve_list)
        if exist(deconvolve_list{k},'var')
            eval([new_name_list{k} '=deconvolve(''' deconvolve_list{k} ''', [],' ...
                deconvolve_list{k} ', fs_fast, setupfilestr, header_version);'])
        end
    end

    % Pressure is special
    P_hres  = deconvolve('P_dP',   P,  P_dP,   fs_slow, setupfilestr, header_version);
    % end of deconvolution
    %
    % Now convert to physical units. The following is a list of possible
    % signals. The last four are for a MicroSquid instrument.
    convert_list = {...
        'T1_hres', 'T2_hres', 'C1_hres', 'C2_hres', 'P', 'P_hres', ...
        'PTV', 'sbt', 'sbc', 'SBT1', 'SBT2', 'SBT3', 'SBC1', 'SBC2', 'SBC3', ...
        'Ax', 'Ay', 'Az', 'Rx', 'Ry', 'Rz', 'Mx', 'My', 'Mz', ...
        'APy', 'APz', ...
        'U', 'V', 'W', 'V_Bat', 'Incl_X', 'Incl_Y', 'Incl_T',...
        'sh1', 'sh2', 'sh3', 'sh4', ...
        'T_dT', 'C_dC', 'C_dC', 'MS_DO', 'Fluo', 'BS', 'Suna', 'Pitch', ...
        'JAC_C', 'JAC_T', 'Turbidity', 'Chlorophyll'};
    name_list = convert_list;
    % use coefficients for non-differentiated channel
    if exist('T1','var'), name_list{1} = 'T1'; else name_list{1} = 'T1_dT1';end;
    if exist('T2','var'), name_list{2} = 'T2'; else name_list{2} = 'T2_dT2';end
    if exist('P' ,'var'), name_list{6} = 'P' ; else name_list{6} = 'P_dP'  ;end
    name_list{3} = 'C1_dC1';
    name_list{4} = 'C2_dC2';

    if exist('SBT1','var'), SBT1 = SBT1 + 1e-12;end
    if exist('SBC1','var'), SBC1 = SBC1 + 1e-12;end

    % name_list gives the names of the channels that contain the conversion
    % coefficients. For example, "P" is a variable and is found in the
    % configuration file, but P_hres is not in the coefficient file. The
    % coefficients for P_hres are in the section with the name "P".
    % Similarly for "T1_hres", etc. Most of the variables in the convert_list
    % have names that are identical to those in the configuration file.

    % Note that sh1 -- sh4 (shear probes) still must be divided by speed
    % squared in order to complete the conversion to physical units.

    for k = 1:length(convert_list)
        if exist(convert_list{k},'var')
            eval([convert_list{k} '=' ...
                'convert_odas(' convert_list{k} ', '''  name_list{k} ''', ' ...
                '''string''' ' , setupfilestr, header_version);']);
        end
    end

    % Generate pressure for fast sampling rate
    P_fast = interp1(t_slow, P_hres,t_fast,'spline','extrap');
    % Generate year-day time vectors
    t_fast_YD = ...
        datenum(Year, Month, Day, Hour, Minute, Second) - ...
        datenum(Year,1,1,0,0,0) +1  + (t_fast - 1)/(3600*24);
    t_slow_YD = ...
        datenum(Year, Month, Day, Hour, Minute, Second) - ...
        datenum(Year,1,1,0,0,0) +1 + (t_slow - 1)/(3600*24);

    % We now generate the velocity and speed vectors for this data file.
    % There are three posibilities:
    % (1) We have recorded in this file the vectors U, V, and W from a
    %   Nortek vector current meter.
    % (2) We have the Vector data but it is located in a separate file
    %   (such as for the Grand Passage data set of 2013).
    % (3) We have no fricking velocity data whatsoever, in which case we
    %   simply take speed=1, U=1, V=0, and W=0, and figure out the actual
    %   speed separately.
    % For case (1) we might have an issue with noise and so it might
    %   be prudent to initially set the velocity to 1 m/s until an
    %   appropriate smoothing is developed. You can force the speed used
    %   for converting the shear and gradient channels to be a constant value
    %   the optional constant_speed parameter in the variable input arguments.

    if exist('U','var') && exist('V','var') && exist('W','var') % they are slow channels
        U_MR_slow  = U; % We need a separate name
        V_MR_slow  = V;
        W_MR_slow  = W;
        speed_slow = sqrt(U.^2 + V.^2 + W.^2);
        U_MR_fast  = interp1(t_slow, U, t_fast, 'pchip','extrap');
        V_MR_fast  = interp1(t_slow, V, t_fast, 'pchip','extrap');
        W_MR_fast  = interp1(t_slow, W, t_fast, 'pchip','extrap');
        speed_fast = interp1(t_slow, speed_slow, t_fast, 'pchip','extrap');
        speed_source = 'Directly recorded U V W';
    elseif exist('vector.mat','file') % we have the data in a file
        load ('vector.mat')
        % Generate time base of MR data in Year-Day units with January 1 equal
        % to 1.

        U_MR_slow = interp1(t_vector, U, t_slow_YD, 'pchip','extrap');
        V_MR_slow = interp1(t_vector, V, t_slow_YD, 'pchip','extrap');
        W_MR_slow = interp1(t_vector, W, t_slow_YD, 'pchip','extrap');

        U_MR_fast = interp1(t_vector, U, t_fast_YD, 'pchip','extrap');
        V_MR_fast = interp1(t_vector, V, t_fast_YD, 'pchip','extrap');
        W_MR_fast = interp1(t_vector, W, t_fast_YD, 'pchip','extrap');

        speed_fast = sqrt(U_MR_fast.^2 + V_MR_fast.^2 + W_MR_fast.^2); %
        speed_slow = sqrt(U_MR_slow.^2 + V_MR_slow.^2 + W_MR_slow.^2); %
        speed_source = 'From file = vector.mat';
    else % We have no vector data
        U =  ones(size(t_slow));
        V = zeros(size(t_slow));
        W = V;
        speed_slow = U;
        speed_fast =  ones(size(t_fast));
        U_MR_slow  = U;
        V_MR_slow  = V;
        W_MR_slow  = V;
        U_MR_fast  =  ones(size(t_fast));
        V_MR_fast  = zeros(size(t_fast));
        W_MR_fast  = zeros(size(t_fast));
        speed_source = 'No source, using 1 m/s';
    end
    if ~isempty(constant_speed) % force a constant speed
        speed_fast =  constant_speed * ones(size(t_fast));
        speed_source = 'Forced to constant speed by input parameter';
    end

    junk = find(speed_fast < min_speed);
    if ~isempty(junk)
        speed_fast(junk) = min_speed;
        % to avoid divide by zero or very small number.
    end
    junk = find(speed_slow < min_speed);
    if ~isempty(junk)
        speed_slow(junk) = min_speed;
    end

    if exist('sh1','var'), sh1 = sh1 ./ speed_fast.^2; end
    if exist('sh2','var'), sh2 = sh2 ./ speed_fast.^2; end
    if exist('sh3','var'), sh3 = sh3 ./ speed_fast.^2; end
    if exist('sh4','var'), sh4 = sh4 ./ speed_fast.^2; end

    % generate gradient signals
    gradient_list        = {'gradT1',  'gradT2',  'gradC1',  'gradC2'};
    gradient_source_list = {'T1_hres', 'T2_hres', 'C1_hres', 'C2_hres'};

    for k = 1:length(gradient_list)
        if exist(gradient_source_list{k},'var')
            eval([gradient_list{k} ' = ' gradient_source_list{k} ';']);
            eval([gradient_list{k} '(2:end) = diff(' gradient_source_list{k} ') *fs_fast;']);
            eval([gradient_list{k} '(1) = ' gradient_list{k} '(2);']);
            eval([gradient_list{k} ' = ' gradient_list{k} './ speed_fast;']);
        end
    end

    obj = setupstr(setupfilestr);
    if exist('Fluo','var')
       offset = str2double(setupstr(obj,'Fluo', 'adc_fs'));
       Fluo = Fluo + offset/2;
    end
    if exist('BS','var')
       offset = str2double(setupstr(obj,'BS', 'adc_fs'));
       BS = BS + offset/2;
    end
    if exist('Suna','var')
       offset = str2double(setupstr(obj,'Suna', 'adc_fs'));
       Suna = Suna + offset/2;
    end

% cleanup the fliers in the magnetometer data
if exist('Mx','var') || exist('My','var') || exist('Mz','var')
    num_points_per_second = round(fs_slow);
    select = 1:round(fs_slow); %
    num_segments = length(Mx)/num_points_per_second;
    threshold = 10; % in uT
    for index = 1:num_segments
        median_x = median(Mx(select));
        median_y = median(My(select));
        median_z = median(Mz(select));
        nx = find(abs(Mx(select) - median_x) >= threshold);
        ny = find(abs(My(select) - median_y) >= threshold);
        nz = find(abs(Mz(select) - median_z) >= threshold);
        if ~isempty(nx), Mx(select(nx)) = median_x; end
        if ~isempty(ny), My(select(ny)) = median_y; end
        if ~isempty(nz), Mz(select(nz)) = median_z; end
        select = select + num_points_per_second;
    end
end

    % Save the converted variables to the mat-file
    eval (['save ' fpath filesep N ...
        ' t_slow t_fast Year Month Day Hour Minute speed_source -append -v6'])
    if exist('accel_shear_info','var')
        eval (['save ' fpath filesep N ' accel_shear_info -append -v6'])
    end

    var_list = {...
        'P', 'P_hres', 'T1_hres', 'T2_hres', 'sbt', 'sbc', ...
        'SBT1', 'SBT2', 'SBT3', 'SBC1', 'SBC2', 'SBC3', ...
        'Ax', 'Ay', 'Az', 'Rx', 'Ry', 'Rz', 'GC_T', 'APy', 'APz', ...
        'W_fast', 'W_slow', 'P_fast', 'sh1', 'sh2', 'sh3', 'sh4', ...
        'gradT1', 'gradT2', ...
        'gradC1', 'C1_hres', 'gradC2', 'C2_hres', 'Mx', 'My', 'Mz', ...
        'speed_fast', 'speed_slow', ...
        'U_MR_slow', 'V_MR_slow', 'W_MR_slow', ...
        'U_MR_fast', 'V_MR_fast', 'W_MR_fast', ...
        't_fast_YD', 't_slow_YD', ...
        'U', 'V', 'W', 'V_Bat', 'V_bat', 'Incl_X', 'Incl_Y', 'Incl_T', 'PTV', ...
        'Fluo', 'BS', 'Suna', 'Pitch', 'JAC_C', 'JAC_T', 'Turbidity', 'Chlorophyl'};

    for k = 1: length(var_list)
        if exist(var_list{k},'var')
            flag = save_odas([fpath filesep N '.mat'], var_list{k},  eval(var_list{k}) );
            if (flag < 0), error(['Could not write ' var_list{k} ' to ' fpath filesep N '.mat']), end
        end
    end


    % End of conversion to physical units' for version 6
end

% This is the end of the section for converting to physical units. This
% section was skipped if the mat-file already existed. If the mat-file is
% corrupted or has a bad conversion, then erase it and run this function
% again

%
% ____________________________________

Year_string        = num2str(Year);
Month_string       = sprintf('%02d',Month);
Day_string         = sprintf('%02d',Day);
Hour_string        = sprintf('%02d',Hour);
Minute_string      = sprintf('%02d',Minute);
profile_num_string = sprintf('%02d',profile_num);

title_string = {[ 'Grand Passage, MR, ' texstr(N) '; '  Year_string '\_' Month_string '\_' Day_string ...
    ', ' Hour_string ':' Minute_string ' UTC'], ['U = ' num2str(mean(speed_slow),3) ' m s^{-1}']};

green = [0.25 0.5 0.25];
gold = 0.85*[1 1 0.25];
xlabel_string =  ['t - ' num2str(YD_0) ' [Year-Day]'];


% Make figures
% Pressure
%
% ---------------------------------------------------------------
fig_num = 0;
if exist('P_hres','var')
    fig_num = fig_num + 1;
    figure (fig_num); clf

    P_mean = mean(P_hres);
    h = plot(t_slow_YD - YD_0,  P_hres - P_mean);grid

    set(h(1), 'linewidth', 2,   'color','b')
    set(gca,'ydir','rev')
    xlabel(xlabel_string)
    ylabel ('[ dBar ]')
    title(title_string)
    legend(['P_{hres} - ' num2str(P_mean,4)], 'location','eastoutside')

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end

% Pitch, Roll
%
% ---------------------------------------------------------------
LP_ADIS = 2; % low-pass filter applied to ADIS Inclinometer
[b,a] =   butter(1,LP_ADIS/(fs_slow/2));% low-pass to show only the gravity signal
if exist('Incl_X','var') && exist('Incl_Y','var')
    fig_num = fig_num + 1;
    figure (fig_num); clf

    fs_ADIS = 482; % Internal sampling rate of the ADIS Inclinometer
    tau_N = 256; % The nuber of samples in the running average internal to the ADIS
    ADIS_delay = (tau_N/2) / fs_ADIS; % time delay of output from ADIS

    h = plot(...
        t_slow_YD - YD_0 - ADIS_delay/(3600*24), [Incl_Y Incl_X], ...
        t_slow_YD - YD_0, 180*asin([-Ax Ay]/9.81)/pi);grid

    set(h(1), 'linewidth', 2,   'color','b')
    set(h(2), 'linewidth', 2,   'color','r')
    set(h(3), 'linewidth', 1,   'color','c')
    set(h(4), 'linewidth', 1,   'color','m')
    ylabel('[ ^{\circ} ]')
    xlabel(xlabel_string)
    title(title_string)
    legend('\theta_Y','\theta_X', '-A_X','A_Y','location','eastoutside')

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end

%
% Dumb bell accelerometers before despiking
%-----------------------------------------------------------------------
% Note if DC accelerometers exist, then they are called Ax, Ay, and Az,
%   while the piezo-accelerometers are then called APy and APz. Else, the
%   piezo acceleromerers are called Ax and Ay.
if (exist('Ax','var') && exist('Ay','var') && ~exist('Az','var')) || ...
        (exist('APy','var') && exist('APz','var'))
    % Then we have piezo-accelerometers, which have no DC response
    if exist('APy','var') && exist('APz','var')
        plot_data = [APy APz];
        legend_string{1} = 'AP_Y';    legend_string{2} = 'AP_Z';
        legend_string{3} = 'AP_Y LP'; legend_string{4} = 'AP_Z LP';
    else
        plot_data = [Ax Ay];
        legend_string{1} = 'A_X';    legend_string{2} = 'A_Y';
        legend_string{3} = 'A_X LP'; legend_string{4} = 'A_Y LP';
    end

    fig_num = fig_num + 1;
    figure (fig_num); clf
    [b,a] =   butter(1,1/(fs_fast/2));% low-pass to show only the gravity signal

    h = plot(t_fast_YD - YD_0, [plot_data filtfilt(b,a,plot_data)]);grid
    set(h(1), 'color','b')
    set(h(2), 'color','r')
    set(h(3), 'color', 'k','linewidth',1);
    set(h(4), 'color','k','linewidth',1)
    ylabel('[counts]')
    xlabel(xlabel_string)
    legend(legend_string,'location','NorthEastOutside');
    title(title_string)
    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end

%
%-----------------------------------------------------------------------
% Dumb bell accelerometers after despiking
spikes_A1 = []; % they must exist
spikes_A2 = [];

despike_A_string = ...
    ', despike_A_1(1), despike_A_1(2), fs_fast, round(despike_A_1(3)*fs_fast));';

if (exist('Ax','var') && exist('Ay','var') && ~exist('Az','var')) || ...
        (exist('APy','var') && exist('APz','var'))
    % Then we have piezo-accelerometers, which have no DC response

        if exist('APy','var') && exist('APz','var')
            A = [APy APz];
            legend_string{1} = 'AP_Y';    legend_string{2} = 'AP_Z';
            legend_string{3} = 'AP_Y LP'; legend_string{4} = 'AP_Z LP';
        else
            A = [Ax Ay];
            legend_string{1} = 'A_X';    legend_string{2} = 'A_Y';
            legend_string{3} = 'A_X LP'; legend_string{4} = 'A_Y LP';
        end
        if despike_A_1(1) ~= inf
            [A(:,1), spikes_A1]  = eval(['despike(A(:,1)' despike_A_string]);
            [A(:,2), spikes_A2]  = eval(['despike(A(:,2)' despike_A_string]);
        end

    fig_num = fig_num + 1;
    figure (fig_num); clf

    [b,a] =   butter(1,1/(fs_fast/2));% low-pass to show only the gravity signal
    h = plot(t_fast_YD - YD_0, [A filtfilt(b, a, A)]);grid
    set(h(1), 'color','b')
    set(h(2), 'color','r')
    set(h(3), 'color', 'k','linewidth',1);
    set(h(4), 'color','k','linewidth',1)
    ylabel('[counts]')
    xlabel(xlabel_string)
    new_title_string = title_string;
    new_title_string{end+1} = ...
            ['with despiking, 1^{st} = ' num2str([length(spikes_A1) length(spikes_A2)])];
    legend(legend_string,'location','NorthEastOutside');
    title(new_title_string)
    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end

%
%-----------------------------------------------------------------------
% Vector current meter data
if exist('U_MR_slow','var') || exist('V_MR_slow','var') || exist('W_MR_slow','var')
    fig_num = fig_num + 1;
    figure (fig_num); clf
    plot_data = [];
    legend_string = [];
    if exist('U_MR_slow','var')
        plot_data = [plot_data U_MR_slow];
        legend_string{end+1} =  'U\_MR';
    end
    if exist('V_MR_slow','var')
        plot_data = [plot_data V_MR_slow];
        legend_string{end+1} =  'V\_MR';
    end
    if exist('W_MR_slow','var')
        plot_data = [plot_data W_MR_slow];
        legend_string{end+1} =  'W\_MR';
    end
    h = plot(t_slow_YD - YD_0, plot_data);grid
    set(h(1), 'color','b',   'linewidth', 2)
    set(h(2), 'color','r',   'linewidth', 2)
    set(h(3), 'color', green,'linewidth', 2);
    ylabel('[m s^{-1}]')
    xlabel(xlabel_string)
    legend(legend_string,'location','NorthEastOutside');
    title( title_string)
    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end
%
%----------------------------------------------------------
%Thermistors
if exist('T1_hres','var') || exist('T2_hres','var') || exist('T_hres','var')
    fig_num = fig_num +1;
    figure (fig_num); clf
    legend_string = [];
    plot_data = [];
    T_count = 0;
    if exist('T1','var')
        T_count = T_count + 1;
        legend_string{T_count} = '\itT\rm_1';
        plot_data = [plot_data T1_hres];
    end
    if exist('T2','var')
        T_count = T_count + 1;
        plot_data = [ plot_data T2_hres];
        legend_string{T_count} = '\itT\rm_2';
    end
    if exist('T_hres','var')
        T_count = T_count + 1;
        plot_data = [ plot_data T_hres];
        legend_string{T_count} = '\itT';
    end
    h = plot(t_fast_YD - YD_0, plot_data);grid
    legend(legend_string,  'location','NorthEastOutside')
    xlabel(xlabel_string)
    ylabel('[  ^{\circ}C ]')
    title(title_string)

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end
%
%
% __________________
% Battery Voltage
if exist('V_Bat','var')
    fig_num = fig_num + 1;
    figure(fig_num); clf
    [b,a]=butter(1,1/(fs_slow/2));
    h = plot(t_slow_YD - YD_0, [V_Bat filtfilt(b,a,V_Bat)]);grid
    %legend('\itV_{Bat}',  'location','NorthEastOutside')
    set(h(2),'linewidth',3,'color','m')
    xlabel(xlabel_string)
    ylabel('\itV_{Bat}\rm [volts]')
    title(title_string)
    legend('V_{Bat}','B_{Bat}-LP','location','eastoutside')
    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        fig2pdf(gcf,    print_file, [], '')
    end
end
%
%--------------------------------------------------------------------------
% Temperature gradient, microconductivity gradient, and shear

% First despike the existing shear probe signals

shear_wish_list = {'sh1', 'sh2', 'sh3', 'sh4'};
shear_list  = [];
shear_count = 0;
SH          = [];
plot_gradient_wish_list = [];
legend_gradient_wish_list = [];
for index = 1:length(shear_wish_list)
    if exist(shear_wish_list{index},'var')
        shear_count = shear_count + 1;
        shear_list{shear_count} = shear_wish_list{index};
        plot_gradient_wish_list{shear_count} = ...
            ['sh' num2str(shear_count) '_HP'];
        legend_gradient_wish_list{shear_count} = ...
            ['sh_' num2str(shear_count) ' HP'];
        junk = eval(shear_wish_list{index});
        SH = [SH junk];
    end
end

gradT_wish_list = {'gradT', 'gradT1', 'gradT2'};
gradT_list = [];
T_count = 0;
GRADT = [];
for index = 1:length(gradT_wish_list)
    if exist(gradT_wish_list{index},'var')
        T_count = T_count + 1;
        gradT_list{T_count} = gradT_wish_list{index};
        plot_gradient_wish_list{shear_count + T_count} = ...
            [gradT_wish_list{index} '_LP'];
        legend_gradient_wish_list{shear_count + T_count} = ...
            ['\partialT_' num2str(T_count) '/\partialx LP'];
        junk = eval(gradT_wish_list{index});
        GRADT = [GRADT junk];
    end
end

gradC_wish_list = {'gradC', 'gradC1', 'gradC2'};
gradC_list = [];
C_count = 0;
GRADC = [];
for index = 1:length(gradC_wish_list)
    if exist(gradC_wish_list{index},'var')
        C_count = C_count + 1;
        gradC_list{C_count} = gradC_wish_list{index};
        plot_gradient_wish_list{shear_count + T_count + C_count} = ...
            gradC_wish_list{index};
        legend_gradient_wish_list{shear_count + T_count + C_count} = ...
            ['\partialC_' num2str(C_count) '/\partialx'];
        junk = eval(gradC_wish_list{index});
        GRADC = [GRADC junk];
    end
end

for index = 1:length(gradC_list)
        C_count = C_count + 1;
        plot_gradient_wish_list{shear_count + T_count + C_count} = ...
            [gradC_wish_list{index} '_LP'];
        legend_gradient_wish_list{shear_count + T_count + C_count} = ...
            ['\partialC_' num2str(C_count) '/\partialx LP'];
end

spikes_1 = []; % they must exist
spikes_2 = [];
spikes_3 = [];
spikes_4 = [];

despike_string = ...
    ', despike_sh_1(1), despike_sh_1(2), fs_fast, round(despike_sh_1(3)*fs_fast));';

for index = 1:length(shear_list)
        if despike_sh_1(1) ~= inf
        eval(['[SH(:,index), spikes_' num2str(index) '] = ' ...
            'despike(SH(:,index)' despike_string]);
        end
end

% Despike the micro-conductivity signals.
spikes_C_1 = [];
spikes_C_2 = []; % assume no more than 2 micro-conductivity channels
despike_string_C = ...
    ', despike_C(1), despike_C(2), fs_fast, round(despike_C(3)*fs_fast));';
for index = 1:length(gradC_list)
        if despike_C(1) ~= inf
            eval(['[GRADC(:,index), spikes_C_' num2str(index) '] = ' ...
                'despike(GRADC(:,index)' despike_string_C]);
        end
end

% Apply gentle high-pass filter
[bh,ah] = butter(1, HP_cut/(fs_fast/2), 'high'); % HP at ~1 cpm
[bl,al] = butter(4, LP_cut/(fs_fast/2));% for profile plotting only

% High-pass the shear data with sensible initial conditions (filtfilt is
%       not good enough with its initial conditions).
SH_HP = filtfilt_odas(bh, ah, SH, 'high');
% Low-pass the shear data with sensible initial conditions.
SH_BP = filtfilt(bl, al, SH_HP);

% Low-pass filter the scalar gradients.
GRADT_LP = [];
GRADC_LP = [];
if T_count>0, GRADT_LP = filtfilt(bl, al, GRADT);end
if C_count>0, GRADC_LP = filtfilt(bl, al, GRADC);end

fig_num = fig_num + 1;
figure(fig_num); clf

% plot_gradient_wish_list = {'sh1_BP', 'sh2_BP', 'sh3_BP', 'sh4_BP', ...
%     'gradT1_LP', 'gradT2_LP',  'gradT_LP', ...
%     'gradC1', 'gradC2', ...
%     'gradC1_LP', 'gradC2_LP'};
% legend_gradient_wish_list = {'sh_1 BP', 'sh_2 BP', 'sh_3 BP', 'sh_4 BP', ...
%     '\partialT_1/\partialz LP', '\partialT_2/\partialz LP', '\partialT/\partialz LP',...
%     '\partialC_1/\partialz', '\partialC_2/\partialz', ...
%     '\partialC_1/\partialz LP','\partialC_2/\partialz LP' };
plot_gradient_list = [];

plot_scalar_wish_list = {'T_hres', 'T1_hres', 'T2_hres', 'C1_hres', 'C2_hres'};
legend_scalar_wish_list = {'T', 'T_1','T_2', 'C_1','C_2' };
plot_scalar_list = [];

index_gradient = 0;
for k = 1:length(plot_gradient_wish_list)
%    if exist(plot_gradient_wish_list{k},'var')
        index_gradient = index_gradient + 1;
        plot_gradient_list{index_gradient} = plot_gradient_wish_list{k};
        legend_list{index_gradient} = legend_gradient_wish_list{k};
%    end
end

% assume all are fast channels have the same sampling rate
%X_gradient_data = zeros(length(sh1),index_gradient);
x_gradient_offset = 0;
x_gradient_step   = 20;

LL = index_gradient; % the number of gradient channels

index_scalar = 0;
for k = 1:length(plot_scalar_wish_list)
    if exist(plot_scalar_wish_list{k},'var')
        index_scalar = index_scalar + 1;
        plot_scalar_list{index_scalar} = plot_scalar_wish_list{k};
        legend_list{LL + index_scalar} = legend_scalar_wish_list{k};
    end
end
% assume all are fast channels with same sampling rate
X_scalar_data = zeros(length(sh1),index_scalar);
x_scalar_offset = 0;

X_gradient_data = [SH_HP GRADT_LP GRADC GRADC_LP];
for k = 1: index_gradient
%   eval(['X_gradient_data(:,k) = '  plot_gradient_list{k} ';'])
   X_gradient_data(:,k) = X_gradient_data(:,k) + x_gradient_offset;
   x_gradient_offset = x_gradient_offset + x_gradient_step;
end

for k = 1: index_scalar
   eval(['X_scalar_data(:,k) = '  plot_scalar_list{k} ';'])
   X_scalar_data(:,k) = X_scalar_data(:,k) + x_scalar_offset;
   x_scalar_offset = x_scalar_offset + 0;
end

if (isempty(X_gradient_data) && isempty(X_scalar_data))
    error ('no vertical profile')
end

X_data = [];
if ~isempty(X_gradient_data)
    X_data = X_gradient_data;
end
if ~isempty(X_scalar_data)
    X_data = [X_data X_scalar_data];
end

plot(t_fast_YD - YD_0, X_data); grid
%set(gca,'ylim',[-15 25])
legend(legend_list, 'location','eastoutside')
xlabel(xlabel_string)
title(title_string)
junk = ['BP = ' num2str(HP_cut) '--' num2str(LP_cut) ' Hz'];
title([title_string junk])

if eps
    print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
    %saveas(fig_num, print_file)
    fig2pdf(gcf,    print_file, [], '')
end
%
%----------------------------------------------------------------
% Frequency spectra
fig_num = fig_num + 1;
figure(fig_num); clf

fft_num = round(fft_length*fs_fast);% fft length in units of samples

range = find((t_fast > t_start) & (t_fast <= t_end));
if isempty(range)
    warning(['The time range of ' num2str(t_start) ' to ' num2str(t_end) ...
        ' was not found in this profile'])
else

t_start = round(t_fast(range(1))); % So the labels give the actual range of data used
t_end   = round(t_fast(range(end)));
Data_fast = []; % not need for the single display spectrum
Data_slow = [];

%info
info.fft_length    = fft_num;
info.diss_length   = length(range);
info.overlap       = info.diss_length;
info.fs_fast       = fs_fast;
info.fs_slow       = fs_slow;
info.speed         = speed_fast(range);
info.T             = T1_hres(range);
info.t             = t_fast(range);
info.P             = P_fast(range);
info.Data_fast     = Data_fast;
info.Data_slow     = Data_slow;
info.fit_order     = fit_order;
info.f_AA          = f_AA;
info.fit_2_isr     = fit_2_isr;
info.f_limit       = f_limit;

num_of_shear_probes = shear_count;

diss = get_diss_odas(SH(range,:), A(range,:), info);

mean_speed = diss.speed;
F = diss.F;
K = diss.K;

% get_diss returns wavenumber spectra. Convert to frequency spectra by
% dividing by the mean speed. The spectra are corrected for the wavenumber
% response of the shear probe.
for index = 1:num_of_shear_probes
    r = num2str(index); % put it into a string
    junk = ['P_' shear_list{index}]; % The name of the spectrum
    eval([junk '       = diss.sh      (:,' r ',' r ') / mean_speed;' ])
    eval([junk '_clean = diss.sh_clean(:,' r ',' r ') / mean_speed;' ])
end

if (exist('Ax','var') && exist('Ay','var') && ~exist('Az','var')) || ...
        (exist('APy','var') && exist('APz','var'))
    % we have piezo-accelerometers so scale-down the acceleration spectra
    % to bring the in line with the shear spectra.
    diss.AA = 1e-7*diss.AA;
end

P_Ax =diss.AA(:,1,1);
P_Ay =diss.AA(:,2,2);
if size(diss.AA,2) == 3, P_Az = diss.AA(:,3,3);end

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
            '[ scalar_vectors ' scalar_wish_list{index} '];'])
        tmp = setupstr(obj, name_wish_list{index}, 'diff_gain');  % Find diff_gain
        diff_gain = [diff_gain str2double(tmp{1})];% Assume it is numeric

    end
end
num_of_scalar_spectra = count;
if count >0, we_have_scalars = true; end
% the number of scalar vectors for which we have a spectrum


% if exist('gradT1','var') && size(gradT1,1)==size(t_fast,1)
%     count = count + 1;
%     scalar_vector_list{count} = 'gradT1';
%     scalar_vectors = [scalar_vectors gradT1(range)];
%     tmp = setupstr(obj, 'T1_dT1', 'diff_gain');  % Find diff_gain
%     diff_gain = [diff_gain str2double(tmp{1})];% Assume it is numeric
% end

scalar_info.fft_length    = fft_num;
scalar_info.spec_length   = length(range);
scalar_info.overlap       = 0;
scalar_info.fs            = fs_fast;
scalar_info.gradient_method = 'first_difference';
scalar_info.diff_gain     = diff_gain;
scalar_info.f_AA          = 98;

scalar_spectra = []; % in case there are no scalars for spectra
if we_have_scalars
    scalar_spectra = ...
    get_scalar_spectra_odas(scalar_vectors(range,:), ...
    info.P, info.t, info.speed, scalar_info);
    scalar_speed = scalar_spectra.speed;
end

for index = 1:num_of_scalar_spectra
    % extract the scalar spectra
    junk = ['P_' scalar_vector_list{index}];
    eval([junk ' = scalar_spectra.scalar_spec(:,index) / scalar_speed;']);
end

epsilon = 1e-8*[1 10 100 1e3 1e4]; %
if exist('T_hres','var'),  nu = visc35(mean(T_hres (range)));end
if exist('T1_hres','var'), nu = visc35(mean(T1_hres(range)));end
[Pn, kn]=nasmyth(epsilon, nu);
fn = kn*mean_speed; Pn = Pn/mean_speed;

plot_wish_list = {...
    'P_Ax', 'P_Ay', 'P_Az', ...
    'P_sh1', 'P_sh2', 'P_sh3', 'P_sh4', ...
    'P_gradT', 'P_gradT1', 'P_gradT2', ...
    'P_gradC1', 'P_gradC2'};
legend_wish_list = {...
    'A_x', 'A_y', 'A_z', ...
    'sh_1', 'sh_2', 'sh_3', 'sh_4', ...
    '\partialT/\partialz', '\partialT_1/\partialz', ...
    '\partialT_2/\partialz', '\partialC_1/\partialz', '\partialC_2/\partialz'};

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

h = loglog(F, Y_data, fn, Pn, 'k');grid
hh=legend(legend_list,'location','NorthEastOutside');

for index = 1: size(Y_data,2)
    set(h(index), 'linewidth',1.5)
end

set(h(1),'linewidth',3,'color','b'); % Ax
set(h(2),'linewidth',3,'color',[0 0.5 0]); % Ay

f_upper = max(F);

y_limit = [1e-8 1e2]; % limits for y-axis
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
    [ num2str(t_start) ' < t < ' num2str(t_end) ' s, f_{HP} = ' num2str(HP_cut) ' Hz'];
title(new_title_string)
if eps
    print_file = ...
        [fpath filesep N '_P_' profile_num_string '_' sprintf('%02d',t_start) '_' ...
        sprintf('%02d',t_end) '_Fig_' sprintf('%02d',fig_num)];
    %saveas(fig_num, print_file)
    fig2pdf(gcf,    print_file, [], '')
end
%
%_____________________________
% Wavenumber Spectra
fig_num = fig_num + 1;
figure(fig_num); clf

% Extract wavenumber spectra of shear
for index = 1:num_of_shear_probes
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
for index = 1:num_of_scalar_spectra
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
for index = 1:num_of_shear_probes
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
    {'P_sh3_clean', '\partial u_3/\partial z_{clean}',  green,  3   }
    {'P_sh3',       '\partial u_3/\partial z',          green,  1.5 }
    {'phi3',        e3_string,                          'k',    1.5 }
    {'P_sh4_clean', '\partial u_4/\partial z_{clean}',  'm',    3   }
    {'P_sh4',       '\partial u_4/\partial z',          'm',    1.5 }
    {'phi4',        e4_string,                          'k',    1.5 }
    {'P_gradT',     '\partial T/\partial z ',           'c',    2   }
    {'P_gradT1',    '\partial T_1/\partial z ',         gold,   2   }
    {'P_gradT2',    '\partial T_2/\partial z ',         'm',    2   }
    {'P_gradC1',    '\partial C_1/\partial z ',         'c',    2   }
    {'P_gradC2',    '\partial C_2/\partial z ',         'c',    2   }};

point_wish_list = {
    % X pos         Y pos                           Label                                           Colour  Size
    {'K_max(1)', 'P_sh1_clean(K_max_index_1)', '[''k_{max} u_1='' num2str(round(K_max(1))) ''cpm'']',  'b',   18  }
    {'K_max(2)', 'P_sh2_clean(K_max_index_2)', '[''k_{max} u_2='' num2str(round(K_max(2))) ''cpm'']',  'r',   18  };
    {'K_max(3)', 'P_sh3_clean(K_max_index_3)', '[''k_{max} u_3='' num2str(round(K_max(3))) ''cpm'']',  green, 18  };
    {'K_max(4)', 'P_sh4_clean(K_max_index_4)', '[''k_{max} u_4='' num2str(round(K_max(4))) ''cpm'']',  'm',   18  }};

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
grid;

set(gca, 'ylim',[1e-6 1e1], 'xlim', [K(1) 300])

xlabel('\itk \rm [cpm]')
ylabel('\Phi (\itk\rm)  [s^{-2} cpm^{-1}]')

legend_list = {};
for plt = plot_list',  legend_list{end+1} = plt{1}{2}; end
for plt = point_list', legend_list{end+1} = plt{1}{3}; end
legend(legend_list,'location','NorthEastOutside');

new_title_string = [title_string ['Method = ' num2str(diss.method')], ...
        sprintf('t = %d - %d s, U = %0.2f m s^{-1}, f_{HP} = %0.2f Hz', ...
        t_start, t_end, mean_speed, HP_cut)];

title(new_title_string)

if eps
    print_file = ...
        [fpath filesep N '_P_' profile_num_string '_' sprintf('%02d',t_start) '_' ...
        sprintf('%02d',t_end) '_Fig_' sprintf('%02d',fig_num)];
    %saveas(fig_num, print_file)
    fig2pdf(gcf,    print_file, [], '')
end
end

%
% Now we calculate a dissipation profile for the entire file

fig_num = fig_num + 1;
figure(fig_num); clf

Data_fast = [];
Data_slow = [];

% Develop the lists of fast and slow vectors in this data file. We will
% only recognize coulmn vectors of the right length.
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

if exist('T1_hres','var'),
    info.T = T1_hres;
elseif exist(T2_hres','var')
    info.T = T2_hres;
else
    info.T = 10*ones(size(t_fast));
end

info.fft_length    = fft_num;
info.diss_length   = round(diss_length*fs_fast);
info.overlap       = round(overlap    *fs_fast);
info.fs_fast       = fs_fast;
info.fs_slow       = fs_slow;
info.speed         = speed_fast;
info.T             = T1_hres;
info.t             = t_fast;
info.P             = P_fast;
info.fit_2_isr     = fit_2_isr;
info.f_AA          = f_AA;
info.Data_fast     = Data_fast;
info.Data_slow     = Data_slow;
info.fast_list     = fast_list;
info.slow_list     = slow_list;


diss = get_diss_odas(SH_HP, A, info);

diss.fast_list  = fast_list;
diss.slow_list  = slow_list;
diss.header     = header(1,:);
diss.spikes_A1 = spikes_A1; % accelerometer spikes, save them for future use
diss.spikes_A2 = spikes_A2;
diss.spikes_1  = spikes_1; % shear probe spikes
diss.spikes_2  = spikes_2; % shear probe spikes
diss.spikes_3  = spikes_3; % shear probe spikes
diss.spikes_4  = spikes_4; % shear probe spikes


% produce spectra of scalar variable, if any.
scalar_info.spec_length   = info.diss_length;
if we_have_scalars
        scalar_spectra = ...
    get_scalar_spectra_odas(scalar_vectors, ...
    info.P, info.t, info.speed, scalar_info);
end

diss.scalar_spectra     = scalar_spectra; % lock into returned structure
diss.scalar_vector_list = scalar_vector_list;
diss.scalar_info        = scalar_info;


%
% Now extract dissipation rates and spectra for plotting
e = diss.e;

% Find the column that has the year-day for each dissipation estimate
index = find(strcmp('t_fast_YD', diss.fast_list));
if ~isempty(index)
    t_YD = diss.Data_fast(index(1),:); % use (1) in case of multiple matches.
else
    t_YD = diss.t; % Plots axis will not look good but we shoud get plots
end

h=  semilogy(t_YD - YD_0, e', '-o', 'linewidth', 1.5); grid on
my_colors = {'b', green, 'r', 'c'};

for index = 1:num_of_shear_probes
    set(h(index), 'markersize', 10, 'markerfacecolor', my_colors{index});
end

y_limits = get(gca,'ylim');
if log10(y_limits(2) / y_limits(1)) < 1
    y_limits(1) = 10^(floor(log10(y_limits(1))));
    y_limits(2) = 10^( ceil(log10(y_limits(2))));
    set(gca,'ylim',y_limits)
end

xlabel (xlabel_string)
ylabel ('\epsilon [W kg^{-1}]')

for index = 1:size(diss.e,1)
    legend_str{index} = ['\epsilon_' num2str(index)];
end

legend (legend_str, 'location','NorthEastOutside');

title(title_string)

if eps
    print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
    %saveas(fig_num, print_file)
    fig2pdf(gcf,    print_file, [], '')
end
diss.ql_info = ql_info;
result = diss;

%
% Now a wavenumber spectragram of shear sh1 and sh2
for index = 1:shear_count
    junk = squeeze(diss.sh_clean(:,index,index,:));
    eval(['sh' num2str(index) '_spec = junk;'])
end

if shear_count > 1
    fig_num = fig_num + 1;
    figure(fig_num); clf

    pcolor(t_YD - YD_0,  diss.K, log10(sh1_spec));grid
    set(gca,'ylim',[-150 150])
    caxis([-6 0])
    title(title_string)
    ylabel('\it k \rm [cpm]')
    xlabel(xlabel_string)
    colorbar('location','eastoutside')
    shading flat

    hold on

    pcolor(t_YD - YD_0, -diss.K, log10(sh2_spec));grid
    shading flat

    h=plot(t_YD - YD_0, diss.K_max(1,:));
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(1,:) == 1)
        index = find(diss.method(1,:) == 1);
        h=plot(t_YD(index) - YD_0, diss.K_max(1,index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    h=plot(t_YD - YD_0, -diss.K_max(2,:));
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(2,:) == 1)
        index = find(diss.method(2,:) == 1);
        h=plot(t_YD(index) - YD_0, -diss.K_max(2,index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    x = t_YD(1) - YD_0;
    y = 140;
    text(x,  y, 'u_1', 'color', 'w', 'backgroundcolor','k');
    text(x, -y, 'u_2', 'color', 'w', 'backgroundcolor','k');

    hold off

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        print('-dpng', '-zbuffer', '-r300', print_file)
    end
end
%
% Now a wavenumber spectragram of shear sh3 and sh4
fig_num = fig_num + 1;
figure(fig_num); clf

if shear_count > 3
    pcolor(t_YD - YD_0,  diss.K, log10(sh3_spec));grid
    set(gca,'ylim',[-150 150])
    caxis([-6 0])
    title(title_string)
    ylabel('\it k \rm [cpm]')
    xlabel(xlabel_string)
    colorbar('location','eastoutside')
    shading flat

    hold on

    pcolor(t_YD - YD_0, -diss.K, log10(sh4_spec));grid
    shading flat

    h=plot(t_YD - YD_0, diss.K_max(3,:));
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(3,:) == 1)
        index = find(diss.method(3,:) == 1);
        h=plot(t_YD(index) - YD_0, diss.K_max(3,index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    h=plot(t_YD - YD_0, -diss.K_max(4,:));
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(4,:) == 1)
        index = find(diss.method(4,:) == 1);
        h=plot(t_YD(index) - YD_0, -diss.K_max(4,index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    x = t_YD(1) - YD_0;
    y = 140;
    text(x,  y, 'u_3', 'color', 'w', 'backgroundcolor','k');
    text(x, -y, 'u_4', 'color', 'w', 'backgroundcolor','k');

    hold off

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        print('-dpng', '-zbuffer', '-r300', print_file)
    end
end
%
% Now a spectragram of shear using frequency, sh1 and sh2
fig_num = fig_num + 1;
figure(fig_num); clf

if shear_count > 1
    % scale to preserve variance
    speed_2 = repmat(diss.speed, 1, size(sh1_spec,1))';

    for index = 1:shear_count
        eval(['sh' num2str(index) '_spec = sh' num2str(index) '_spec ./ speed_2;'])
    end

    pcolor(t_YD - YD_0, diss.F, log10(sh1_spec));grid
    shading flat

    set(gca,'ylim', 350*[-1 1])
    caxis([-6 0])

    ylabel('\it f \rm [Hz]')
    xlabel(xlabel_string)
    title(title_string)

    colorbar('location','eastoutside')

    hold on

    pcolor(t_YD - YD_0,  -diss.F, log10(sh2_spec));grid
    shading flat

    h=plot(t_YD - YD_0, diss.K_max(1,:) .* diss.speed');
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(1,:) == 1)
        index = find(diss.method(1,:) == 1);
        h=plot(t_YD(index) - YD_0, diss.K_max(1,index) .* diss.speed(index)', '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    h=plot(t_YD - YD_0, -diss.K_max(2,:).* diss.speed');
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(2,:) == 1)
        index = find(diss.method(2,:) == 1);
        h=plot(t_YD(index) - YD_0, -diss.K_max(2,index).* diss.speed(index)', '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    x = t_YD(1) - YD_0;
    y = 140;
    text(x,  y, 'u_1', 'color', 'w', 'backgroundcolor','k');
    text(x, -y, 'u_2', 'color', 'w', 'backgroundcolor','k');

    hold off

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        print('-dpng', '-zbuffer', '-r300', print_file)
    end
end

%
% Now a spectragram of shear sh3 and sh4 using frequency
fig_num = fig_num + 1;
figure(fig_num); clf

if shear_count > 3
    pcolor(t_YD - YD_0,  diss.F, log10(sh3_spec));grid
    shading flat

    set(gca,'ylim', 350*[-1 1])
    caxis([-6 0])
    title(title_string)
    ylabel('\it f \rm [Hz]')
    xlabel(xlabel_string)
    colorbar('location','eastoutside')

    hold on

    pcolor(t_YD - YD_0, -diss.F, log10(sh4_spec));grid
    shading flat

    h=plot(t_YD - YD_0, diss.K_max(3,:).* diss.speed');
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(3,:) == 1)
        index = find(diss.method(1,:) == 1);
        h=plot(t_YD(index) - YD_0, diss.K_max(3,index).* diss.speed(index)', '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    h=plot(t_YD - YD_0, -diss.K_max(4,:).* diss.speed');
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(4,:) == 1)
        index = find(diss.method(4,:) == 1);
        h=plot(t_YD(index) - YD_0, -diss.K_max(4,index).* diss.speed(index)', '*');
        set(h(1),'linewidth', 2, 'color','w')
    end

    x = t_YD(1) - YD_0;
    y = 140;
    text(x,  y, 'u_3', 'color', 'w', 'backgroundcolor','k');
    text(x, -y, 'u_4', 'color', 'w', 'backgroundcolor','k');

    hold off

    if eps
        print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
        %saveas(fig_num, print_file)
        print('-dpng', '-zbuffer', '-r300', print_file)
    end
end
%
% Now a spectragram of acceleration
fig_num = fig_num + 1;
figure(fig_num); clf

Ax_spec = squeeze(diss.AA(:,1,1,:));
Ay_spec = squeeze(diss.AA(:,2,2,:));

pcolor(t_YD - YD_0, diss.F, log10(Ax_spec));grid
shading flat

set(gca,'ylim', 350*[-1 1])
 caxis([1 9])

ylabel('\it f \rm [Hz]')
xlabel(xlabel_string)
title(title_string)

colorbar('location','eastoutside')

hold on

pcolor(t_YD - YD_0,  -diss.F, log10(Ay_spec));grid
shading flat

x = t_YD(1) - YD_0;
y = 140;
text(x,  y, 'A_z', 'color', 'w', 'backgroundcolor','k');
text(x, -y, 'A_y', 'color', 'w', 'backgroundcolor','k');

hold off

if eps
    print_file = [fpath filesep N '_P_' profile_num_string '_Fig_' sprintf('%02d',fig_num)];
    %saveas(fig_num, print_file)
    print('-dpng', '-zbuffer', '-r300', print_file)
end

% That is all fokes. :-)

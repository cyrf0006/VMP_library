function varargout = parameter_list(instrument,varargin);

% PARAMETER_LIST - Return values of constants/parameters for ODAS or MSS data processing
% [X1,X2,...] = parameter_list(instrument,'X1','X2',....);
% 
% Given the instrument type (one of 'default', 'odas', or 'mss'), return the value of
% the desired constants and parameters.
%
% Examples of usage:
% Fs = parameter_list('mss','Fs');
% [coef,deconvolve] = parameter_list('odas','coef','deconvolve');
%
% Isabelle Gaboury, 10 Mar. 2005
% 28 Mar. 2005 - trimmed down

if isempty(instrument), instrument = 'default'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Definitions of various constants, instrument-specific settings, %%%%
%%%%%% coefficients, etc.  Divided into various sections, based on the %%%%
%%%%%% type of instrument involved.                                  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants, used by all instruments %%

grav = 9.81;        % gravitational constant
karman = 0.4;       % van Karman constant (may be changed as desired)
mixeff = 0.2;       % mixing efficiency (gamma) (may be changed as desired)
twopi = 6.28318530717959;   % value of 2*pi

%% Values specific to ODAS instruments %%

if strcmpi(instrument,'odas')
    
    % Fast sampling frequency
    Fs = 512;
    
    % Default coefficients for compressibility and thermal expanion (from
    % Sea-Bird calibration sheet, May 2004)
    CPcor = -9.5700e-8; 
    CTcor = 3.2500e-6; 
    
    % Differentiator gains (for use by deconvolve routines)
    T_gain = 1;
    C_gain = 5/pi;
    P_gain = 20.6;  % 33 for newer instruments
    
    % Default channel numbers and names for vertical profilers
    ch_nums = [0:12 16:19 255]; 
    ch_names = {'gnd','Ax','Ay','Az','T1','T1_dT1','T2','T2_dT2','sh1','sh2',...
        'P','P_dP','C_dC','sbtE','sbtO','sbcE','sbcO','sp_char'};
    
    % Locations of various values in the header matrix:
    ii_nrows = 31;      % # of rows
    ii_profiler = 63; % profiler flag
    header_length = 64; % assumed length of the header
    
    % Default calibration coefficients:
    % Temperature calibration coefficients: g,h,i,j,f_0,f_ref,p_count_set
    coef.sbt = [4.34248180e-3,6.39925941e-4,2.13053636e-5,1.79260071e-6,1000.0,24e6,128]; % SBE3 SN-4439, calibrated 2004-09-21
    % Conductivity calibration coefficients: g,0,h,i,j,f_ref,p_count_set
    coef.sbc = [-1.01419820e1,0,1.40329950e0,-3.70968211e-5,6.32783997e-5,24e6,128];  % SBE4 SN-2847, calibrated 2004-09-04
    % Pressure calibration coefficients: polynomial coefficients
    coef.P = [4.03,0.049997,1.344e-8];        % Hartmut Peters instrument
    coef.P_dP = [4.03,0.049997,1.344e-8];     % Hartmut Peters instrument
    % Accelerometer calibration coefficients: bias, sensitivity
    coef.ax = [256.5,13146.5];        % Hartmut Peters instrument (may be wrong, &/or be linear coefficients)
    coef.ay = [586.0,13166.0];
    coef.az = [0,1];
    % Ground voltage is not converted
    coef.gnd = [0,1];
    % Propeller flowmeter: 4 polynomial coefficients, 0, scale (currently
    % ignored), fref, period. Note that this is not generally used; will
    % probably have to be updated when brought back into service.
    coef.u = [0.001,3.389e-2,0,0,0,256,1e7,16];
    
    % Calibration coefficients for shear processing
    ADC = [5/2^16 0]; % Polynomial coefficients for A/D counts-volts conversion (slope, offset)
    she_gain = [1.00 1.00];     % amplifier gain
    she_cal = [0.0835 0.0828];  % shear probe sensitivities (in channel order);
    
%% Values specific to MSS instruments %%

elseif strcmpi(instrument,'mss')
    
    % Sampling frequency
    Fs = 1024;
    
    % First few letters of "dummy" channels
    dummylabel = 'dum';
    
    % Names of the columns in the .mrd data file and in the .prb setup file
    col_names = {'data_type','profiler_id','sensor_num','calc_type','var_name',...
        'units','coef'};
    col_names_prb = {'sens_no','calc_type','var_name','units','coef'};
    
    % Buffer region used by msscleanup.m & rmdrop.m
    rmdrop_buf = 2;
    
    % Polynomial coeffs for A/D counts to volts conversion
    ADC = [9.155423e-5 4.577707e-5]; 
    
    % Calibration coefficients for shear processing
    she_gain = 11;      % amplifier gain
    S186 = 0.0856;      % shear probe sensitivities
    S189 = 0.0765;
    S191 = 0.0784;
    she_cal = [S186 S189];  % shear calibration coefficients (in channel order)
         
%% Default values, used by routines that are not specificaly written %%
%% for one of the two above instruments                              %%

elseif strcmpi(instrument,'default')
    
    % Sampling frequency
    Fs = 1024;
               
else
    error('Unknown instrument type');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Evaluate the parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(varargin)
    varargout{ii} = eval(eval(['varargin{' num2str(ii) '};']));
end
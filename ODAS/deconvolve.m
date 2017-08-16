%% deconvolve
% Deconvolve signal plus derivative to produce high resolution data
%%
% <latex>\index{Type A!deconvolve}</latex>
%
%%% Syntax
%   xHires = deconvolve( dataType, X, XdX, fSample, gain, ver, xmpSerNum )
%
% * [dataType] String describing the data type to be deconvolved. For legacy 
%       ODAS and XMP, it can be one of 'pressure', 'temperature', or 
%       'conductivity'.  For ODAS v6 and up, this string identifies the channel
%       to be deconvolved - a channel containing a 'diffGain' key.  The channel
%       is identified through the 'name' parameter or, when not found, the
%       section title.
% * [X] Low-resolution signals.  Can be empty for data types other than 
%       pressure.
% * [XdX] Pre-emphasized signals.
% * [fSample] Sampling rate of the data in Hz.
% * [gain] Differentiator gain as either a numeric value or as the configuration 
%       string from which it can be found.
% * [ver] ODAS header version.
% * [xmpSerNum] String identifying the XMP serial number.  Only used with XMP
%       instruments. Should be of the form 'xmp_nnnnn' where nnnnn is the 
%       instrument's serial number with leading zeros. For example, 'xmp_00012'.
% * []
% * [xHires] Deconvolved, high-resolution signal.
%
%%% Description
% Deconvolve vector of pre-emphasized data (temperature, conductivity, or 
% pressure) to yield high-resolution data. The pre-emphasized signal 
% (x+gain*dx/dt) is low-pass filtered using appropriate initial conditions 
% after the method described in Mudge and Lueck, 1994.
%
%%% Examples
%
%    >> pHires = deconvolve( 'pressure', P, P_dP, 32, 20.5 )           (legacy)
%
%    >> tHires = deconvolve( 'temperature', [], T_dT, 512, 1)          (legacy)
%
%    >> pHires = deconvolve( 'pressure', P, P_dP, fs_slow, setupfilestr, 
%                            header_version, xmp_ser_num)            (xmp only)
%
%    >> t1Hires = deconvolve( 'T1_dT1', [], T1_dT1, fs_fast, setupfilestr, 
%                             header_version);                 (odas v6 and up)
%
% For temperature and conductivity, the initial conditions are estimated from 
% the pre-emphasized signal only. Therefore, the vector X does not have to be
% passed to the function and can be given by '[ ]'.
%
% For pressure, you must pass both 'X' and 'X_dX' to this function. Both 
% vectors are needed to make a good estimate of the initial conditions for 
% filtering. Initial conditions are very important for pressure because the 
% gain is usually ~20 seconds, and errors caused by a poor choice of initial
% conditions can last for several hundred seconds! In addition, the
% high-resolution pressure is linearly adjusted to match the low-resolution
% pressure so that the factory calibration coefficients can later be used to
% convert this signal into physical units.
%
% The gains for all signal types is provided in the calibration report of 
% your instrument.
%
% Mudge, T.D. and R.G. Lueck, 1994: Digital signal processing to enhance 
% oceanographic observations, _J. Atmos. Oceanogr. Techn._, 11, 825-836.
%
% @image @images/pressure_deconvolution1 @Deconvolution example. @The green 
% curve is from the normal pressure channel, and the blue curve is derived from 
% the enhanced pressure channel. This data is from a profiler that has impacted 
% the soft bottom of a lake. Both signals are shown with full bandwidth (0 - 32 Hz)
% without any smoothing. The full-scale of the pressure transducer is 500 dBar.
%
% @image @images/pressure_deconvolution2 @Deconvolution example two. @Same 
% as previous figure but with zoom-in on only the high-resolution pressure. Again,
% full bandwidth without any smoothing.
%
% @image @images/pressure_deconvolution3 @The rate of change of pressure derived 
% from the normal pressure signal and the high-resolution pressure signals using 
% the gradient function. @Full bandwidth signals. The normal pressure signal 
% produces a fall-rate that is useless without extensive smoothing. It even has 
% negative values, which means that the normal pressure record is not even 
% monotonic. The rate of change of the high-resolution pressure is smooth, always 
% positive, and the high-resolution pressure, itself, can be used directly for 
% plotting other variables as a function pressure (or depth). Notice that the 
% high-resolution rate of change of pressure has been multiplied by 10 for visual 
% clarity. The fall-rate is about 0.17 m/s.

% *Version History:*
%
% * 2005-02-08 (IG) based on OTL/RGL routines pressure_all, thermistor_all, 
%        pressure_vel
% * 2010-01-14 (AWS) odas v6 related changes
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2010-02-28 (RGL) reverted back to simple initial conditions for the case of 
%        pressure. The previous method was not right and can give very erroneous 
%        results for the case of a pressure transducer on a near-surface mooring 
%        with significant waves.
% * 2012-04-11 (WID) calls to inifile_with_instring changed to setupstr
% * 2012-11-05 (WID) documentation update

function X_hires = deconvolve(data_type,X,X_dX,f_s,arg5,ver,xmp_ser_num)

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5                 % Default differentiator gains
    if strcmp(data_type,'pressure'),         diff_gain=20.5; 
    elseif strcmp(data_type,'temperature'),  diff_gain = 1;
    elseif strcmp(data_type,'conductivity'), diff_gain = 1;  % (manufacturer's value)
    else error('Unknown data type');
    end
else
    if ischar(arg5)
        %we are dealing with the setup file str being passed for arg5
        cfg = setupstr(arg5);

        if nargin < 6
            error('odas header version is needed');
        end

        header_version = ver;

        if(header_version < 6)
            error(['incorrect header version: ' header_version ]);
        end

        tmp = setupstr(cfg, '', 'xmp');

        if(~isempty(tmp)) % We have an XMP
            if(nargin < 7)
                error('no xmp serial number provided in the form xmp_nnnnn, where nnnnn is the instrument ''s serial number with leading zeroes');
            end
            switch data_type
                case 'pressure'
                    diff_gain_key = 'P_dP_diff_gain';
                case 'temperature'
                    diff_gain_key = 'T_dT_diff_gain';
            end
            
            tmp = setupstr(cfg, xmp_ser_num, diff_gain_key);
            if(isempty(tmp))
                error(['cannot find the following key: ' diff_gain_key '']);
            else
                diff_gain = str2double(tmp);
            end

        else
            diff_gain_key = 'diff_gain';
            
            tmp = setupstr(cfg, data_type, diff_gain_key);
            if(isempty(tmp))
                error(['cannot find the following key: ' diff_gain_key '']);
            else
                diff_gain = str2double(tmp);
            end          
        end
        
    else
        diff_gain = arg5;
    end
end



f_c = 1/(2*pi*diff_gain);   % Cut-off frequency for filtering
if strcmp(data_type,'pressure') || strcmpi(data_type,'dpres')
   t_ave = 8; % Averaging time scales (for initial conditions)
else 
    t_ave = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Prepare for filtering: coefficients, initial conditions %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = butter(1,f_c/(f_s/2));                          % Low-pass filtering coeffs.
if strcmp(data_type,'pressure') || strcmpi(data_type,'dpres')    % For the pressure, take into account the rate of change of pressure (slope)
    p = polyfit(X,X_dX,1);              %   from the low-resolution record
    if p(1)<-0.5, X_dX=-X_dX; end       % A few instruments have X_dX inverted with respect to X
%    N = (1:t_ave*f_s)'; t = (N-1)/f_s;
%    p = polyfit(t(N),X(N),1);                 
%    z = filtic(b,a,mean(X_dX(N))-diff_gain*p(1)-p(1)*t_ave/2,X_dX(1));
    z = filtic(b,a,X(1),X_dX(1));
else
%    z = filtic (b,a,mean(X_dX(1:t_ave*f_s)),X_dX(1));   % For temperature, initial conditions based on first portion of the record
    z = filtic (b,a,mean(X_dX(1:round(t_ave*f_s))),X_dX(1));   %For t/c, initial conditions based on first portion of the record
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Deconvolve to obtain the high-resolution data %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_hires= filter(b,a,X_dX,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% For the pressure, regress against the low-resolution vector to %%%%%
%%%%%% remove the small offset in the derivative circuit.             %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(data_type,'pressure') || strcmpi(data_type,'dpres')
    p = polyfit(X_hires,X,1);
    X_hires = polyval(p,X_hires);
end

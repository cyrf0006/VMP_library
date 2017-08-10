%% convert_odas
% Convert from raw counts to physical units
%%
% <latex>\index{Type A!convert\_odas}</latex>
%
%%% Syntax
%  [phys, units] = convert_odas(var, sensor, coefSource, coef, ver, xmpSerNum)
%
% * [var] ODAS variable(s) in counts.  Either a single vector or a matrix with a 
%       single variable in each column.
% * [sensor] String array of names for the included ODAS variables. Names should
%       match the channel names given in the setup / configuration file.
% * [coefSource] Location of the coefficients (optional).  Set to 'default' (use
%       values within the routine), 'file' (load from setup file) or 'string' 
%       (load from configuration string).  Most v6 users should use 'string' 
%       while v0 users will use 'file' or 'default'.
% * [coef] Coefficient values.  When coefSource is set to 'file', this is the 
%       name of the setup file.  When coef_source is 'string', this is the
%       variable that contains the configuration string (e.g., 'setupfilestr').
% * [ver] Header version - always 1 for ODAS versions prior to v6.
% * [xmpSerNum] String identifying the XMP serial number.  Only used with XMP
%       instruments. Should be of the form 'xmp_nnnnn' where nnnnn is the 
%       instrument's serial number with leading zeros. For example, 'xmp_00012'.
% * []
% * [phys] Variable(s) converted into physical units.
% * [units] String representation of the converted variable units.
%
%%% Description
% Converts ODAS data from recorded raw counts into physical units.  Requires
% data as raw counts (var), the channel type used to collect the data (sensor),
% and the calibration coefficients needed to perform the conversion.  For v6 
% data files, the calibration coefficients are found within the configuration 
% string embedded within the data file.  For previous v1 data files, one can use
% the default calibration coefficients found within this function, or
% calibration coefficients contained within a setup file.  However, only a
% limited number of signals can be converted when using v1 data files.
%
% Returns the data converted into physical units as double values along with the
% observed units as string values.
%
%%% Examples (ODAS >= v6)
%
%    >> T2 = convert_odas(T2, 't2', 'string', setupfilestr, 6.0)
%
% Convert data within 'T2' using coefficients from section 't2' found within
% the configuration file 'setupfilestr'.  Notice that the second argument, 't2',
% will work with or without capital letters.
%
%    >> P = convert_odas(p, 'P', 'string', setupfilestr, 6.0, 'xmp_00012')
%
% Convert a pressure vector from an XMP style instrument into physical units.
%
% For more examples and detailed information regarding conversion from Raw data
% to physical units, please see section named 'Converting into Physical Units'.
%
%%% Examples (ODAS == v1)
%
%    >> SBT1 = convert_odas( sbt1, 'sbt' )
%
% Convert vector sbt1 into physical units using the Seabird temperature formula.
% Coefficients are not supplied so default values from within the function are
% used.
%
%    >> [SBC1, cUnits] = convert_odas( sbc1, 'sbc', 'file', 'setup.txt' )
%
% Convert vector sbc1 into physical units using the Seabird conductivity 
% formula.  Coefficients will be extracted from the provided 'setup.txt' file.
% Both the converted values and the units of these values are returned.
%
%    >> [SBT1,SBC1,units] = convert_odas([sbt1 sbc1], {'sbt','sbc'}, 'default')
%
% Multiple vectors can be supplied and converted in one call. The input vectors 
% '[sbt1 sbc1]' are processed using their respective 'sbt' and 'sbc' types.  
% Default coefficients from the convert_odas.m file are used for the conversion.

% *Version History:*
%
% * 2005-02-15 (IG) Based on RGL/OTL routines including plot_tomi, sb_t, sb_c, 
%               accel_all, and pressure_all.
% * 2005-02-21 (IG) now case-insensitive
% * 2005-05-01 (IG) corrected a small typo in the coefficient loading section
% * 2005-07-12 (IG) Provided an option to find the start of coefficient lines 
%               that will work with pre-R14 versions of MATLAB.  Made commas the
%               only acceptable delimiters
% * 2006-04-17 (RGL) Added 3-axis magnetometer section.
% * 2007-06-16 (RGL) Added recognition for Alec Electronics EM Current meter, 
%               Altimeter and Fluorometers
% * 2007-11-05 (RGL)  Added SBE43F oxygen sensor. Only gives the frequency and 
%               not the proper disolved oxygen concentration.
% * 2010-01-14 (AWS) odas v6 related changes
% * 2010-09-24 (AWS) added adis function, inclinometer x, y, t conversion
% * 2010-10-04 (AWS) added documentation tags for matlab publishing
% * 2010-12-03 (AWS) added conversion for thermistor and micro conductivity
% * 2011-03-24 (AWS) changed conversion for voltage type to include sensor Gain
% * 2011-11-19 (RGL) Corrected voltage type so that it divides by the gain 
%               instead multiplying by the gain.
% * 2012-03-04 (RGL) big changes to make this whole damn process more
%               rational. You can now use a sensor (or channel) name to 
%               find the type and convert it to physical units. It is no 
%               longer necessary to know the name of the section in which
%               the channel and its coefficients are identified. Brother!!
% * 2012-03-21 (RGL) changed thermistor conversion to avoid taking log of
%               negative numbers if thermistor is broken.
% * 2012-03-23 (RGL) changed the way micro-conductivity is calculated.
% * 2012-03-28 (RGL) Expanded the XMP section so that it can habdle all
%               channels.
% * 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% * 2012-04-30 (WID) fix for incl?? channel conversions with pre-v6 data files
% * 2012-05-07 (WID) update to documentation 
% * 2012-11-02 (WID) update to documentation 
% * 2013-02-26 (WID) support for beta_1 and beta_2 coefficients on therm channels

function varargout = convert_odas(X,sensor,coef_source, arg4, ver, xmp_ser_num)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Deal with the input parameters, set constants %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% By default, use the coefficients in this M-file
if nargin<3, coef_source = 'default'; end
if nargin<4, coef_file = 'setup.txt'; end
if nargin<5
    if strcmp(coef_source,'string')
        error('odas version is required');
    end
end

if exist('ver','var')
    header_version = ver;
else
    header_version = 1; %assume legacy odas version
end

% Check the sizes of the "X" and "sensor" variables
if ischar(sensor), sensor={sensor}; end
if size(X,1)>1 && size(X,2)>1   
    if size(X,2)~=length(sensor)
        error('The number of columns of X must match the length of the sensor list');
    end
end

% Constants
g = 9.81;               % Gravitational constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Deal with calibration coefficients. If using default calibration %%%
%%%%%% coefficients, the values listed in this M-file are used. These %%%%%
%%%%%% may be modified to reflect instrument setup. Note that changes %%%%%
%%%%%% to the names (other than changing the suffixes on the SeaBird %%%%%%
%%%%%% variables) will affect how they are calibrated; see the series %%%%%
%%%%%% of IF statements at the end of this section.                  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(coef_source,'default')    % Default calibration coefficients
    % Temperature calibration coefficients: g,h,i,j,f_0,f_ref,p_count_set
    coef.sbt = [4.34248180e-3,6.39925941e-4,2.13053636e-5,1.79260071e-6,1000.0,24e6,128]; % SBE3 SN-4439, calibrated 2004-09-21
    % Conductivity calibration coefficients: g,0,h,i,j,f_ref,p_count_set
    coef.sbc = [-1.01419820e1,0,1.40329950e0,-3.70968211e-5,6.32783997e-5,24e6,128];  % SBE4 SN-2847, calibrated 2004-09-04
    % Pressure calibration coefficients: polynomial coefficients
    coef.pres = [4.03,-0.049997,1.344e-8];        % Hartmut Peters instrument
    % Accelerometer calibration coefficients: bias, sensitivity
    coef.ax = [256.5,13146.5];        % Hartmut Peters instrument (may be wrong, &/or be linear coefficients)
    coef.ay = [586.0,13166.0];
    coef.az = [0,1];
    % Ground voltage is not converted
    coef.gnd = [0,1];
    coef.Mx = [0,0.1];
    coef.My = [0,0.1];
    coef.Mz = [0,0.1];
    coef.O2 = [0,0,0,0,0,24e6,128];

elseif strcmp(coef_source,'file')   % Load coefficients from file
    % Obtain the setup file name from the user if not provided in the inputs
    
    if isempty(arg4)
        coef_file = input('Enter the setup file name: ','s'); 
    else
        coef_file = arg4;
    end
    
    % Open the file
    fid_setup = fopen(coef_file,'r');
    if fid_setup<3
       error(sprintf('Unable to open the setup file %s\n',coef_file))
    end
    % Extract the data from the file
    done=0;
    while ~done
        junk=fgetl(fid_setup);
        if junk==-1, done=1;
        elseif (length(junk)>1&~strcmp(junk(1),'#'))        % ignore commented lines
            if ~isempty(findstr(junk,'channel'))            % get the calibration data
                if str2num(char(version('-release')))>=14
                    start_ii = min(find(isstrprop(junk,'digit')==1));  % Find the first number in the line (usually a channel number)
                else
                    for ii=1:length(junk)
                        isdigit(ii) = ~isempty(str2num(junk(ii)));
                    end
                    start_ii = min(find(isdigit==1));
                end 
%                 delim_ii = find(isspace(junk) | junk==',' | junk==';' ==1);
                delim_ii = find(junk==',' ==1); % find delimiters (only commas currently)
                delim_ii = [delim_ii(find(delim_ii>=start_ii)) length(junk)+1];
                for jj=1:length(delim_ii)-2
                    coef.(deblank(junk(delim_ii(1)+1:delim_ii(2)-1)))(jj) = str2num(junk(delim_ii(jj+1)+1:delim_ii(jj+2)-1)); 
                end
            end
        end
    end
    clear start_ii delim_ii
    fclose(fid_setup);   

elseif strcmp(coef_source, 'string')  %arg4 will be interpreted as the setup file string (ini file format)
    header_version = ver;
    cfg = setupstr(arg4);

    if(header_version < 6)
        error(['incorrect header version: ' header_version ]);
    end
    isXMP = 0;
    tmp = setupstr(cfg, '', 'xmp');
    if(~isempty(tmp))
        isXMP = 1;
        if(nargin < 6)
            error('no xmp serial number provided - cannot continue');
        end
    end

else
    error('Unknown coefficient source');
end

if (header_version < 6)
%    %% legacy odas version

    % Figure out what types of sensors has been provided. Note that this is not
    % case-sensitive.
    for ii=1:length(sensor)
        if length(sensor{ii})>=3&strcmpi(sensor{ii}(1:3),'sbt'), sensor2{ii} = 'sbt'; 
        elseif length(sensor{ii})>=3&strcmpi(sensor{ii}(1:3),'sbc'), sensor2{ii} = 'sbc'; 
        elseif strcmpi(sensor{ii},'pres')|strcmpi(sensor,'P')|strcmpi(sensor,'P_dP'), sensor2{ii} = 'pres'; 
        elseif length(sensor{ii})>=3&strcmpi(sensor{ii}(1:3),'gnd'), sensor2{ii} = 'gnd'; 
        elseif strcmpi(sensor{ii},'incl_x') || strcmpi(sensor{ii},'incl_y'), sensor2{ii} = 'inclxy'; 
        elseif strcmpi(sensor{ii},'incl_t'), sensor2{ii} = 'inclt'; 
        elseif strcmpi(sensor{ii},'ax')|strcmpi(sensor{ii},'ay')|strcmpi(sensor{ii},'az')|...
                strcmpi(sensor{ii},'pitch')|strcmpi(sensor{ii},'roll'), sensor2{ii} = 'accel'; 
        elseif strcmpi(sensor{ii},'Mx')|strcmpi(sensor{ii},'My')|strcmpi(sensor{ii},'Mz')|...
                strcmpi(sensor{ii},'compass')|strcmpi(sensor{ii},'mag'), sensor2{ii} = 'magn';
        elseif strcmpi(sensor{ii},'F_R')|strcmpi(sensor{ii},'F_N')|strcmpi(sensor{ii},'F_C')|...
                strcmpi(sensor{ii},'fluoro')|strcmpi(sensor{ii},'trans'), sensor2{ii} = 'fluorometer';
        elseif strcmpi(sensor{ii},'U')|strcmpi(sensor{ii},'V'), sensor2{ii} = 'Alec_EMC';
        elseif strcmpi(sensor{ii},'Alt')|strcmpi(sensor{ii},'Altimeter'), sensor2{ii} = 'altimeter';
        elseif strcmpi(sensor{ii},'Alt_error')|strcmpi(sensor{ii},'Altimeter_Error'), sensor2{ii} = 'altimeter_error';
        elseif strcmpi(sensor{ii},'Ux')|  strcmpi(sensor{ii},'Uy'),   sensor2{ii} = 'voltage';
        elseif strcmpi(sensor{ii},'Vbat')|strcmpi(sensor{ii},'v_bat'), sensor2{ii} = 'voltage';
        elseif strcmpi(sensor{ii},'O2'), sensor2{ii} = 'oxygen';
        else
            error(['Unknown sensor type: ' sensor{ii}]);
        end
    end

%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %%%%%% For each variable requested, extract the appropriate %%%%%%
%    %%%%%% coefficients and apply based on the sensor type.     %%%%%%
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii=1:length(sensor)

        % Extract the calibration coefficients from the structure created above
        % (case-insensitive)
        junk = fieldnames(coef);
        jj = find(strcmpi(fieldnames(coef),sensor{ii})==1);
        if ~isempty(jj), coeffs = coef.(junk{jj});
        else error(['Unknown sensor: ' sensor{ii}]);
        end

        switch sensor2{ii}

            case 'sbt'      % Sea-Bird temperatures
                units = 'deg C';
                f = coeffs(7)*coeffs(6)./X(:,ii);
                c = log(coeffs(5)./f);
                X_phys = 1./polyval(fliplr(coeffs(1:4)),c) - 273.15;

            case 'sbc'      % Sea-Bird conductivity
%                %%% Note that this version does not correct for thermal
%                %%% expansion/compressibility; this must be done separately.
                units = 'mS/cm';
                f = coeffs(7)*coeffs(6)./X(:,ii)/1000;
                X_phys = polyval(fliplr(coeffs(1:5)),f);

            case 'pres'     % Pressure
                units = 'dBar';
                X_phys = polyval(fliplr(coeffs(1:3)),X(:,ii));

            case 'gnd'      % Ground voltage
                units = 'V';
                X_phys = polyval(fliplr(coeffs(1:2)),X(:,ii));

            case 'accel'    % Accelerometers
                units = 'm sec^-2';
                X_phys = g*(X(:,ii)-coeffs(1))/coeffs(2);

            case 'magn'    % Magnetometers
                units = 'micro-Tesla';
                X_phys =   (X(:,ii)-coeffs(1))/coeffs(2);

          case 'fluorometer'    % Fluoromter
                units = 'volts';
                X_phys = polyval(fliplr(coeffs(1:2)),X(:,ii));

            case 'Alec_EMC'    % Alec EM Current Meter
                units = 'm/s';
                X_phys = polyval(fliplr(coeffs(1:2)),X(:,ii));

            case 'altimeter'    % Altimeter signal
                units = 'm';
                X_phys = polyval(fliplr(coeffs(1:2)),X(:,ii));

            case 'inclxy'    % ADIS inclinometer x or y incline
                units = 'deg';
                [raw, old_flag, error_flag] = adis(X(:,ii));
                X_phys = polyval(fliplr(coeffs(1:2)),raw);        

            case 'inclt'    % ADIS inclinometer Temperature
                units = 'deg C';
                [raw, old_flag, error_flag] = adis(X(:,ii));
                X_phys = polyval(fliplr(coeffs(1:2)),raw);        
            

            case 'altimeter_error'    % Altimeter error signal
                units = ' ';
                X_phys = polyval(fliplr(coeffs(1:2)),X(:,ii));        
            
            case 'voltage'    % A generic linear voltage signal
                units = ' ';
                X_phys = polyval(fliplr(coeffs(1:2)),X(:,ii));

            case 'oxygen'      % Sea-Bird oxygen sensor
%                %%% Note that this version does not correct for thermal
%                %%% expansion/compressibility and for oxugen saturatipon;
%                %%% this must be done separately. It only outputs the frequency
                units = 'Hz';
                f = coeffs(7)*coeffs(6)./X(:,ii);
                X_phys = f;

            otherwise
                error(['Unknown sensor type: ' sensor2{ii}]);
                
        end

        outparams{ii}=X_phys; 
        units_all{ii}=units;
        
    end

elseif isXMP
%   %% odas v6 XMP
    for ii=lower(sensor)
        switch ii{1}
          case {'pres'}
            coef.pres(1) = str2double(setupstr(cfg,xmp_ser_num,'P_coef0'))';
            coef.pres(2) = str2double(setupstr(cfg,xmp_ser_num,'P_coef1'))';
              
          case {'therm'}
            coef.therm(1) = str2double(setupstr(cfg,xmp_ser_num,'T_coef0'))';
            coef.therm(2) = str2double(setupstr(cfg,xmp_ser_num,'T_coef1'))';
          
          case {'shear1'}
            coef.shear1(1) = str2double(setupstr(cfg,xmp_ser_num,'sh1_sens'));
            coef.shear1(2) = str2double(setupstr(cfg,xmp_ser_num,'sh1_diff_gain'));
             
          case {'shear2'}
            coef.shear2(1) = str2double(setupstr(cfg,xmp_ser_num,'sh2_sens'));
            coef.shear2(2) = str2double(setupstr(cfg,xmp_ser_num,'sh2_diff_gain'));
            
          case {'pitch'}
            coef.pitch(1) = str2double(setupstr(cfg,xmp_ser_num,'pitch_coef0'))';
            coef.pitch(2) = str2double(setupstr(cfg,xmp_ser_num,'pitch_coef1'))';
          
          case {'v_bat', 'gnd'}
            coef.(ii{1}) = [1;0];
          
          otherwise
            warning('Invalid sensor type: %s\n', ii{1});
        end
    end

    % Figure out what types of sensor has been provided. Note that this is not
    % case-sensitive.
    for ii=1:length(sensor)
        if strcmpi(sensor{ii},'pres') || strcmpi(sensor,'P') || strcmpi(sensor,'P_dP'), sensor2{ii} = 'pres'; 
        elseif length(sensor{ii})>=3 && strcmpi(sensor{ii}(1:3),'gnd'),   sensor2{ii} = 'gnd'; 
        elseif length(sensor{ii})>=5 && strcmpi(sensor{ii}(1:5),'shear'), sensor2{ii} = 'shear'; 
        elseif length(sensor{ii})>=5 && strcmpi(sensor{ii}(1:5),'therm'), sensor2{ii} = 'therm'; 
        elseif strcmpi(sensor{ii},'pitch') || strcmpi(sensor{ii},'roll'), sensor2{ii} = 'tilt'; 
        elseif strcmpi(sensor{ii},'v_bat'),                               sensor2{ii} = 'voltage';
        else
            error(['Unknown sensor type: ' sensor{ii}]);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% For each variable requested, extract the appropriate %%%%%%
    %%%%%% coefficients and apply based on the sensor type.     %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii=1:length(sensor)

        % Extract the calibration coefficients from the structure created above
        % (case-insensitive)
        junk = fieldnames(coef);
        jj = find(strcmpi(fieldnames(coef),sensor{ii})==1);
        if ~isempty(jj), coeffs = coef.(junk{jj});
        else error(['Unknown sensor: ' sensor{ii}]);
        end
        ADC_FS = 4.096;
        ADC_Bits = 16;
        
        switch sensor2{ii}

            case 'pres'     % Pressure
                units = '[ dBar ]';
                % coeffs(1) is the keller "Zero" value in mV
                % coeffs(2) is the Keller "SENS" value in mv/bar
                % this assumes a gain of 6 in the electronics and an
                % excitation current of 0.8533 mA
                a0 = -11.719*coeffs(1) / coeffs(2);
                a1 = 0.12208 / coeffs(2);
                X_phys = polyval([a1 a0],X(:,ii));

            case 'gnd'      % Ground voltage
                units = '[ count ]';
                X_phys = X(:,ii);

            case 'tilt'    % Analog Devices ADXL103 tilt sensor
                units = '[ ^{\circ} ]';
                % coeffs(1) is the maximum reading 
                % coeffs(2) is the minimum reading
                a0 = (coeffs(1) + coeffs(2))/2; % offset in counts
                a1 = (coeffs(1) - coeffs(2))/2; % sensitivity to gravity               
                X_phys = (180/pi)*asin((X(:,ii)-a0)/a1);
            
           case 'therm'    % Thermistor
                units = '[ ^{\circ}C ]';
                Gain = 1; % bridge gain
                E_B = 4.096; % Bridge excitation voltage
                a = 0; % offset
                b = 1; % slope                
                Z = ((X(:,ii) - a)/b) * (ADC_FS / 2^ADC_Bits)*2/(Gain*E_B);
                n = find(Z >  0.6); if ~isempty(n), Z(n) =  0.6; end; % prevents taking log of negative numbers
                n = find(Z < -0.6); if ~isempty(n), Z(n) = -0.6; end                
                X_phys = (1-Z) ./ (1+Z); %this is the resistance ratio
                X_phys = 1/coeffs(1) + (1/coeffs(2))*log(X_phys);    X_phys = (1 ./X_phys) - 273.15; %temperature in deg C
            
           case 'shear'    % Shear probes. Output must still be divided by
               %             speed-squared to get it into physical units.
                units = '[ s^{-1} ]';
                % coeffs(1) is shear probes sensitivity 
                % coeffs(2) is the differentiator gain
                X_phys = (ADC_FS / 2^ADC_Bits) * X(:,ii); 
                X_phys = X_phys ./(2*sqrt(2)*coeffs(1)*coeffs(2));
            
            case 'voltage'    % Battery voltge in XMP
                units = '[ V ]';
                X_phys = ADC_FS + X(:,ii) * 2 * ADC_FS / 2^ADC_Bits;

            otherwise
                error(['Unknown sensor type: ' sensor2{ii}]);
                
        end

        outparams{ii}=X_phys;
        units_all{ii}=units;
        
    end

else
%   %% odas V6

    % Copy the parameters from the setup file into a "coef" structure to allow
    % for easier access.
    for individual = sensor
      [S,K,V] = setupstr(cfg, individual{1}, '', '');
      for i = 1:length(K)
        try, V{i} = eval(V{i}); catch, end
        try, coef.(individual{1}).(K{i}) = V{i}; catch, end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% For each variable requested, extract the appropriate %%%%%%
    %%%%%% coefficients and apply based on the sensor type.     %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii=1:length(sensor)

        % coefficient names correspond to the parameter names used in the setup file
        params = coef.(sensor{ii}); %= eval(['coef.' sensor{ii}]);

        odas_type = setupstr(cfg, sensor{ii}, 'type');
        switch char(odas_type{1})

            case 'sbt'      % Sea-Bird temperatures
                units ='deg C';
                
                f = params.coef6*params.coef5./X(:,ii);
                c = log(params.coef4./f);
                X_phys = 1./polyval([params.coef3 params.coef2 params.coef1 params.coef0],c) - 273.15;                

            case 'sbc'      % Sea-Bird conductivity
%               %%% Note that this version does not correct for thermal
%               %%% expansion/compressibility; this must be done separately.
                units = 'mS/cm';
                f = params.coef6*params.coef5./X(:,ii)/1000;
                X_phys = polyval([params.coef4 params.coef3 params.coef2 params.coef1 params.coef0],f);

            case 'o2_43f'      % Sea-Bird oxygen sensor
                units = '?';
                f = params.coef3*params.coef2./X(:,ii);
                X_phys = 100*params.coef0 + params.coef1*f;                

            case 'poly'     % generic polynomial
                units = ' ';
                X_phys = polyval([params.coef3 params.coef2 params.coef1 params.coef0],X(:,ii));

            case 'gnd'      % Ground voltage
                units = 'V'; %shouldn't this be in counts?
                X_phys = X(:,ii);

            case 'accel'    % Accelerometers
                units = 'm sec^-2';
                X_phys = g*(X(:,ii)-params.coef0)/params.coef1;

            case 'magn'    % Magnetometers
                units = 'micro-Tesla';
                X_phys =   (X(:,ii)-params.coef0)/params.coef1;

            case 'Alec_EMC'    % Alec EM Current Meter
                units = 'm/s';
                X_phys = polyval([params.coef1 params.coef0],X(:,ii));

            case 'altimeter'    % Altimeter signal
                units = 'm';
                X_phys = polyval([params.coef1 params.coef0],X(:,ii));

            case 'inclxy'    % ADIS inclinometer x or y incline
                units = 'deg';
                [raw, old_flag, error_flag] = adis(X(:,ii));
                X_phys = polyval([params.coef1 params.coef0],raw);        

            case 'inclt'    % ADIS inclinometer Temperature
                units = 'deg C';
                [raw, old_flag, error_flag] = adis(X(:,ii));
                X_phys = polyval([params.coef1 params.coef0],raw);        
            
            case 'altimeter_error'    % Altimeter error signal
                units = ' ';
                X_phys = polyval([params.coef1 params.coef0],X(:,ii));        
            
            case 'voltage'    % A generic linear voltage signal
                units = 'V';
%               X_phys = X(:,ii) *params.G*params.adc_fs/2^params.adc_bits;
                X_phys = X(:,ii) *(params.adc_fs/2^params.adc_bits) / params.g;
               
            case 'therm'     % thermistor
                units = 'deg C';
                Z = ((X(:,ii) - params.a)/params.b) * (params.adc_fs/2^params.adc_bits)*2/(params.g*params.e_b);
                n = find(Z >  0.6); if ~isempty(n), Z(n) =  0.6; end; % prevents taking log of negative numbers
                n = find(Z < -0.6); if ~isempty(n), Z(n) = -0.6; end
                X_phys = (1-Z) ./ (1+Z); %this is the resistance ratio
                Log_R = log(X_phys);
                if isfield(params, 'beta')
                    X_phys = 1/params.t_0 + (1/params.beta)*Log_R;
                elseif isfield(params, 'beta_1')
                    X_phys = 1/params.t_0 + (1/params.beta_1)*Log_R;
                else
                    warning('No beta or beta_1 parameter for this thermistor')
                end
                if isfield(params, 'beta_2')
                    X_phys = X_phys + (1/params.beta_2)*Log_R.^2;
                    if isfield(params, 'beta_3')
                        X_phys = X_phys + (1/params.beta_3)*Log_R.^3;
                    end
                end
                X_phys = 1 ./X_phys - 273.15; %temperature in deg C
                
            case 'ucond'     % micro conductivity
                units = 'mS/cm';
%               X_phys = (params.adc_fs / 2^params.adc_bits) * X(:,ii) + params.adc_fs/2; %turn into voltage at uC board
                X_phys = (params.adc_fs / 2^params.adc_bits) * X(:,ii) ; %turn into voltage at uC board
                X_phys = (X_phys - params.a) ./ params.b; %it is now in units of conductance [S]
                X_phys = X_phys ./ params.k; %it is now in units of conductivity [S/m]
                X_phys = 10*X_phys; %and now in [mS/cm]

            case 'shear'     
                units = 's^{-1}';
                X_phys = (params.adc_fs / 2^params.adc_bits) * X(:,ii); 
                X_phys = X_phys ./(2*sqrt(2)*params.diff_gain*params.sens);


            otherwise
                error(['Unknown sensor type: ' char(setupstr(cfg, sensor{ii}, 'type'))]);
                
        end

        outparams{ii}=X_phys;
        
        if ~isfield(params, 'units')
            % if the units parameter is not set in the sensor's section,
            % use the default one as specified in the above case statements
            params.units = units;
        end
        
        units_all{ii}=params.units;
        
    end %for end

end %if end
   
outparams{length(outparams)+1} = units_all;

varargout = outparams;


% =========================================================================
%> @brief Helper function to convert inclinometer data into meaningful raw counts
%>
%> [X_Out, old_flag, error_flag] = adis(X_In) 
%>
%> The ADIS16209 inclinometer will produce data that needs to be converted into 
%> raw counts prior to conversion into physical units. This Function will extract
%> the 14 least significant bits for X and Y inclination and the least significant 
%> 12 bits for the Temperature channel. The remaining upper bits are not part of
%> the raw data but do contain some useful information. The data are 16 bit 
%> values with bit 15 indicating new data and bit 14 indicating an error. If the 
%> sensor is working properly, it will have the MS-bit set to 1 and the second 
%> MS-bit set to zero. The remaining 14-bits are data in 2s-compliment for the X- and the
%> Y-inclination, but the temperature data is only 12-bits, unsigned.
%> In the bulk-reading of the data files, all data are treated as
%> 2s-compliment signed numbers, which means (among other things) that all data
%> with the MS-bit equal to one are considered negative numbers.
%>
%> @param X_In - a vector of data from the inclinometer.
%>
%> @retval X_Out - a vector of data in units of counts (not in physical units)
%>   with the upper 2 bits properly stripped,
%> @retval old_flag - a vector of indicies to the data elements that are not new.
%> @retval error_flag - a vector of indicies to the data elements that have
%>   their error-flag set.
%>
%> Normally, error_flag and old_flag will be empty vectors.
%>
%> <b>Version History:</b>
%> - 2010-06-07 (RGL) original version.
% =========================================================================
function [X_Out, old_flag, error_flag] = adis(X_In) 

error_flag = [];
old_flag = find(X_In >= 0); % MS-bit is not set so data is not new.

n = find(X_In < -2^14); % these data points have the MS-bit set and are new.
if ~isempty(n)
    X_In(n) = X_In(n) + 2^15; % shift to clear MS-bit
end

n = find(X_In >= 2^14); % Second MS-bit is set, which indicates an error. 
if ~isempty(n)
    error_flag = n;
    X_In(n) = X_In(n) - 2^14;
end

% The data are 2s-compliment for inclination. So we find the upper half of the 
% 14-bit range. This is the data with negative values. So we shift it down
% to get a proper 2s-compliment conversion. The temperature data is 12-bit
% only and will automatically not get shifted.
n = find(X_In >= 2^13);
if ~isempty(n)
    X_In(n) = X_In(n) - 2^14;
end
X_Out = X_In;

    


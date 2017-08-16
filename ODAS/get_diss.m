%% get_diss
% Calculate dissipation rates over an entire profile.
%%
% <latex>\index{Type B!get\_diss}</latex>
%
%%% Syntax
%   diss = get_diss( SH, A, info )
%
% * [SH] Shear probe signals in physical units. It is assumed to be high-pass 
%        filtered to remove low-frequency wobble. Cut-off frequency should be 
%        roughly twice the fft-length, equivalent. Values are formatted as a 
%        matrix where each column is a shear probe.
% * [A] Matrix of acceleration signals. Does not need to be in physical units.
% * [info] Structure containing additional elements.
% * []
% * [diss] Structure of results including the dissipation rate of TKE and 
%          related information. Results formatted as matrices with a row for 
%          each dissipation estimate.
%
%%% Description
%
% Calculate dissipation rates over an entire profile. 
%
%%% Input Structure
%
% * [info.fft_length] length of each fft segment in units of samples. For example
%      info.fft_length = 1024
% * [info.diss_length] the length of data used for each dissipation estimate.
%      diss_length should be at least 2 times longer than fft_length. I
%      recommend no less than a factor of 3. 
% * [info.overlap] the number of samples by which this function moves forward to
%      make the next dissipation estimate. Should be larger than or equal to
%      fft_length, and smaller than or equal to diss_length. I recommend
%       diss_length/2. If you do not want any overlap, set overlap =
%       diss_length.
% * [info.fs] sampling rate in Hz.
% * [info.speed] a scalar or vector of profiling speed to derive wavenumber spectra.
%      If a vector, its length must equal the number of rows in SH.
% * [info.T] a vector of temperature in degree C for calculating kinematic
%      viscosity. Length must match SH. Salinity = 35 is assumed.
% * [info.P] a vector of pressure for plotting the dissipation rates. Must match SH.
% * [info.Tolerance] optional logarithmic tolerance condition on the average ratio
%      of the observed spectrum divided by the Nasmyth spectrum. 
%      The ratio is calculated from the lowest wavenumber to K = 0.04*k_s
%      (a position just past the peak of the spectrum), where K is the 
%      wavenumber in cpm (cyles per metre) and k_s is the Kolmogorov wavenumber,
%      k_s = (epsilon/nu^3)^{1/4} in units of radians per metre. 
%      The default value is 0.3, which corresponds to a factor of 2 = (10^{0.3}).
% * [info.fit_order] the order of the polynomial fit to the spectrum in log-log
%       space used to estimate the spectral minimum wavenumber and this wavenumber
%       is the inital estimate of the upper limit of integration for the estimation
%       of the rate of dissipation of TKE. Optional - default = 6.
% * [info.f_AA] the cut-off frequency of the anti-aliasing filter. Default = 98 Hz.
% * [info.fit_2_Nasmyth] boolean value indicating a forced fit to the Nasmyth
%       spectrum in the +1/3-range is desired.  Default value is 0.
%
%%% Output Structure
%
% A structure of dissipation rate of TKE and related information. The
% number of rows equals the number of dissipation estimates. The elements
% are:
%
% * [diss.e] the rate of dissipation of TKE [W/kg]. One column for every shear
%      probe.
% * [diss.K] wavenumbers [cpm, or cycles per metre] with one column for every dissipation estimate.
% * [diss.sh_clean] the wavenumber spectrum for each shear probe signal at each
%      dissipation estimate. Every row is a 3-D matrix, where the highest
%      dimension is the wavenumber index. The second and third dimensions are the
%      cross-spectral matrix of the shear probes. The diagonal elements are the
%      auto-spectra.
% * [diss.sh] same as sh_clean but without cleaning by the Goodman coherent noise
%       removal algorithm.
% * [diss.AA] the cross-spectral matrix of acceleration signals. The diagonal
%      elements are the auto-spectra. 
% * [diss.UA] the cross-spectral matrix of shear probe and acceleration signals.
% * [diss.F] frequencies, with one column for every dissipation estimate.
%      Currently every column is identical.
% * [diss.speed] a column vector of the mean speed at each dissipation estimate.
% * [diss.nu] column vector of kinematic viscosity [m^2/s] at each dissipation estimate.
% * [diss.P] column vector of pressure at each dissipation estimate.  Derived directly
%      from the input.
% * [diss.T] column vector of temperature at each dissipation estimate.  Derived directly
%      from the input.
% * [diss.Nasmyth] the Nasmyth spectrum for the estimated dissipation rate,
%      evaluated at the same wavenumbers as K

% *Version History:*
% 2012-05-27 (RGL) started at Narita airport.
% 2012-06-26 (RGL) added while loop to adjust upper limit of integration if
%            spectra agree poorly with Nasmyth for wavenumber just beyond the
%            spectral peak.
% 2012-07-11 (RGL and WID) changed detection of spectral minimum. Added
%            mordern parser.
% 2012-10-30 (RGL) Temporary modification to make K_max go to at least 100 cpm
% 2012-10-31 (RGL) Added feature to force a fit to the Nasmyth spectrum in the 
%            +1/3 slope range, if requested.
% 2012-11-08 (WID) Replaced call to ismatrix() with the logic equivalent - to 
%            allow the function to run on earlier versions of Matlab.
% 2012-11-09 (WID) Documentation update

function diss = get_diss(SH, A, varargin)

% Default values for optional fields
default_Tolerance = 0.3;
default_fit_order = 6;
default_f_AA      = 98;
default_fit_2_Nasmyth = 0;

p = inputParser;

val_Tolerance  = @(x) isnumeric(x) && (x >= 0.1) && (x <= 1);
val_fft_length = @(x) isnumeric(x) && isscalar(x) && (x >= 2);
val_positive   = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
val_matrix     = @(x) isnumeric(x) && (size(x,1) > 1) && (size(x,2) >= 1);
val_speed      = @(x) isnumeric(x) && isvector(x) && (x >= 0);

addRequired(  p, 'SH',          val_matrix);
addRequired(  p, 'A',           val_matrix);
addParamValue(p, 'fft_length',  val_fft_length);
addParamValue(p, 'diss_length', val_fft_length);
addParamValue(p, 'overlap',     val_fft_length);
addParamValue(p, 'fs',          val_positive);
addParamValue(p, 'speed',       val_speed);
addParamValue(p, 'T',           val_matrix);
addParamValue(p, 'P',           val_matrix);
addParamValue(p, 'Tolerance',   default_Tolerance, val_Tolerance);
addParamValue(p, 'fit_order',   default_fit_order, val_positive);
addParamValue(p, 'f_AA',        default_f_AA,      val_positive);
addParamValue(p, 'fit_2_Nasmyth', default_fit_2_Nasmyth, val_positive);

% Parse the arguments.
parse(p, SH, A, varargin{:});

% Perform last stages of input validation.
if p.Results.diss_length < 2*p.Results.fft_length,
  error('Invalid size for diss_length - must be greater than 2 * fft_length.');
end

if (size(p.Results.SH,1) ~= size(p.Results.A,1) || ...
    size(p.Results.SH,1) ~= size(p.Results.T,1) || ...
    size(p.Results.SH,1) ~= size(p.Results.P,1)),
  error('Same number of rows required for SH, A, T, P.');
end

if ~isscalar(p.Results.speed) && size(p.Results.SH,1) ~= size(p.Results.speed,1),
  error ('speed vector must have the same number of rows as shear');
end

if p.Results.diss_length > size(p.Results.SH,1)
  error('Diss_length cannot be longer than the length of the shear vectors');
end

fft_length  = p.Results.fft_length;
diss_length = p.Results.diss_length;
overlap     = p.Results.overlap;
fs          = p.Results.fs;
speed       = p.Results.speed;
T           = p.Results.T;
P           = p.Results.P;
Tolerance   = p.Results.Tolerance;
fit_order   = p.Results.fit_order;
f_AA        = p.Results.f_AA;
fit_2_Nasmyth = p.Results.fit_2_Nasmyth;

% end of input argument checking.
select = (1:diss_length)';

% Calculate the number of dissipation estimates that can be made.
number_of_rows = 1 + floor((size(SH,1) - diss_length) / overlap);
F_length = 1 + floor(fft_length/2); % size of frequency and wavenumber vectors

% pre-allocate matrices
diss.e            = zeros(number_of_rows,size(SH,2));
diss.K_max        = diss.e;
                                          % third dimension is frequency index
diss.Nasmyth_spec = zeros(number_of_rows, size(SH,2), F_length);
                                          % fourth dimension is frequency index
diss.sh           = zeros(number_of_rows, size(SH,2), size(SH,2), F_length);
diss.sh_clean     = diss.sh;
diss.AA           = zeros(number_of_rows, size(A,2),  size(A,2),  F_length);
diss.UA           = zeros(number_of_rows, size(SH,2), size(A,2),  F_length);
diss.F            = zeros(number_of_rows,                         F_length); 
diss.K            = diss.F;
diss.speed        = zeros(number_of_rows,1);
diss.nu           = diss.speed;
diss.P            = diss.speed;
diss.T            = diss.speed;

Nasmyth_spectrum  = zeros(size(SH,2),F_length);
index = 1;
while select(end) <= size(SH,1)
    [P_sh_clean, AA, P_sh, UA, F] = clean_shear_spec(A(select,:), SH(select,:), fft_length, fs);
% Note, P_sh_clean and P_sh are [M M N] matrices where M is the number of
% columns in SH (the number of shear probe signals) and N is the length of
% the frequency vector F (usually fft_length/2 + 1).
% If there is only a single shear probe, these matrices are [N 1] in size,
% and I use the "squeeze" function to remove the redundant dimensions.
    
    % Convert frequency spectra to wavenumber spectra
    W = mean(abs(speed(select)));
    
    K = F/W;
    
    junk = repmat(K,[1 size(P_sh,2),size(P_sh,2)]);
    junk = permute(junk,[2 3 1]); % move wavenumbers into third dimension
    junk = squeeze(junk); % remove extra dimensions if SH is only a vector
    
    P_sh_clean = P_sh_clean * W; P_sh = P_sh * W;
    P_sh_clean = P_sh_clean .* (1 + (junk / 48).^2); % Wavenumber correction after Macoun & Lueck
    P_sh       =       P_sh .* (1 + (junk / 48).^2);
    
    % Now we find a good value for the upper limit of spectral integration
    % for the computation of shear variance and the rate of dissipation of
    % TKE.
    % Find first minimum of spectrum in log-log space. Look only at
    % wavenumbers smaller than 150 cpm. This is a little complicated.
    e = zeros(1,size(SH,2)); % pre-assign, one column per shear probe
    K_max = e;

    mean_T = mean(T(select));
    nu     = visc35(mean_T);
    
    for column_index = 1:size(SH,2)
        if size(SH,2) == 1 % we have only a sinle shear probe
            shear_spectrum = P_sh_clean;
        else % the auto-spectra are on the diagnonal
            shear_spectrum = P_sh_clean(column_index, column_index, :);
        end
        shear_spectrum = squeeze(shear_spectrum);
        % Do not make a fit beyond 150 cpm or f_AA/W
        valid_shear = find( K <= fit_in_range(f_AA/W, [0 150]) );
        Index_limit = length(valid_shear);
        y = log10(shear_spectrum(2:Index_limit));
        x = log10(K(2:Index_limit));
        % Note K(1) is always zero and has problems with logarithms. Now we fit a
        % polyninomial to the spectrum
                
        % Keep the fit_order reasonable. Any value from 4 to 8 could work
        fit_order = fit_in_range(fit_order, [3 8]);
        
        % p   - polynomials of a polynomial fit of data
        % pd1 - first degree derivative of polynomials
        % pr1 - roots of first degree derivative polynomials
        p = polyfit(x, y, fit_order);
        pd1 = polyder(p);
        pr1 = sort( roots(pd1) );
        pr1 = pr1(~imag(pr1)); % remove complex roots
        % Filter roots so that only the minima above 10cpm remain.
        if ~isempty(pr1)
            pr1 = pr1( polyval( polyder(pd1), pr1 ) > 0 );      % minima only
            pr1 = pr1( pr1 > log10(10) );                       % 10cpm
            % Fit root within a given range.
        end
        K_limit = fit_in_range( pr1, [log10(10) log10(150)], log10(150) );
        Range = find( K <= 10^(K_limit) + 0.5 ); %Integration range for shear probes.
        
        % Now check that the integration limits are reasonable by comparing the
        % Nasyth spectrum against the measured spectrum for non-dimensional
        % wavenumbers up to 0.04 (when the wavenumber is in units of cpm). This
        % point is just above the spectral maximum. If the measured spectrum is
        % consistently above the nasmyth spectrum, then increase the upper limit of
        % integration. Else do the opposite. However, the range must never be
        % less than 10 cpm and never more than 150 cpm.
        
        
        % Here, we adjust the limit of integration in order to get a better fit
        % to the Nasmyth spectrum
        mean_error = 1; % just a starting value to force the next while loop
        loop_count = 0;
        while abs(mean_error) > Tolerance,
            if loop_count > 10, break, end
            loop_count = loop_count + 1;
            e(column_index) = 7.5*nu*trapz(K(Range), shear_spectrum(Range));
            k_s = sqrt(sqrt((e(column_index)/nu^3))); % the Komogorov wavenumber
            check_range = find(K <= 0.04*k_s);
            Nasmyth_values = nasmyth(e(column_index),nu,check_range);
            mean_error = mean(log10(shear_spectrum(check_range) ./ Nasmyth_values));
            if mean_error < -Tolerance
                Range = Range(1:round(0.75*Range(end))); %shrink the range to 3/4
                if K(Range(end)) < 10 % cannot shrink below 10 cpm. Use 10 cpm and give up.
                    Range = find(K <= 10);
                    mean_error = 0; % That is all we can do. Exit this loop.
                end
            elseif mean_error > Tolerance
                upper_limit = 1.33*K(Range(end)); % try to expand the range to 4/3
                if upper_limit > 150,    upper_limit = 150;   mean_error = 0;end
                if upper_limit > K(end), upper_limit = K(end);mean_error = 0;end
                Range = find(K <= upper_limit); %expand the range to 4/3
            end
            if mean_error == 0, % that is all that we can do. Exit this loop.
                e(column_index) = 7.5*nu*trapz(K(Range), shear_spectrum(Range));
            end
            % Now we force a fit to the Nasmyth spectrum, in the +1/3-range, if requested
            if fit_2_Nasmyth
                for Nasmyth_loop_count = 1:2
                    k_s = sqrt(sqrt((e(column_index)/nu^3))); % the Komogorov wavenumber
                    NF_Range = find(K <= 0.75*f_AA / W); % fit only up to 75% of F_AA
                    if (K(NF_Range(end)) > 0.015*k_s) % Do not go beyond the +1/3 slope range
                        NF_Range = find(K <= 0.015*k_s); % fit only up to 75% of F_AA
                    end
                    Nasmyth_values = nasmyth(e(column_index),nu,NF_Range);
                    N_fit_error = mean(log10(shear_spectrum(NF_Range) ./ Nasmyth_values));
                    if N_fit_error < 0, N_fit_error = 0; end % Measured spectrum is below the Nasmyth spectrum
                    e(column_index) = e(column_index) * 10^(1.5*N_fit_error);% scale up epsilon
                    % Now do it all again, one time only, in case the change was large
                end
            end
                 
            % Now a small correction for the lowest wavenumbers, K(1)=0 and
            % K(2)= lowest non-zero wavenumber. Trapazoidal integration equals 
            % only 2/3 of the correct integral.
            junk = e(column_index);
            e(column_index) = junk*(1 + 15*nu*junk^(2/3)*K(2)^(4/3));
            K_max(column_index) = K(Range(end));
            Nasmyth_spectrum(column_index,:) = nasmyth(e(column_index),nu,K);
        end
    end
    
    diss.e(index,:) = e;
    diss.Nasmyth_spec(index,:) = Nasmyth_spectrum(:);
    diss.K_max(index,:) = K_max;
    diss.sh_clean(index,:) = P_sh_clean(:);
    diss.sh(index,:) = P_sh(:);
    diss.AA(index,:) = AA(:);
    diss.UA(index,:) = UA(:);
    diss.F(index,:) = F;
    diss.K(index,:) = K;
    diss.speed(index,:) = W;
    diss.nu(index,:) = nu;
    diss.T(index,:) = mean_T;
    diss.P(index,:) = mean(P(select));
    
    index = index + 1;

    select = select + overlap;
end

end

function result = fit_in_range(value, range, default)
% Fit first element of value into range.  Returns scalar result.
% - range is of the form [min max]
% - if value is empty, use default value

if ~isempty(value), 
  result = value(1);
else
  result = default(1);
end

if result < range(1), result = range(1); end
if result > range(2), result = range(2); end

end
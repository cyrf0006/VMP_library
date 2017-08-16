%% clean_shear_spec
% Remove acceleration contamination from all shear channels.
%%
% <latex>\index{Type A!clean\_shear\_spec}</latex>
%
%%% Syntax
%   [cleanUU, AA, UU, UA, F] = clean_shear_spec( A, U, nFFT, rate )
%
% * [A] Matrix of acceleration signals usually represented by [Ax Ay Az] 
%       where Ax Ay Az are column vectors of acceleration components.
% * [U] Matrix of shear probe signals with each column representing a 
%       single shear probe.
% * [nFFT] Length of the fft used for the calculation of auto- and cross-spectra.
% * [rate] Sampling rate of the data in Hz
% * []
% * [cleanUU] Matrix of shear probe cross-spectra after the coherent 
%             acceleration signals have been removed. The diagonal has the 
%             auto spectra of the clean shear.
% * [AA] Matrix of acceleration cross-spectra. The diagonal has the auto-spectra.
% * [UU] Matrix of shear probe cross-spectra without the removal of coherent
%        acceleration signals. The diagonal has the auto spectra of the original 
%        shear signal.
% * [UA] Matrix of cross-spectra of shear and acceleration signals.
% * [F] Column vector of frequency for the auto- and cross-spectra.
%
%%% Description
% Remove components within a shear probe signal that are coherent with 
% signals reported by the accelerometers. It is very effective at removing
% vibrational contamination from shear probe signals.
% 
% It works only in the spectral domain. 
%
% Developed by Louis Goodman (University of Massachusetts) and adapted by RSI.
%
% For best results and statistical significance, the length of the fft (nFFT) 
% should be several times shorter than the length of the vectors, [size(U,1)], 
% passed to this function. That is, several ffts have to be used to form an 
% ensemble-averaged auto- and cross-spectrum.
% 
% @image @images/clean_shear_spectra.pdf @Coherent noise removal applied to shear
% probe spectrum. @Original spectrum (thin green line) and after noise removal 
% (thick green line). Data collected with Remus-600 in Buzzards Bay MA, courtesy 
% of Tim Boyd , Scottish Association for Marine Science. Only shear probe 2 was 
% installed, so, the spectrum for probe 1 falls off the deep end of the plot. 
%

% Version History
%
% * circa 2003 (Lou Goodman)
% * 2004-05-03 (RGL) tweaks
% * 2006-05-05 (RGL) "squeezed" some of the output vectors
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-15 (WID) updated documentation for publishing
% * 2012-07-11 (ILG) replaced csd_rolf and psd_rolf function calls with calls
%                    to csd_odas.m; for loops take advantage of additional 
%                    outputs available from csd_odas.m

function  [clean_UU, AA, UU, UA, F] = clean_shear_spec(A, U, n_fft, rate)

% do some checking
if ((length(A(:)) == 1) | (length(U(:)) == 1))
    error('Acceleration and shear matrices must contain vectors')
end
if (size(A,2) > size(A,1))
    error('the vectors do not seem to go down the columns of the acceleration matrix')
end
if (size(U,2) > size(U,1))
    error('The vectors do not seem to go down the columns of the shear matrix')
end
if (size(A,1)~=size(U,1))
    error('Acceleration and shear matrices must have the same number of rows')
end
if (size(A,1) <= 2*n_fft)
    error('Your fft length is too long for the length of the vectors')
end
if (rate <= 0)
    error('Sampling rate is negative')
end
% end of checking

% Pre-allocate shear matrices
AA = zeros(size(A,2),size(A,2),n_fft/2+1);
UU = zeros(size(U,2),size(U,2),n_fft/2+1);
UA = zeros(size(U,2),size(A,2),n_fft/2+1);
clean_UU = zeros(size(U,2),size(U,2),n_fft/2+1);

% Compute auto- and cross-spectra
for k = 1:size(U,2)
    % Shear-probe and acceleration cross-spectra, auto-spectra
    for m = 1:size(A,2)
        [UA(k,m,:), F, UU(k,k,:), AA(m,m,:)] = csd_odas(U(:,k),...
            A(:,m),n_fft,rate,[],n_fft/2,'parabolic');
    end
    % Shear probe cross-spectra
    for n = k+1:size(U,2)
        UU(k,n,:) = csd_odas(U(:,k),U(:,n),n_fft,rate,[],n_fft/2,'parabolic');
        UU(n,k,:) = conj(UU(k,n,:));
    end
end
% Acceleration cross-spectra
for m = 1:size(A,2)
    for n = m+1:size(A,2)
        AA(m,n,:) = csd_odas(A(:,m), A(:,n), n_fft,rate,[],n_fft/2,'parabolic');
        AA(n,m,:) = conj(AA(m,n,:));
    end
end

% Clean the shear spectra
for ii = 1:length(F)
   clean_UU (:,:,ii) = UU(:,:,ii) - (UA(:,:,ii)/AA(:,:,ii)) *conj(UA(:,:,ii)).';
end

% Remove extra dimensions
clean_UU = squeeze(clean_UU);clean_UU = real(clean_UU);
UU = squeeze(UU);
AA = squeeze(AA); 
UA = squeeze(UA);%UA = UA';

   

function  [clean_UU, F, AA, UU, UA] = clean_shear_spec_fc(A, U, rate)
%
%function [clean_UU, AA, UU, UA, F] = clean_shear_spec(A, U, rate)
%
% function to remove acceleration contamination from all shear channels
% using all available acceleration channels and other signals related shear contamination.
% Uses full and partial coherency. 
% A = [ax ay az ...] are the accelerations (can be more than 3),
% U = [sv sw ...] are the shear probe signals (can be more than 2).
% The time series for each vector must progress down the columns of the
% matrices.
% n_fft length of fft for calculation of auto- and cross-spectra
% rate is the sampling rate of the time series.
%
% The outputs are:
% the cleaned up cross-spectrum of shear, clean_UU,
% the acceleration cross-spectra, AA
% the original (dirty) shear cross-spectrum, UU,
% the cross-spectrum of dirty shear and acceleration, UA, and
% the frequency of the cross-spectra, F
% If you only need the cleaned-up cross-spectra of shear use
% clean_UU = clean_shear_spec(A,U,n_fft,rate);
%
% Written by Lou Goodman, circa 2003.
% Dicked with by R. Lueck 2004-05-03.
% "squeezed" some of the output vectors. R. Lueck 2006-05-05

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
% $$$ if (size(A,1) <= 2*n_fft)
% $$$     error('Your fft length is too long for the length of the vectors')
% $$$ end
if (rate <= 0)
    error('Sampling rate is negative')
end
% end of checking

% Now the acceleration cross-spectra
%AA = zeros(size(A,2),size(A,2),n_fft/2 + 1);% pre-allocated the
%keyboard                                            % matrix
for k = 1:size(A,2)
    %[AA(k,k,:), F] = psd_rolf(A(:,k),       n_fft,rate); % -FC
    [AA(k,k,:), F] = pwelch(A(:,k), [], [], [],rate);
    for m= k+1:size(A,2)
        %AA(k,m,:) = csd_rolf(A(:,k),A(:,m),n_fft,rate); % -FC
        AA(k,m,:) = cpsd(A(:,k), A(:,m), [], [], [], rate);
        AA(m,k,:) = conj(AA(k,m,:));
    end
end

% Now the shear probe cross-spectra
%UU = zeros(size(U,2),size(U,2),n_fft/2 + 1);
for k = 1:size(U,2)
    %UU(k,k,:) = psd_rolf(U(:,k),n_fft,rate); % -FC
    UU(k,k,:) = pwelch(U(:,k), [], [], [], rate);
    for m = k+1:size(U,2)
        %UU(k,m,:) = csd_rolf(U(:,k),U(:,m),n_fft,rate);
        UU(k,m,:) = cpsd(U(:,k), U(:,m), [], [], [], rate);
        UU(m,k,:) = conj(UU(k,m,:));        
    end
end

% Now the shear probe and acceleration cross-spectra
%UA = zeros(size(U,2),size(A,2),n_fft/2 + 1);
for k = 1:size(U,2)
    for m = 1:size(A,2)
        % UA(k,m,:) = csd_rolf(U(:,k),A(:,m),n_fft,rate);
        UA(k,m,:) = cpsd(U(:,k), A(:,m), [], [], [], rate);
    end
end

for ii = 1:length(F)
   clean_UU (:,:,ii) = UU(:,:,ii) - (UA(:,:,ii)/AA(:,:,ii)) *conj(UA(:,:,ii)).';
end
clean_UU = squeeze(clean_UU);
clean_UU = real(clean_UU);
UU = squeeze(UU);
AA = squeeze(AA); 
UA = squeeze(UA);
UA = UA';


   
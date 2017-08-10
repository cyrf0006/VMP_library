function [chi,zfft] = quickchi_vmp(dtdz, A, z, W, visc, shearFreqHz, pressFreqHz, zbin)

% function [chi, zfft] = quickchi_vmp(dtdz, A, z, W, visc, shearFreqHz, pressFreqHz, zbin)
% 
% usage ex: [chi1, p_chi1] = chi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1);
% called in chi_profile.m
%
% Pretty similar to chi_vmp.m, but only compute chi, no fit to
% Batchelor is done! Only for a quick check of chi properties.
%
% FCyr - March 2012
%--------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------- PART 1: INPUT & PARAM SETUP ----------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ % ---- Handle input arguments ---- %
% $$$ allowNewFields = 0;
% $$$ isCaseSensitive = 0;
% $$$ defaultVals.zbin = 1;
% $$$ defaultVals.iplt = 0;
% $$$ defaultVals.pfilename = [];
% $$$ defaultVals.explicit = 0; % explicit on plotting -FC 
% $$$ 
% $$$ parValStruct = parvalpairs(defaultVals,varargin{:},[allowNewFields isCaseSensitive]);
% $$$ zbin = parValStruct.zbin;
% $$$ % --------------------------------- %                          

% zbin makes BIG changes!
% zbin = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------- PART 2: data preparation (bins) --------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- Compute good indexes ---- %
% interp from slow freq. (64Hz) to fast freq. (512Hz) (zs = zshear)
zs=interp1([1:length(z)],z,[1:length(dtdz)]'*(pressFreqHz/shearFreqHz));
if isempty(zs)==1 % no shear data to process
   epsilon = [];
   zfft = [];
   return;
end

z1 = ceil(min(abs(z)));  % top of profile
z2 = floor(max(abs(z)));  % bottom of profile
zfft = z1+zbin/2:zbin:z2-zbin/2; %reg. z prof. for FFT

if isempty(zfft), % no shear data to process
   epsilon = [];
   zfft = [];
   return;
end % if

I=find(zs>=z1 & zs<=z2); % we keep only full bins
dtdz=dtdz(I);
zs = zs(I); 
A = A(I,:); % if shear is adjusted, size of A must be too
% ----------------------------- %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------ PART 3: Spectral integration ----------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---- MAIN LOOP! ---- %%
for i = 1:length(zfft)

    % ----- compute indexes on which FFT will be applied ----- % 
    I = find(zs>zfft(i)-zbin/2 & zs<zfft(i)+zbin/2); % I = Imicrosctruc
    Ifine = find(z>zfft(i)-zbin/2 & z<zfft(i)+zbin/2); % Ifinesctruc

    if (isempty(I)==1)
        chi(i) = NaN;
        continue
    end
    
    II = find(isnan(dtdz(I))==1);
   
    if length(II)>round(length(I)/4)
        disp(sprintf('skip bin %d', zfft(i)))
        chi(i) = NaN;
        continue % skip if more than 25% NaNs in the bin
    else
        dtdz(I(II))=0; % pad NaNs with zeros
    end
    % -------------------------------------------------------- %
    
    % --- Compute FFT & clean spectrum ---- %
    [ps0, f0] = pwelch(dtdz(I), [], [], [], 512); % Raw Power Spectrum
    J = find(f0 ~= 0); % remove zeros
    ps0 = ps0(J);
    f0 = f0(J);
    % ------------------------------------- %

    % ---- Conversion Freq -> Wavenumber ---- %
    Wm = mean(W(Ifine)); % mean falling speed for this bin 
    k0 = f0/Wm;
    %---------------------------------------- %
    
    % ----------- chi calculation ------------------- %
    nois=noise(k0); %get the noise at each value of k      ---> FC: noise must be adjusted
    dk=k0(3)-k0(2);
    Dt=.00000014; % Heat molec. diffus.
    chi_i = 6*Dt*sum((ps0'-nois).*dk); % (Eqn. 9)
    if chi_i > 0
        chi(i) = chi_i;
    else
        chi(i) = NaN;
    end
    
    % ------------------------------------------------ %
    
end

%% ---- END MAIN LOOP! ---- %%
   
% ------------------------ End of function ------------------------ %



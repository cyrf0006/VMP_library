function [chi,epsilon,zfft, ratio] = chi_vmp(dtdz, A, z, W, visc, shearFreqHz, pressFreqHz, zbin)

% function [chi,epsilon,zfft, ratio] = chi_vmp(dtdz, A, z, W, visc, shearFreqHz, pressFreqHz, zbin)
% 
% usage ex: [chi1, eps1, p_chi1, ratio] = chi_vmp(dt1dz,[ax ay az],P, W, nu, fs, FS, 1);
% called in chi_profile.m
%
% Compute chi, epsilon, ratio, etc. by fitting the batchelor
% spectrum using MLE technique (Ruddick et al. 2000)
%
% FCyr - Feb./March 2012, modified from B. Ruddick
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

z1 = ceil(min(abs(z)));  % top of profile <-------might want to adjust this (round-up bin limits)
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
        disp('skip bin')
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
    %    [ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A(I,:), dtdz(I), 512); % Clean Power Spectrum 
    % ------------------------------------- %

    % ---- Conversion Freq -> Wavenumber ---- %
    Wm = mean(W(Ifine)); % mean falling speed for this bin 
    
    k0 = f0/Wm;
    %k1 = f1/Wm;
   
    
% $$$     % corrected spectrum according to Macoun & Lueck 04
% $$$     ps0 = respcorr(ps0, k0);
% $$$     ps1 = respcorr(ps1, k1);
    % --------------------------------------- %

% $$$     disp('chi_vmp')
% $$$     keyboard
    
    % ---- Dissipation rate calculation ---- %
    visco = mean(visc(Ifine));
    %epsilon(i) = DAN'S script!!!              <-------------------------------- HERE DANIEL!!!!    
    %epsilon(i) = dewey_opti(ps1, k1, f1, visco);    

% $$$     figure(1)
% $$$     clf
    [chi_i, eps_i, likelihoodratio,likelihood,kb,f11] =  fit_kb(k0,ps0, 0);
% $$$     loglog(k0, ps0)
% $$$     hold on
% $$$     loglog(k0, f11, 'r')
% $$$     title(sprintf('z = %0.5g \n K_b = %0.5g', zfft(i), kb));
 
    %    keyboard
    chi(i) = chi_i;
    epsilon(i) = eps_i; 
    ratio(i) = likelihoodratio;
    
    
end


%% ---- END MAIN LOOP! ---- %%
   
% ------------------------ End of function ------------------------ %



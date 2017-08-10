function [epsilon,zfft] = epsilon_vmp_fc(shear, A, z, W, visc, shearFreqHz, pressFreqHz, varargin)
% function [epsilon,zfft] = epsilon_vmp_fc(shear, A, z, W, visc, shearFreqHz, pressFreqHz, varargin)
% 
% 
%
%
%--------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------- PART 1: INPUT & PARAM SETUP ----------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Handle input arguments ---- %
allowNewFields = 0;
isCaseSensitive = 0;
defaultVals.zbin = 4;
defaultVals.iplt = 0;
defaultVals.pfilename = [];
defaultVals.explicit = 0; % explicit on plotting -FC 

parValStruct = parvalpairs(defaultVals,varargin{:},[allowNewFields isCaseSensitive]);
zbin = parValStruct.zbin;
% --------------------------------- %                          


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------- PART 2: data preparation (bins) --------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- Compute good indexes ---- %
% interp from slow freq. (64Hz) to fast freq. (512Hz) (zs = zshear)
zs=interp1([1:length(z)],z,[1:length(shear)]'*(pressFreqHz/shearFreqHz));
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
shear=shear(I);
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
        epsilon(i) = NaN;
        continue
    end
    
    II = find(isnan(shear(I))==1);
   
    if length(II)>round(length(I)/4)
        disp('skip bin')
        epsilon(i) = NaN;
        continue % skip if more than 25% NaNs in the bin
    else
        shear(I(II))=0; % pad NaNs with zeros
    end
    % -------------------------------------------------------- %
    
    % SHOULD DETREND HERE!!
    
    % --- Compute FFT & clean spectrum ---- %
    [ps0, f0] = pwelch(shear(I), [], [], [], 512); % Raw Power Spectrum
    [ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A(I,:), shear(I), 512); % Clean Power Spectrum 
    % ------------------------------------- %

    % ---- Conversion Freq -> Wavenumber ---- %
    Wm = mean(W(Ifine)); % mean falling speed for this bin 
    
    k0 = f0/Wm;
    k1 = f1/Wm;
    
    % corrected spectrum according to Macoun & Lueck 04
    ps0 = respcorr(ps0, k0);
    ps1 = respcorr(ps1, k1);
    % --------------------------------------- %

    
    % ---- Dissipation rate calculation ---- %
    visco = mean(visc(Ifine));
    %epsilon(i) = DAN'S script!!!              <-------------------------------- HERE DANIEL!!!!
    epsilon(i) = dewey_opti(ps1, k1, f1, visco);  
    
    % Instrument resolution ~1e-10
    if epsilon(i)<1e-10;
        epsilon(i) = NaN;
    end
    % -------------------------------------- %
    
    % Test eroneous epsilon
    if epsilon(i) > 1e-5
        figure(1)
        clf
        
        subplot(121)
        loglog(k0, ps0)
        hold on
        loglog(k1, ps1, 'r')        
        title(sprintf('z = %0.1f m   e = %d', zfft(i), epsilon(i)))
        xlim([1e0 1e3])
        [phi,k] = nasmyth(epsilon(i), visco, length(ps1));
        loglog(k, phi, 'k')
        hold off
        %axis([1e0 1e3 1e-8 1e-2])

        subplot(122)
        ind1 = I(1)-4*length(I);
        ind2 = I(end)+4*length(I);

        if ind1<1
            ind1 = 1;
        elseif ind2 > length(zs)
            ind2 = length(zs);
        end
        plot(shear(ind1:ind2), zs(ind1:ind2))
        hold on
        plot(shear(I), zs(I), 'r')
        hold off
        set(gca, 'ydir', 'reverse')
        
        disp('Epsilon is very high. Do you want a check this?')
        answer=0;
        while answer==0
            R1 = input('Do you accept the profile? (y/n) ', 's');
            if strcmp('y',R1)==1 
                answer=1;
            elseif strcmp('n',R1)==1
                epsilon(i) = NaN;
                answer=1;
            else
                disp('bad answer! please type y or n')
            end
        end
    end % ----------  end test ----- %

           
    
    
end
%% ---- END MAIN LOOP! ---- %%
   
% ------------------------ End of function ------------------------ %



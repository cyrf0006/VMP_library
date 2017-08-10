function epsilon = patch_epsilon(shear, A, Wm, visco)

% function epsilon = patch_epsilon(shear, A)
%   
% Calculate epsilon on patch of variable size
%
    
% remove NaNs
II = find(isnan(shear)==1);
if length(II)>round(length(shear)/4)
    disp('skip bin')
    epsilon = NaN;
    return % skip if more than 25% NaNs in the bin
else
    shear(II)=0; % pad NaNs with zeros
end



% --- Compute FFT & clean spectrum ---- %
[ps0, f0] = pwelch(shear, [], [], [], 512); % Raw Power Spectrum
[ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A, shear, 512); % Clean Power Spectrum 
                                                                    % ------------------------------------- %
% ---- Conversion Freq -> Wavenumber ---- %
k0 = f0/Wm;
k1 = f1/Wm;

% corrected spectrum according to Macoun & Lueck 04
%ps0 = respcorr(ps0, k0);
ps1 = respcorr(ps1, k1);
% --------------------------------------- %


% ---- Dissipation rate calculation ---- %
epsilon = dewey_opti(ps1, k1, f1, visco);    

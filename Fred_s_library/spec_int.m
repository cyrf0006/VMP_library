function   epsilon = spec_int(shear, A, Wm, visco, fs);

% Spectral integration of shear to obtain dissip. rate of TKE. This
% function is used in Lo_Lt.
%    
% ex: eps1 = spec_int(sh1, [ax(I) ay(I) az(I)], Wm, visco, fs);
%
% --------------------------------------------------------------- %
    

%  ---- clean shear ---- %    
II = find(isnan(shear)==1);

if length(II)>round(length(shear)/4)
    epsilon = NaN;
else
    shear(II)=0; % pad NaNs with zeros and continue
% -------------------------------------------------------- %

    % --- Compute FFT & clean spectrum ---- %
    [ps0, f0] = pwelch(shear, [], [], [], fs); % Raw Power Spectrum
    [ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A, shear, fs); % Clean Power Spectrum 
    % ------------------------------------- %
    
    k0 = f0/Wm;
    k1 = f1/Wm;

    % corrected spectrum according to Macoun & Lueck 04
    ps0 = respcorr(ps0, k0);
    ps1 = respcorr(ps1, k1);
    % --------------------------------------- %


    % ---- Dissipation rate calculation ---- %
    epsilon = dewey_opti(ps1, k1, f1, visco);    

    
% $$$     % plot result
% $$$     loglog(k0, ps0)
% $$$     hold on
% $$$     loglog(k1, ps1, 'r')        
% $$$     title(sprintf('e = %d', epsilon))
% $$$     xlim([1e0 1e3])
% $$$     [phi,k] = nasmyth(epsilon, visco, length(ps1));
% $$$     loglog(k, phi, 'k')
% $$$     hold off
% $$$     pause
    
end


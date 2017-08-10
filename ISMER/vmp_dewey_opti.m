function eps0 = vmp_dewey_opti(ps1, k1, f1, viscosity)
% function epsilon = dewey_opti(ps, k, visco)
%
% Portion of original epsilon_VMP.m from R. Dewey, where the
% optimization for the best fit for spectral integration is
% done. This function returns epsilon for a specific bin.
    

% RKD (0707) a few changes below to correct and improve eps calcualtion
% first estimate small variance epsilon
eps0=1e-9; % RKD 07/07
lost=4.0; % start assuming 75% lost variance
nasfft=length(ps1);
df = f1(2)-f1(1);


for iloop=1:4, % RKD 07/07 increased the search for ideal integration limits

    % Note: These wavenumber limits KKSL/U are NON-dimensional (k/k_s)
    % RKD 07/07 modified selections below
    if eps0<1e-9,
        kksl=1.5e-2;kksu=3e-2; % lowest epsilon, shink integration band
    elseif eps0>=1e-9 & eps0<1e-8,
        kksl=1.5e-2;kksu=5e-2; % low epsilon, widen integration band
    elseif eps0>=1e-8 & eps0<1e-7,
        kksl=8e-3;kksu=7e-2; % high epsilon, widen integration band
    elseif eps0>=1e-7 & eps0<1e-6
        kksl=3e-3;kksu=9e-2; % higher epsilon, widen integration band
    elseif eps0>=1e-6;
        kksl=1e-3;kksu=1e-1; % highest epsilon, widen integration band
    end

    for ie=1:4,
        
        ks0=(eps0*lost/viscosity^3)^0.25;
        kl0=kksl*ks0;
        ku0=kksu*ks0; % find approximate integration limits 1e-3< k/ks <1e-1
        il=max([1 max(find(k1 < kl0))]); % index to lower integration limit
        iu=min([min(find(k1 > ku0)) length(k1)]);   % index to upper integration limit
        variance=sum(ps1(il:iu))*df; % integrate power spectrum from kl to ku m^-1
        eps0=7.5*viscosity*variance; % area under spect between limits
        epsn=eps0*lost;
        % RKD 07/07 modified determination of % lost below
        lost=1;
        for in=1:2,
            [uphi,uk]=nasmyth(epsn,viscosity,nasfft);
            uil=max([0 max(find(uk(1:end-2) < kl0))]) + 1;
            uiu=min([min(find(uk(3:end) > ku0)) length(uk)+1]) - 1;
            uvar=sum(uphi);

            % Here, should do the cumsum up to ks (we overshoot
            % here). However, the extra variance added is negligible
            sp=(cumsum(uphi)/uvar);  % WARNING: SOMETIMES DIVIDES BY ZEROS!!! -FC
            epsnan = 0;
            if uil==uiu | uil>=nasfft | uiu <=1 
                disp('prob. in dewey_opti.m')
                %keyboard; 
                disp(' Return epsilon = NaN')
                epsnan = 1;

                if uiu <=1 % problem reported 2012-11-05
                           % (2012_10_25, profile 20)
                    uiu = 1;
                end 
            end % what's happening...
            pll=sp(uil);
            plu=1-sp(uiu);
            lostn=pll + plu; % this is the percentage lost outside integration limits
            lost=1/(1-lostn);
            epsn=eps0*lost;  % updated estimate or epsilon from Nasmyth
        end % in=1:2,
            % RKD 07/07 modified code above
        lost=lost - (lost-2)/2;
        if epsnan == 1;
            eps0 = NaN;
        else            
            eps0=7.5*viscosity*variance*lost;  % final estimate of Epsilon
        end
    end % for ie=1:4
end % for iloop=1:4

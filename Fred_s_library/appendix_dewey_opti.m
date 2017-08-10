% To be run anywhere (like in /media/cyrf0006/Seagate2TB/cyrf0006_2014-03-26/WINDEX/data_processing/dewey_opti)

clear

%% Original load (for thesis):
%load /home/cyrf0006/WINDEX/data_processing/Robin_stuff/profile_2011-09-22_030.mat

%% Test done for Daniel on 2016-03-07:
%load /media/cyrf0006/Seagate2TB/cyrf0006_2014-03-26/WINDEX/data_processing/Robin_stuff/profile_2011-09-22_030.mat
load /media/cyrf0006/Seagate2TB/cyrf0006_2014-03-26/WINDEX/data_processing/Robin_stuff/profile_2012-10-25_019.mat

% ---> From vmp_shear_analysis.m
iplt = 0; %   iplt==0: No plotting done.

% Angle and viscosity
ax = ang2acc(pitch); ay = ang2acc(roll);

DENS = sw_dens(SBS, SBT, P);
[nu, mu] = viscosity(SBS, SBT, DENS);

% round frequencies
fs = round(fs);
FS = round(FS);

sh1 = vmp_shear_preanalysis(shear1, w, pitch, roll, p, 5);

% [eps1,p_eps1] = vmp_spectral_integration(sh1,[ax ay az],P, W, nu, fs, FS,'iplt',iplt,'zbin',zbin, 'explicit', explicit);
% function [epsilon,zfft] = vmp_spectral_integration(shear, A, z, W, visc, shearFreqHz, pressFreqHz, varargin)
A = [ax ay az];
zbin = 1;
shear = sh1;
z = P;
shearFreqHz = fs;
pressFreqHz = FS;
visc = nu;

% ---> From vmp_spectral_intergration.m

zs=interp1([1:length(z)],z,[1:length(shear)]'*(pressFreqHz/shearFreqHz));

z1 = ceil(min(abs(z)));  % top of profile
z2 = floor(max(abs(z)));  % bottom of profile
zfft = z1+zbin/2:zbin:z2-zbin/2; %reg. z prof. for FFT
epsilon = nan(size(zfft));

I=find(zs>=z1 & zs<=z2); % we keep only full bins
shear=shear(I);
zs = zs(I); 
A = A(I,:); % if shear is adjusted, size of A must be too


figure_appendix
I = [];
indices(1) = find(zfft==zzVec(1));
indices(2) = find(zfft==zzVec(2));
indices(3) = find(zfft==zzVec(3));


for i = indices

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
    ps10 = ps1;
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
    viscosity = mean(visc(Ifine));
    
    % ------> From dewey_opti.m
    %function eps0 = dewey_opti(ps1, k1, f1, viscosity)
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
            eps0_fc = eps0;

            ks0=(eps0*lost/viscosity^3)^0.25; % <--- 1/Lk (inverse of kolmogorov lengthscale)
                                              % !!!!!!!!! SHOULD WE PUT 2PI FACTOR??
            %ks0 = ks0/(2*pi);
            kl0=kksl*ks0;
            ku0=kksu*ks0; % find approximate integration limits 1e-3< k/ks <1e-1
            il=max([1 max(find(k1 < kl0))]); % index to lower integration limit
            iu=min([min(find(k1 > ku0)) length(k1)]);   % index to upper integration limit
            variance=sum(ps1(il:iu))*df; % integrate power spectrum from kl to ku m^-1
            eps0=7.5*viscosity*variance; % area under spect between limits
            epsn=eps0*lost;
            eps1_fc = epsn; % -FC (cf. thesis Appendix I)
            
            % plot_Dewey
            % RKD 07/07 modified determination of % lost below
            lost=1;
            for in=1:2,
                [uphi,uk]=nasmyth(epsn,viscosity,nasfft);               
                uil=max([0 max(find(uk(1:end-2) < kl0))]) + 1;
                uiu=min([min(find(uk(3:end) > ku0)) length(uk)+1]) - 1;
                %uil = max(find(uk<=kl0)); % -FC
                %uiu = min(find(uk>=ku0));
                uvar=sum(uphi);
                
                sp=(cumsum(uphi)/uvar);  % WARNING: SOMETIMES DIVIDES BY ZEROS!!! -FC
                epsnan = 0;
                if uil==uiu | uil>=nasfft | uiu <=1 
                    disp('prob. in dewey_opti.m')
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
                epsn=eps0*lost;
                          
                 % Honestly, I dont underestand the next line, but I've used it during my thesis
                lost=lost - (lost-2)/2;
                eps2_fc = 7.5*viscosity*variance*lost; % -FC (cf. thesis Appendix I)
                plot_Dewey
                %[iloop, ie, in]
                %eps2_fc
                %pause
            end % in=1:2,
                % RKD 07/07 modified code above
            
            
            %lost=lost - (lost-2)/2; % <------------- why this???????????????
            if epsnan == 1;
                eps0 = NaN;
            else            
                eps0=7.5*viscosity*variance*lost;  % final estimate of Epsilon
            end            

        end % for ie=1:4
    end % for iloop=1:4
    
    
    % lostn
    epsilon(i) = eps0;
    disp('You might want to rename Dewey_algo.eps')
    keyboard
% $$$     % Instrument resolution ~1e-10
% $$$     if epsilon(i)<1e-10;
% $$$         epsilon(i) = NaN;
% $$$     end
% $$$     
% $$$     % Test eroneous epsilon
% $$$     if epsilon(i) > 1e-5
% $$$         figure(1)
% $$$         clf
% $$$         
% $$$         subplot(121)
% $$$         loglog(k0, ps0)
% $$$         hold on
% $$$         loglog(k1, ps1, 'r')        
% $$$         title(sprintf('z = %0.1f m   e = %d', zfft(i), epsilon(i)))
% $$$         xlim([1e0 1e3])
% $$$         [phi,k] = nasmyth(epsilon(i), visco, length(ps1));
% $$$         loglog(k, phi, 'k')
% $$$         hold off
% $$$         %axis([1e0 1e3 1e-8 1e-2])
% $$$ 
% $$$         subplot(122)
% $$$         ind1 = I(1)-4*length(I);
% $$$         ind2 = I(end)+4*length(I);
% $$$ 
% $$$         if ind1<1
% $$$             ind1 = 1;
% $$$         elseif ind2 > length(zs)
% $$$             ind2 = length(zs);
% $$$         end
% $$$         plot(shear(ind1:ind2), zs(ind1:ind2))
% $$$         hold on
% $$$         plot(shear(I), zs(I), 'r')
% $$$         hold off
% $$$         set(gca, 'ydir', 'reverse')
% $$$         
% $$$         disp('Epsilon is very high. Do you want a check this?')
% $$$         answer=0;
% $$$         while answer==0
% $$$             R1 = input('Do you accept the profile? (y/n) ', 's');
% $$$             if strcmp('y',R1)==1 
% $$$                 answer=1;
% $$$             elseif strcmp('n',R1)==1
% $$$                 epsilon(i) = NaN;
% $$$                 answer=1;
% $$$             else
% $$$                 disp('bad answer! please type y or n')
% $$$             end
% $$$         end
% $$$     end % ----------  end test ----- %
% $$$ 
% $$$            
    
    
end


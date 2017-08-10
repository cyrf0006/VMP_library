function compute_vars(no_profile)

% function compute_vars(no_profile)
%
% usage ex: compute_vars(11)
% 
% where no_profile is the number of profile to analyse
%
% This function will save a .MAT for each profile with different variables
% in it. This function is easy to modify to compute new variables. For now,
% it uses profilexxx.mat and eps_profilexxx.mat to compute:
%
% - P_vars: Pressure vector corresponding to computed vaiables
% - eps: mean dissipation of TKE
% - N2: buoyancy freq.
% - SIGb: binned sigmat_t
% - K_rho: dissip. coeff.
%
% all variables are binned following the dissipation of the TKE wich is in 
% eps_profilexxx.mat and will be saved in binned_varsxxx.mat
%
% author: F. Cyr - 25/01/2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % some constants
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.035e3;%kg/m^3
GAMMA = 0.2; %mixing efficiency

for profile = 1:no_profile

    % the name of the profile
    if profile<10
        data_fname = sprintf('profile00%d', profile);
        eps_fname = sprintf('eps_profile00%d', profile);
        disp(sprintf('profile00%d', profile));
        outfile = sprintf('binned_vars00%d', profile);
    else
        if profile<100
            data_fname = sprintf('profile0%d', profile);
            eps_fname = sprintf('eps_profile0%d', profile);
            disp(sprintf('profile0%d', profile));
            outfile = sprintf('binned_vars0%d', profile);          
        else %profile>100
            data_fname = sprintf('profile%d', profile);
            eps_fname = sprintf('eps_profile%d', profile);
            disp(sprintf('profile%d', profile));
            outfile = sprintf('binned_vars%d', profile);
        end
    end
 
    load(data_fname)
    load(eps_fname)
    
    
    % -- Combination of the 2 eps_profiles... sort of average ignoring NaN
    p_eps = nanmean([p_eps1; p_eps2]);
    eps = nanmean([eps1; eps2]);  
    
    % -- Computing density -- %
    DENS = sw_dens(SBS, SBT, P); % density profile
    SIG_T = DENS-1000;

    % -- Sorting density -- %
    [SIG_SORTED, I] = sort(SIG_T'); % must be row vector   

    % -- average bins -- %
    Pb = p_eps; % Keep the same discretization as eps
    dP = p_eps(2)-p_eps(1);
    
    for i = 1:length(Pb)
        I = find( P >= Pb(i)-dP/2 &  P <= Pb(i)+dP/2); %find data both side of the bin center
        SIGb(i) = mean(SIG_SORTED(I));
        %Tb(i) = mean(SBT(I));
        %Sb(i) = mean(sal(I));
        
    end
    
    
    % -- Brunt-Vaisala -- %
     drhodz = diff(SIGb);
    
    N2 = g/rho_0*drhodz/dP;

    J = 1:length(N2); %indice of variable (length(N2) = length(eps)-1) 

    EPS = eps(J);
    K_rho = GAMMA.*eps(J)./N2;
    P_vars = p_eps(J);
    SIGb = SIGb(J);    
    
    save(outfile, 'P_vars', 'EPS', 'N2', 'SIGb', 'K_rho')
    
    clear P_vars N2 SIGb K_rho Tb drhodz J P_vars Pb eps p_eps DENS SIG_T

end
        
    
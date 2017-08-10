function diffus_model_noflux(temp_file, sal_file, time0, diffus_value, diffus_method, forcing, dt, DT, L)

% function diffus_model(T0, S0, P0, diffus_value, diffus_method, forcing, dt, DT, L)
%
% Small diffusivity model which compute the evolution of
% temperature and salinity profiles in time. The model consider
% also the heat flux computed at surface. 
% K_model which compute evolution of variable D0 for constant K.
% D0 can be density, salinity or temperature (just need to adjust
% xlim in plot).
% P0 is the corresponding pressure, test_K is the constant K value,
% dt, DT and L are respectively the time step of the model for
% calculation, the time step for the outputplot and the total time
% of the simulation (all in seconds)
%
% usage ex: >> diffus_model('T_climato_04.dat', 'S_climato_04.dat' , datenum(999, 4, 22), 5e-5, 'constant_k', 'daily_forcing.dat', 100, 3600*24*30, 3600*24*30*6)
%           >> monthly_clim_model([4 5 6 7 8 9 10 11], 'T_diffus_daily_5p5e-5.dat', 'N_diffus_daily.dat', 'T_diffus_daily_Kobsmin.dat', 'T_diffus_daily_Kobsmax.dat')

% When diffus_method is variable_k or variable_e, diffus_value must
% be in the format [P, K] and the size Nx2
% Note: daily_heatflx starts at april 22th and finishes at nov. 18th
    
    
% Author: F. Cyr (feb. 2011, inspired from diffus_model.m)
%
% Modifications: 
%  - may 5: Add a computation of integrated heat content on the CIL
% ------------------------------------------------------------ %

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- Preliminaries ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------- Temp and salinity -------- %
temp = load(temp_file); % temp = [P T_ave (T_error or 95%CI, optionnal)]
if size(temp,2)~=2 & size(temp,2)~=3 & size(temp,2)~=4 
    error('myApp:argChk','temperature file not in good format: [P temp]')
end
T0 = temp(:,2);
P0 = temp(:,1);

sal = load(sal_file); 
if size(sal,2)~=2 & size(sal,2)~=3 & size(sal,2)~=4
    error('myApp:argChk','salinity file not in good format: [P sal]')
end
S0 = sal(:,2);
% -------------------------------------- %


% ------- some constants -------- %    
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.035e3;%kg/m^3
GAMMA = 0.2; %mixing efficiency
freeze_pt = -1.8; 
cp = 3.99; %Kj/Kg/K
CIL_def = 1; % max Temp for CIL definition

dz = P0(2)-P0(1);

C = dt/dz^2; % coefficient for model

% pressure limits
Zmin = round(min(P0));
Zmax = round(max(P0));

P_top = 1; % top boundary of the model (where Q is forced)
% -------------------------------------- %

% --------- surface forcing---------- %
forcing = load(forcing); % forcing = [n SST SSS]
nf = forcing(:,1);
SST = forcing(:,2);
SSS = forcing(:,3);

%Make sure time are corresponding
I = find(nf>=time0 & nf<=time0+L/3600/24); 
if ~isempty(find(isnan(SST(I))==1))==1
    disp('There is at least one NaN in model forcing. Are you sure times are corresponding?')
    disp('  Extrapolation to closest value to remove NaNs')
    %disp('Type "return" to continue anyway')
    
    % find NaNs in forcing
    J = find(isnan(SST(I))==1);
    JJ = find(~isnan(SST(I))==1);
    % replace 1st value
    if J(1)<JJ(1)
        SST(I(J(1)))=SST(I(JJ(1)));
        SSS(I(J(1)))=SSS(I(JJ(1)));
    end
    % replace last value
    if J(end)>JJ(end)
        SST(I(J(end)))=SST(I(JJ(end)));
        SSS(I(J(end)))=SSS(I(JJ(end)));
    end
    
    % remove NaNs in forcing
    J = find(isnan(SST(I))==1);
    SST(I(J))=[];
    SSS(I(J))=[];
    nf(I(J))=[];
end



interp_time = time0:dt/3600/24:time0+L/3600/24;
SST_forc = interp1(nf, SST, interp_time);
SSS_forc = interp1(nf, SSS, interp_time);

% -------------------------------------- %

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- How to diffuse ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Determine which method -- %
if(diffus_method=='constant_k')
    K(1:length(P0)-1) = diffus_value;
    K=K';
    
elseif(diffus_method=='variable_k')
    
    % interpolate K at P0 pressure level
    % extrapolation when K value are missing
    P_to_itp = diffus_value(:,1);
    K_to_itp = diffus_value(:,2);
    if P_to_itp(1)>P0(1) %clip K at beginning
        P_to_itp = [P0(1); P_to_itp];
        K_to_itp = [K_to_itp(1); K_to_itp];
    end
    if P_to_itp(end)<P0(end) %clip K at the end
        P_to_itp = [P_to_itp; P0(end)];
        %K_to_itp = [K_to_itp; K_to_itp(end)];
        K_to_itp = [K_to_itp; min(K_to_itp)];
    end
    
    if size(diffus_value,2)==2
        %K = interp1(diffus_value(:,1), diffus_value(:,2),P0(1:length(P0)-1));
        K = interp1(P_to_itp, K_to_itp,P0(1:length(P0)-1));
    else
        error('myApp:argChk','diffus_value not in good format [pressure value]')
    end

elseif(diffus_method=='constant_e')
    
    EPS(P_top:length(P0)-1) = diffus_value;
    EPS=EPS';
    
elseif(diffus_method=='variable_e')
    
  
    % interpolate EPS at P0 pressure level
    if size(diffus_value,2)==2        
        EPS = interp1(diffus_value(:,1), diffus_value(:,2) , P0(P_top:length(P0)-1));

        if ~isempty(find(isnan(EPS)==1)) %remove NaN
            I = find(~isnan(EPS)==1);
            EPS(1:I(1))=EPS(I(1));
            EPS(I(end):end)=EPS(I(end));
        end
        
    else
        error('myApp:argChk','diffus_value not in good format [pressure value]')
    end

else
    error('myApp:argChk','diffus_method must be either ''constant_k'',''variable_k'',''constant_e'' or ''variable_e''')
end

% -- initial profiles -- %
T=T0(P_top:end);
S=S0(P_top:end);
P=P0(P_top:end);
[DENS, I] = sort(sw_dens(S, T, P));
T=T(I);        
S=S(I);

% -- Brunt-Vaisala -- %
drho = diff(DENS);
N2 = g/rho_0*drho/dz;

% -- K computation if epsilon method -- %
if(diffus_method=='constant_e' | diffus_method=='variable_e')
    K = GAMMA.*EPS./N2; %K from epsilon
end        
    
% -- boundary conditions -- %
K1 = [K; 0]; % no flux at bottom
K2 = [0; K]; % no flux at top            
%K2 = [1; K]; % heatflx at top   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- fisrt time step ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 1;
N(count) = time0;
Tmat(:,count) = T;
Smat(:,count) = S;

count = count+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- time loop / core of model----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Entering the model')

for t = 1:L/dt %loop on time
    
    %keyboard
    %T_mk = [atm_heatflx(t); T(1:length(T)-1)]; % T_k-1
    T_mk = [0; T(1:length(T)-1)]; % T_k-1
    T_pk = [T(2:length(T)); 0]; % T_k+1    
    T_av = T + C.*(K1.*(T_pk - T) - K2.*([0; T(2:end)] - T_mk));
    T_av(1) = SST_forc(t+1); % impose T at 1st b in
    
    %S_mk = [atm_salinflx(t); S(1:length(S)-1)]; % T_k-1
    S_mk = [0; S(1:length(S)-1)]; % T_k-1
    S_pk = [S(2:length(S)); 0]; % T_k+1    
    S_av = S + C.*(K1.*(S_pk - S) - K2.*([0; S(2:end)] - S_mk));
    S_av(1) = SSS_forc(t+1); % impose S at 1st b in
  
    
    if t*dt/DT == count-1; % plot every DT (s)
        
% $$$         figure(1)
% $$$         clf
% $$$         plot(T_av, P);
% $$$         %axis([1020 1030 Zmin Zmax])
% $$$         xlabel('temprature')
% $$$         ylabel('depth(m)')
% $$$         axis([-1 6 Zmin Zmax])
% $$$         set(gca, 'ydir', 'reverse')
% $$$         title(sprintf('output %d', count))
% $$$         pause(0.1)  
% $$$         
% $$$         figure(2)
% $$$         clf
% $$$         plot(S_av, P);
% $$$         %axis([1020 1030 Zmin Zmax])
% $$$         xlabel('salinity')
% $$$         ylabel('depth(m)')
% $$$         axis([26 35 Zmin Zmax])
% $$$         set(gca, 'ydir', 'reverse')
% $$$         title(sprintf('output %d', count))
% $$$         pause(0.1)  
        
        Tmat(:,count) = T_av;
        Smat(:,count) = S_av;
        N(count) = N(count-1)+DT/3600/24;

        count = count +1;
 
    end
        
    T = T_av;
    S = S_av;
    
    % -- K re-computation if epsilon method -- %
    if(diffus_method=='constant_e' | diffus_method=='variable_e')
        [DENS, I] = sort(sw_dens(S, T, P));
        T=T(I);
        S=S(I);
        
        drho = diff(DENS);
        N2 = g/rho_0*drho/dz;
         
        K = GAMMA.*EPS./N2; %K from epsilon
        
        K1 = [K; 0]; % no flux at bottom
        K2 = [0; K]; % no flux at top               
        %K2 = [1; K]; % heatflx at top   
    end        
    
        
end


% $$$ % Output idea
% $$$ figure(2)
% $$$ clf
% $$$ imagesc(N, P, Tmat)
% $$$ datetick('x', 3)
% $$$ 
% $$$ figure(3)
% $$$ clf
% $$$ imagesc(N, P, Smat)
% $$$ datetick('x', 3)


% Save result which will be read in post-treatment
dlmwrite('T_diffus_daily.dat', Tmat,'delimiter',' ','precision',6)
dlmwrite('S_diffus_daily.dat', Smat,'delimiter',' ','precision',6)
dlmwrite('N_diffus_daily.dat', N,'delimiter',' ','precision',6)


function diffus_model_O2(O2_file, time0, diffus_value, diffus_method, hflx, dt, DT, L, topcondition)
    
% Here we use a 1D diffusion model to double check Bourgault's
% diffusion model on dissolved oxygen. 
%
% see exectute_Kmodel.m for more info...
%
% Author: F. Cyr (nov. 2011)
%
% ------------------------------------------------------------ %

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- Preliminaries ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------- Temp and salinity -------- %
tracer = load(O2_file); % temp = [P TEMP]
if size(tracer,2)~=2
    error('myApp:argChk','temperature file not in good format: [P temp]')
end
T0 = tracer(:,2);
P0 = tracer(:,1);
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


% --------- atmospheric fluxes ---------- %

heat = load(hflx); % hflx = [N HFLX STD_HFLX]
                   %H = heat(:,2)/cp/rho_0/1000; %/1000 because cp is in kj
H = heat(:,2)*-1;
heat_time = heat(:,1);                                     
interp_time = [time0:dt/3600/24:time0+L/3600/24]';
atm_heatflx = interp1(heat_time, H, interp_time);

atm_heatflx = atm_heatflx/(P_top*dz)/1000; %1000L/m3

% Make sure time are corresponding
if ~isempty(find(isnan(atm_heatflx)==1))==1
    disp('There is at least one NaN in atm_heatflx. Are you sure times are corresponding?')
    disp('Type "return" to continue anyway')
    keyboard
end

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
        K_to_itp = [K_to_itp; K_to_itp(end)];
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
P=P0(P_top:end);

% -- K computation if epsilon method -- %
if(diffus_method=='constant_e' | diffus_method=='variable_e')
    K = GAMMA.*EPS./N2; %K from epsilon
end        
    
% -- boundary conditions -- %
K1 = [K; 0]; % no flux at bottom
K2 = [1; K]; % heatflx at top   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- fisrt time step ----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;
N(count) = time0;
Tmat(:,count) = T;
%Smat(:,count) = S;

count = count+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- time loop / core of model----------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:L/dt %loop on time
    
    T_mk = [atm_heatflx(t); T(1:length(T)-1)]; % T_k-1
    T_pk = [T(2:length(T)); 0]; % T_k+1    
    T_av = T + C.*(K1.*(T_pk - T) - K2.*([0; T(2:end)] - T_mk));

    T_av(end) = topcondition; % impose T at 1st b in

    
    if t*dt/DT == count-1; % plot every DT (s)
        
% $$$         figure(1)
% $$$         clf
% $$$         plot(T_av, P);
% $$$         %axis([1020 1030 Zmin Zmax])
% $$$         %xlabel('temprature')
% $$$         %ylabel('depth(m)')
% $$$         %axis([-1 6 Zmin Zmax])
% $$$         set(gca, 'ydir', 'reverse')
% $$$         title(sprintf('output %d', count))
% $$$         pause(0.25)  
        
        Tmat(:,count) = T_av;
        N(count) = N(count-1)+DT/3600/24;

        count = count +1;
 
    end
        
    T = T_av;
    
end


% Output idea
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 16 10])
O2min = 60;
O2max = 200;
V = O2min:5:O2max;
ucurrent = 0.01; %m/s

% change variable t -> x
X = abs(N-N(end))*86400*ucurrent/1000;
P = flipud(P);

%subplot(3, 1, 1);
contourf(X, P, Tmat)
hold on
[c, h] = contour(X, P, Tmat, V, 'color', 'k');
set(h,'ShowText','on', 'textstep', get(h,'LevelStep')*2);
hold off
set(gca, 'ydir', 'reverse')
cmap = lbmap(128,'RedBlue');
colormap(flipud(cmap));
caxis([O2min O2max])

keyboard
save test_a.mat X P Tmat


% After you can run >> contour_O2_testcases
% to see output!


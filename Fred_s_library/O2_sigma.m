function O2_sigma(NAMELIST)

% function O2_sigma(NAMELIST)
%
% 2D advection-diffusion model in topography following (sigma)
% coordinates for any passive tracer c. The model solves the
% following equation for advection-diffusion in sigma coordinates
% with constant depth H in time (no free surface) :
%
% dc/dt + H d(uc)/dx + d(w'c)/dsigma - H d(Kx dc/dx)/dx    (1)
%      - 1/H d(Kz dc/dsigma)/dsigma = 0
%
% where Kx and Kz are the horizontal and vertical diffusivities, u
% is the horizontal velocity and w' (named wp in equations) is the
% "vertical" velocity in sigma corrdinate (unit s^{-1}) that we
% obtain from continuity:
%
% d(Hu)/dx + d(Hw')/dsigma = 0;                            (2)
%
% The model uses wether true horizontal diffusion (Zängl
% 2002) or diffusion along sigma lines (not recommended).
% Model options could be edited through NAMELIST.
%
% usage ex: O2_model_sigma(config.nml)
% where 'config.nml' looks like:
%
% $$$ % ------------------------------------- %
% $$$ % Spatiotemporal grid
% $$$ dt = 120;%sec
% $$$ dx = 1000; %m
% $$$ dsigma = .02;
% $$$ z0 = 150; %m (rigid lid depth)
% $$$ 
% $$$ % display figure every 'disp_freq' timesteps    
% $$$ DISP_FREQ = 500;
% $$$     
% $$$ % bathymetry file
% $$$ load('H_of_x.mat');
% $$$ 
% $$$ % Diffusivities
% $$$ Kx = 1600; %m2/s
% $$$ Kz = 1e-4; %m2/s
% $$$ 
% $$$ % Bottom and water column respiration
% $$$ Rwc = -19.6; %mmol/m3/yr
% $$$ Rsed = -3540; %mmol/m2/yr
% $$$ 
% $$$ % Horizontal velocity
% $$$ u0 = -.003; %m/s at x=L
% $$$ 
% $$$ % Initial conditions
% $$$ INIT_COND = 170; %mmol/m3;
% $$$ c_xL = 170; %mmol/m3 (b.c. at x = L)
% $$$ c_z0 = 170; %mmol/m3 (b.c. at z = z0);
% $$$ 
% $$$ % Use true horizontal Diffusion or along sigma diffusion
% $$$ % This is a known problem in sigma-coordinate (0 or 1). 
% $$$ % we recommand TRUE_DIFFUSION=1
% $$$ TRUE_DIFFUSION = 1;
% $$$ 
% $$$ % Hybrid diffusion =0 for none, =1, =2, =3, ... for the number of
% $$$ % cell above the bottom where it will be applied
% $$$ HYBRID_DIFFUSION = 1; 
% $$$ 
% $$$ % Since computing true diffusion is time consuming (there is an
% $$$ % interpolation of points between cartesian and sigma grid at each
% $$$ % time step), we developped an alternative scheme which converge
% $$$ % faster. The idea is that we interpolate only once at the
% $$$ % beginning, and the wrelative contribution of the nearest terms
% $$$ % are kept in memory. 
% $$$ % Will have no effect is TRUE_DIFFUSION = 0;
% $$$ % we recommand FASTER_SCHEME = 1;
% $$$ FASTER_SCHEME = 1;
% $$$ % ---------------------------------------------------- %
%
% The output will be a .mat file with final O2 field, together with
% corresponding X and Z. Visualization of the result can be done
% with a simple pcolor or contourf function (imagesc will show
% results in sigma coordinate, not z):
% > pcolor(X,Z,c_t)
%
% Model grid (in sigma domain, positive downward):
%
%          x=0                              x=L
% s=1    ---^-----^-----^-----^---        ---^---
% |      >  c  >  c  >  c  >  c  >        >  c  >
% |      ---^-----^-----^-----^---  (...) ---^---
% |      >  c  >  c  >  c  >  c  >        >  c  >
% |      ---^-----^-----^-----^---        ---^---
% |
% |                  (...)
% |
% s=0    >  c  >  c  >  c  >  c  >        >  c  >
%        ^-----^-----^-----^-----^        ---^---
%
% c: tracer
% ^: vertical velocity (w=0 at s=0 and s=1)
% >: horizontal velocity (u = u0 at x=L)
% s: sigma level (s=0 at surface s=1 at bottom)
%
%
% Author: F. Cyr - October 2012
%   for any bugs or comment: frederic.cyr@uqar.qc.ca
% ------------------------------------------------------------- %


% ---------------------------- Initialization ------------------------------------ %
% model namelist (topography is also loaded through this)
fid = fopen(NAMELIST);
tline = fgetl(fid);
while ischar(tline)
    eval(tline);
    tline = fgetl(fid);
end
fclose(fid);

if TRUE_DIFFUSION~=0;
    TRUE_DIFFUSION=1; % make it default
end

if ~exist('FASTER_SCHEME') == 1
    FASTER_SCHEME = 0;
end


x_u = x0(1):dx:x0(end); % u velocity grid
x_c = x_u(1:end-1)+dx/2; % tracer grid
H = round(interp1(x0, H0, x_u))-z0; %bathymetry
H(2)=H(1); % Because of Neumann condition, same depth at i=1,2

sigma = 0:dsigma:1; % sigma coord.
Hc = min(min(H));

for i = 1:length(H) %sigma levels
    z(:,i) = H(i).*sigma;
end

Rwc = Rwc/(365*86400); %mmol/m3/s
Rsed = Rsed/(365*86400); %mmol/m2/s

% Compute U on the grid
% d(UH)/dx = 0 -> UH = cst
UH = u0*H(end);
U = UH./H; %U(x) = cst/H(x)
for i =1:length(x_u)
    u_grid(1:length(sigma)-1, i) = U(i);
end

[Fx, Fy] = gradient(z);
Hz = Fy./dsigma;

% initializing the grid (c is the tracer)
c_grid = zeros(length(sigma)-1, length(x_c));
for i = 1:length(x_c)
    c_grid(1:length(sigma)-1, i) = INIT_COND;
end

% ----- > w calculation
u_ip1 = [u_grid(:,2:end) zeros(length(sigma)-1, 1)]; %spatial shift
Hz_ip1 = [Hz(:,2:end) zeros(length(sigma), 1)];
Hz_sp1 = [Hz(2:end,:); zeros(1, length(x_u))];


format long
% d(Hu)/dx + d(Hw')/dsigma = 0;                            (2)

wp_grid = zeros(size(c_grid,1)+1, size(u_grid,2)); % w =0 will remains at bottom
for s = 1:length(sigma)-1
    dHudx          = (Hz_ip1(s,:).*u_ip1(s,:) - Hz(s,:).*u_grid(s,:))./dx;
    wH             = Hz(s,:) .* wp_grid(s,:);
    wp_grid(s+1,:) = (dHudx * dsigma + wH) ./ Hz_sp1(s,:);    
end
wp_grid(:,end) = 0;

% bring w to c_grid
wp_grid = (wp_grid(:, 2:end) + wp_grid(:, 1:end-1))./2;


% ----> cell "volume"
Dz = diff(z, 1, 1);
Dz = (Dz(:,1:end-1)+Dz(:,2:end))./2;
Dx = diff(x_u);
Dx = meshgrid(Dx, sigma);
cellVolume = Dx(1:end-1,:).*Dz; % in m^2


% Because of time step used, we need to keep 2 timestep in memory
c_t = c_grid;
c_tm1 = c_grid;

% Grid for plotting in z space
X = meshgrid(x_c, c_t(:,1));
dz = diff(z(1:2,:));
Z = zeros(size(X));
for i = 1:size(X,2)
    Z(:,i) = 1:dz(i):z(end,i); %start at 1
end

% bring all grid in same dimensions
w_t = (wp_grid(2:end, :) + wp_grid(1:end-1, :))./2;
Hz = (Hz(:,2:end) + Hz(:,1:end-1))./2;
Hz = (Hz(2:end, :) + Hz(1:end-1, :))./2;
u_t = (u_grid(:, 2:end) + u_grid(:, 1:end-1))./2;
% --------------------------- end initialization --------------------------------------- %

% ------------- few initialization plots -------------- %
figure(1)
clf
plot(x_u/1000, z+150, 'k')
set(gca, 'ydir', 'reverse')
title('Sigma-contours')
ylabel('depth (m)')
xlabel('Seaward distance from STN 25 (km)')

figure(2)
clf
pcolor(X/1000, Z+150, u_t)
shading interp
set(gca, 'ydir', 'reverse')
title('Horizontal velocity (m/s)')
ylabel('depth (m)')
xlabel('Seaward distance from STN 25 (km)')
colorbar

% W in z space (to verify if the field looks good)
dzdx = diff(z, 1, 2)./dx;
W = Hc*wp_grid(1:end-1,:)+u_grid(:,1:end-1).*dzdx(2:end,:);

figure(3)
clf
pcolor(X/1000, Z+150, W)
shading interp
set(gca, 'ydir', 'reverse')
title('Vertical velocity (m/s)')
ylabel('depth (m)')
xlabel('Seaward distance from STN 25 (km)')
colorbar

figure(4)
clf
pcolor(X/1000, Z+150, c_t)
shading interp
set(gca, 'ydir', 'reverse')
title('Initial condition for tracer')
ylabel('depth (m)')
xlabel('Seaward distance from STN 25 (km)')
colorbar

disp(['Please verify initalization plots and press any key to ' ...
    'continue'])
pause

figure(1)
clf
figure(2)
clf
close(figure(3))
close(figure(4))

disp('simulation started')
% --------------------------------------------------------------- %

% -------------------- Discretization of constant terms -------------------------------- %
% spatial shift (c_t_im1 = c(t, i+1, s), etc.)
u_t_im1 = [zeros(length(sigma)-1, 1) u_t(:,1:end-1)];
u_t_ip1 = [u_t(:,2:end) zeros(length(sigma)-1, 1)];
w_t_sp1 = [w_t(2:end,:); zeros(1, length(x_c))];
w_t_sm1 = [zeros(1, length(x_c)); w_t(1:end-1,:)];

% for us, Hz is constant in time, but in case...
% for simplicity and consistency (at t=t):
H_t = Hz;
H_t_im1 = [zeros(length(sigma)-1, 1) Hz(:,1:end-1)];
H_t_ip1 = [Hz(:,2:end) zeros(length(sigma)-1, 1)];
H_t_sm1 = [zeros(1, length(x_c)); Hz(1:end-1,:)];
H_t_sp1 = [Hz(2:end,:); zeros(1, length(x_c))];
H_tm1 = Hz;
H_tm1_ip1 = H_t_ip1;
H_tm1_im1 = H_t_im1;
H_tm1_sp1 = H_t_sp1;
H_tm1_sm1 = H_t_sm1;

wgrid_t_sm1 = wp_grid(1:end-1,:);
wgrid_t_sp1 = wp_grid(2:end,:);


% For true diffusion scheme (Zängl, 2002)
if TRUE_DIFFUSION
    zc = diff(z,1)/2 + z(1:end-1, :);
    zc = diff(zc,1,2)/2 + zc(:,1:end-1);
    zc_ip1 = [zc(:,2:end) zc(:,end)];
    zc_im1 = [zc(:,1) zc(:,1:end-1)];
    
    if FASTER_SCHEME % Initialization of faster scheme
        CHCK_ERROR = 1;
    
        % This section creates two index matrices kFront and kBack and two
        % weight matrices weightBack and weightFront. These matrices
        % contain the vertical index of the closest level on each side of 
        % each cell and the weight to apply for the linear interpolation
        % carried out below in the main time loop. This is required for the 
        % truly horizontal diffusion scheme.
        
        [kmax imax] = size(zc);
        
   
        % Initialization
        for i = 1:imax
            kFront(1:kmax,i) = [1:kmax];
            kBack(1:kmax,i)  = [1:kmax];
        end
        weightBack(1:kmax,1:imax)  = 1.0;
        weightFront(1:kmax,1:imax) = 1.0;
        
        for i = 2:imax-1
            for k = 2:kmax
                
                % Find all indices associated with lower depth on each side
                % of the cell.
                kBackTmp    = find(zc(:,i-1) < zc(k,i));
                kFrontTmp   = find(zc(:,i+1) < zc(k,i));
                
                % Only keep the last index. We know from these matrices
                % that the interpolation should be done between 
                % kBack(i,k) and kBack(i,k + 1) 
                % and, equivalently, between 
                % kFront(i,k) and kFront(i,k + 1)
                % These matrices are used below within the main time loop. 
                kBack(k,i)  = kBackTmp(end);
                kFront(k,i) = kFrontTmp(end);

                kb = kBack(k,i);
                if kb == kmax
                    weightBack(k,i) = NaN;
                else
                    dzBack           = zc(kb + 1, i - 1) - zc(kb, i - 1);
                    dzBackUp         = zc(k,i)           - zc(kb, i - 1);
                    weightBack(k,i)  = 1.0 - dzBackUp / dzBack;
                end

                kf = kFront(k,i);
                if kf == kmax
                    weightFront(k,i) = NaN;
                else
                    dzFront          = zc(kf + 1,i + 1) - zc(kf,i + 1);
                    dzFrontUp        = zc(k,i)          - zc(kf,i + 1);
                    weightFront(k,i) = 1.0 - dzFrontUp / dzFront;
                end
                
            end
        end
                
    end %if faster_scheme
end %if true_diffusion
% ------------------------ end discretization of constant terms ------------------------------ %

% ------------------------------- MAIN LOOP ----------------------------------- %
no_cv = 1;
timeStep = 1;
M = INIT_COND;
tic


timeCounter = 0;

dsigma = -dsigma; % Inverse z axis

while no_cv
    
    timeCounter = timeCounter + 1;
    
    % --- time-variable terms discretization (only remains c, tracer concetration) --- %
    %spatial shift (c_t_im1 = c(t, i+1, s), etc.)
    c_t_im1   = [zeros(length(sigma)-1, 1) c_t(:,1:end-1)];
    c_t_ip1   = [c_t(:,2:end) zeros(length(sigma)-1, 1)];
    c_t_sm1   = [zeros(1, length(x_c)); c_t(1:end-1, :)];
    c_t_sp1   = [c_t(2:end,:); zeros(1, length(x_c))];
    c_tm1_im1 = [zeros(length(sigma)-1, 1) c_tm1(:,1:end-1)];
    c_tm1_ip1 = [c_tm1(:,2:end) zeros(length(sigma)-1, 1)];
    c_tm1_sm1 = [zeros(1, length(x_c)); c_tm1(1:end-1, :)];
    c_tm1_sp1 = [c_tm1(2:end,:); zeros(1, length(x_c))];
    % -------------------------------------------------------------------------------- %
    
    % --- terms calculation with B.C.(A = advec, D = diffus) --- %    
    % 1) advective terms (centered scheme)
    A1        = - (H_t_ip1.*u_t_ip1.*c_t_ip1 - H_t_im1.*u_t_im1.*c_t_im1)./(2*dx);
    % not good at i = end, but doesnt matter since c_i=end = 170
    A2        = - c_t.*H_t.*(wgrid_t_sp1 - wgrid_t_sm1)./(2*dsigma);
    A3        = - w_t.*H_t.*(c_t_sp1 - c_t_sm1)./(2*dsigma);
    A3(1,:)   = - w_t(1,:).*H_t(1,:).*(c_t_sp1(1,:) - c_t(1,:))./(2*dsigma); %forward
    A3(end,:) = - w_t(end,:).*H_t(end,:).*(c_t(end,:) -  c_t_sm1(end,1))./(2*dsigma);%backward 
    A = A1+A2+A3;
    
    % 2) diffusion terms (two options for horizontal...)
    % 2.1) true horizontal diffusion
    if TRUE_DIFFUSION
        cstar_t = c_tm1;
        
        if FASTER_SCHEME % faster scheme, but hard to understand...
            
            % Initialization. Only needed in the first iteration.   
            if timeCounter == 1
              cstar_t_im1 = c_tm1;
              cstar_t_ip1 = c_tm1;
            end
            
            % This is where the linear interpolation is carried out for
            % the truly horizontal diffusion scheme. 
            for i = 2:imax - 1
                for k = 1:kmax

                    % Intermediate variables kb and kf only created 
                    % for clarity.
                    kb = kBack(k,i);
                    kf = kFront(k,i);
                    
                    % If kb == kmax means that the cell intersects the
                    % bottom. Then simply set to NaN. The code below check
                    % for these NaNs and poses the one-sided diffusion. 
                    if kb == kmax
                        cstar_t_im1(k,i) = NaN;
                    else
                        % Interpolation done here using the weight
                        % functions calculated once above, outside the main
                        % time loop. 
                        cstar_t_im1(k,i) = weightBack(k,i)     * c_tm1(kb, i - 1) ...
                                       + (1 - weightBack(k,i)) * c_tm1(kb + 1, i - 1);
                    end
                    
                    if kf == kmax
                        cstar_t_ip1(k,i) = NaN;
                    else
                        cstar_t_ip1(k,i) = weightFront(k,i)       * c_tm1(kf,i + 1) ...
                                         + (1 - weightFront(k,i)) * c_tm1(kf + 1,i + 1);
                    end                
                end
            end
                        
            
        else % normal (slow) scheme
            cstar_t_ip1 = zeros(size(cstar_t));
            cstar_t_im1 = zeros(size(cstar_t));
            
            for i = 1:size(cstar_t,2)
                cstar_t_ip1(:,i) = interp1q(zc_ip1(:,i), c_tm1_ip1(:,i), zc(:,i)); % will have NaNs at bottom
                cstar_t_im1(:,i) = interp1q(zc_im1(:,i), c_tm1_im1(:,i), zc(:,i));
            end
        end
        
        I = find(isnan(cstar_t_ip1(1,:))==1); % remove nans in 1st line
        cstar_t_ip1(1,I) = c_z0;
        I = find(isnan(cstar_t_im1(1,:))==1); % remove nans in 1st line
        cstar_t_im1(1,I) = c_z0;
        
        I = find(cstar_t_ip1(1,:)==0); % remove zeros in 1st line
        cstar_t_ip1(1,I) = c_z0;
        I = find(cstar_t_im1(1,:)==0); % remove zeros in 1st line
        cstar_t_im1(1,I) = c_z0;
        
        trueDx_right = cstar_t_ip1 - cstar_t;
        trueDx_left = cstar_t - cstar_t_im1;
        
        
        
        % 2.1.1) True horizontal,but along sigma where we hit topography
        if HYBRID_DIFFUSION % keep NaN (can't diffuse throught bottom)

            D1 = zeros(size(trueDx_right));
            D2 = zeros(size(trueDx_right));
            D3 = 0; % Always zero...

            % record where we can't diffuse
            II = find(isnan(trueDx_right) == 1 | isnan(trueDx_left) == 1);            
            III = find(isnan(trueDx_right) == 1 & isnan(trueDx_left) == 1); 
            
            I = find(isnan(trueDx_left)==1);
            trueDx_left(I) = 0.0;
            I = find(isnan(trueDx_right)==1);
            trueDx_right(I) = 0.0;
            trueDx = H_t.*(Kx/(dx.^2)).*(trueDx_right - trueDx_left);
                       
            % Half&Half when can't true diffuse both side
            D1(II) = 0.5*Kx./(2*dx).^2.*(H_tm1_ip1(II) - H_tm1_im1(II)).*(c_tm1_ip1(II) - c_tm1_im1(II));
            D2(II) = 0.5*Kx.*H_tm1(II)./dx.^2.*(c_tm1_ip1(II) - 2*c_tm1(II) + c_tm1_im1(II));
            trueDx(II) = 0.5*trueDx(II);          
            
            % Only sigma when true impossible either side
            D1(III) = Kx./(2*dx).^2.*(H_tm1_ip1(III) - H_tm1_im1(III)).*(c_tm1_ip1(III) - c_tm1_im1(III));
            D2(III) = Kx.*H_tm1(III)./dx.^2.*(c_tm1_ip1(III) - 2*c_tm1(III) + c_tm1_im1(III));
            trueDx(III) = 0.0; 
            
             
            % One-side true / One-side sigma
% $$$             D1_right = (H_tm1_ip1 - H_tm1).*(c_tm1_ip1 - c_tm1);
% $$$             D1_left = (H_tm1 - H_tm1_im1).*(c_tm1 - c_tm1_im1);
% $$$             D2_right = c_tm1_ip1 - c_tm1;
% $$$             D2_left = c_tm1 - c_tm1_im1;
% $$$             trueDx = zeros(size(zc));
% $$$             D1 = zeros(size(zc));
% $$$             D2 = zeros(size(zc));
% $$$             
% $$$             % NaNs both sides (along sigma only)
% $$$             I = find(isnan(trueDx_right)==1 & isnan(trueDx_left)==1);
% $$$             trueDx(I) = 0.0;
% $$$             D1(I) = Kx./(2*dx).^2.*(D1_right(I) - D1_left(I));
% $$$             D2(I) = Kx.*H_tm1(I)./dx.^2.*(D2_right(I) - D2_left(I));
% $$$             
% $$$             % NaNs right side only (half and half)
% $$$             I = find(isnan(trueDx_right)==1 & isnan(trueDx_left)==0);
% $$$             trueDx(I) = H_t(I).*(Kx/(dx.^2)).*(0 - trueDx_left(I));
% $$$             D1(I) = Kx./(2*dx).^2.*(D1_right(I) - 0);
% $$$             D2(I) = Kx.*H_tm1(I)./dx.^2.*(D2_right(I) - 0);
% $$$             
% $$$             % NaNs left side only (half and half)
% $$$             I = find(isnan(trueDx_right)==0 & isnan(trueDx_left)==1);
% $$$             trueDx(I) = H_t(I).*(Kx/(dx.^2)).*(trueDx_right(I) - 0);
% $$$             D1(I) = Kx./(2*dx).^2.*(0 - D1_left(I));
                                                       % $$$             D2(I) = Kx.*H_tm1(I)./dx.^2.*(0 - D2_left(I));
% $$$             
% $$$             % no NaNs (true diffusion only)
% $$$             I = find(isnan(trueDx_right)==0 & isnan(trueDx_left)==0);
% $$$             trueDx(I) = H_t(I).*(Kx/(dx.^2)).*(trueDx_right(I) - trueDx_left(I));
% $$$             D1(I) = 0;
% $$$             D2(I) = 0;
% $$$            
% $$$             D3 = 0; % Always zero...
         
        % 2.1.2) True horizontal, but one-sided where we hit topography
        else % replace NaNs with zeros (no diffusion through bottom)
            D1 = 0;
            D2 = 0;
            D3 = 0;
            
            I = find(isnan(trueDx_left)==1);
            trueDx_left(I) = 0.0;
            I = find(isnan(trueDx_right)==1);
            trueDx_right(I) = 0.0;
            trueDx = H_t.*(Kx/(dx.^2)).*(trueDx_right - trueDx_left);
        end
        
    % 2.2) diffusion along z contours
    else
        trueDx = 0;
        
        D1 = (Kx./((2*dx).^2)).*(H_tm1_ip1 - H_tm1_im1).*(c_tm1_ip1 - c_tm1_im1);
        
        D2 = (Kx.*H_tm1./(dx.^2)).*(c_tm1_ip1 - 2*c_tm1 + c_tm1_im1);
        % not good at i = end, but doesnt matter since c_i=end = 170
        
        D3 = (Kz./((2*dsigma).^2)).*(1./H_tm1_sp1 - 1./H_tm1_sm1).*(c_tm1_sp1 - c_tm1_sm1);
        if sum(sum(diff(H_t),1)) < 1e-9 % which is zero
            D3(:,:) = 0; % Remove discontinuity because Ht constant and ddz(1/H)=Inf
        end
        
    end
    
    % vertical diffusion
    D4        = Kz./H_tm1./dsigma^2.*(c_tm1_sp1 - 2*c_tm1 + c_tm1_sm1);
    D4(1,:)   = Kz./H_tm1(1,:)./dsigma^2.*(c_tm1_sp1(1,:) - c_tm1(1,:)); % bottom (rigid)
    D4(end,:) = Kz./H_tm1(end,:)./dsigma^2.*(-c_tm1(end,:) + c_tm1_sm1(end,:)); % top (rigid)

    D = D1+D2+D3+D4+trueDx;
    
    % 3) Respiration term
    R = Rwc.*H_t; %in mmol/m2/s (units of [A] and [D])
    
    % 4) Sediment flux
    sedFlux = Rsed./Dz(end,:); % mmol/m3/s
    
    % ***** main equation ***** %
    % c(t+1) = c(t-1) + 2dt(X_adv(t) + X_diff(t-1))
    deltaC = 2*dt*(A + D + R)./H_t; % here in mmol/m3
    
    deltaC(end, :) = deltaC(end, :) + sedFlux*2*dt; %add sedimentflux
    
    % deltaC is the concentration variation at each grid point
    c_tp1 = c_tm1 + deltaC; % advance in time
    % ************************** %
    
    % Neumann condition
    c_tp1(:,1) = c_tp1(:,2);
    
    % Dirichlet conditions
    c_tp1(:,end) = c_xL;
    c_tp1(1, :) = c_z0;
    
    % Reinitialize concentration matrices
    %c_tm1 = c_t; % <----- keep Leap-frog stepping scheme
    c_tm1 = c_tp1; % <--- Gauss-Seidel (converge faster but may be instable)
    c_t = c_tp1;
    timeStep = timeStep + 1;
    
    
    if mod(timeStep,DISP_FREQ) == 0
        figure(1)
        clf
        set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 40 15])
        pcolor(X,Z,c_t)
        shading interp
        colorbar
        title(sprintf('timestep %d', timeStep))
        set(gca, 'ydir', 'reverse')
        pause(.1)
        
        % convergence treshold ( M = mean(c) )
        M = [M; sum(sum(c_t.*cellVolume))./sum(sum(cellVolume))];
        figure(2)
        plot(M)
        title('convergence')
        
        % convergence test (d[c]/dt<1e-10)
        if abs(M(end)-M(end-1))./(DISP_FREQ*dt) <= 1e-10
            no_cv = 0;
        end
        save O2.mat c_t X Z
    end
end

% ---------------------------- END MAIN LOOP ----------------------------- %

disp('converged!')
toc

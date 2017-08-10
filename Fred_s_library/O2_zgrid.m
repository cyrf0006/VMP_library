% This script was developed to duplicate and test the model proposed by
% Lefort et. al (2012) for hypoxia in the St. Lawrence System.
%
% This script solves numerically the 2D advection-diffusion equation for 
% dissolved oxygen over irregular bathymetry:
%
%   d[O]/dt + u d[O]/dx + w d[O]/dz 
%           = d/dx(Kx d[O]/dx) + d/dz(Kz d[O]/dz) - R               (1)
%
% where
% t:    time (s)
% x:    Horizontal axis, positive towards the right (m)
% z:    Vertical axis, positive downward (m)
% [O]:  Dissolved oxygen concentration (mmol/m^3)
% u:    Horizontal velocity (m/s)
% w:    Vertical velocity (m/s)
% Kx:   Horizontal diffusivity (m^2/s)
% Kz:   Vertical diffusivity (m^2/s)
% R:    Pelagic respiration (mmol/m^3/s)
% 
% The velocity field is obtained using the depth-averaged steady-state
% continuity equation assuming that the horizontal flow has no 
% depth-dependance and given the flow at the seaward boundary. This is:
%
%   d(Hu)/dx = 0
%
% which gives 
% 
%   u(x,z) = u0 H0 / H(x)                                                (2)
% 
% where H(x) is the thickness of the layer subject to hypoxia, H0 is the
% layer thickness at the seaward boundary and u0 the velocity at the 
% seaward boundary.
%
% Once u is determined from (2), the vertical velocity w is obtained from 
% the 2D continuity equation
%
%   du/dx + dw/dz = 0.                                               (3) 
%
% The model uses a C-grid defined as:
%
%  k=1 x~~~~~~~~~~~~~x~~~~~~~~~~~~~x~~~~~~~~~~~~~x~~~~~~~~~~~~~x~~~~~~ zeta(1) = 0
%    1 +------*------+------*------+------*------+------*------+------* z(1)
%    2 x-------------x-------------x-------------x-------------x        zeta(2)
%    2 +------*------+------*------+------*------+------*------+------* z(2)
%    3 x-------------x-------------x-------------x-------------x        zeta(3)
%    3 +------*------+------*------+------*------+------*------+------* z(3)
%    4 x-------------x-------------x-------------x-------------x
%  ... +------*------+------*------+------*------+------*------+------*
%  ... x-------------x-------------x-------------x-------------x
%  ... +------*------+------*------+------*------+------*------+------*
%  ... x-------------x-------------x-------------x-------------x
% kmax +------*------+------*------+------*------+------*------+------* z(kmax)
% kmax x-------------x-------------x-------------x-------------x        zeta(kmax)
%    i=1      1      2      2     ...    ...  imax-1  imax-1  imax  imax
%            x(1)          x(2)        x(...)       x(imax-1)       x(imax)
%   chi(1)        chi(2)       chi(...)     chi(imax-1)     chi(imax)
%
%
% where
%      *:       u
%      +:       [O]
%      x:       w
%      x(i):    Vector that defines the * grid
%      chi(i):  Vector that defines the + grid
%      zeta(k): Vector that defines the vertical position of the x grid points
%      z(k):    Vector that defines the vertical position of the * and + grid points
%
%
% The numerical scheme is fully-explicit using a mixed leapfrog-Euler 
% forward time-stepping scheme. The advective terms can be treated with 
% centered differences or an first-order upstream scheme.
%
% In order to mimic the model of Lefort et. al. results this model must be 
% ran until steady-state is obtained, in which case equation (1) reduces 
% to the equation used in Lefort et. al. 
%
% The topography must be carefully dealt with. The file is loaded from the
% file H_of_x.mat that must be in the curent folder. The topography is then
% adjusted to match the grid.
%
% Authors: 
% First version: Daniel Bourgault (April 2012)


%% SIMULATION PARAMETERS

% Output filename
%output_fname = 'O2_dx_1km.mat';
output_fname = 'O2_dx_10km.mat';

% Pause for graph
graph = true;       % Set to true to see updated graphs while the model runs
pause_in_sec = 0.1;

% Advecton scheme
% upstream = false means that the centered advection scheme is used.
upstream = false;

% Horizontal and vertical eddy diffusivity
Kx = 1600.0;
Kz = 1.d-4;

% Forcing at seaward boundary (m/s)
u0 = -0.003;

% Initial and boundary conditions (mmol/m^3)
O_bnd = 170.0;
O_ini = 170.0;

% Sediment oxygen demand 
SOD = 5500./(365*86400); % 5500 mmol/m^2/yr converted to 5500 mmol/m^2/s

% Pelagic respiration (mmol/m^3/s)
R = 0.0;

% Grid dimension and size (m)
dx = 10000.0;
dz = 5.0;
%dx = 2000.0;
%dz = 2.0;

% Domain maximum depth (Hmax in m), length (L in m) and level of the 
% rigid lid (z0 in m)
z0 = 150.0;
Hmax = 550;
L = 820000;
Hmax = Hmax - z0;

% Simulationm period (s)
T = 3.0*365*86400;

% Timestep. Automatically determined given dx and Kx in order to respect 
% the stability criterion.
dt_x = (1/8)*(dx.^2)/Kx;
dt_z = (1/8)*(dz.^2)/Kz;
dt = round(min(dt_x,dt_z));

% Maximum number of iteration
nmax = T/dt;  % Maximum number of iteration


%% Make the grid
chi =  [0:dx:L];
zeta = [0:dz:Hmax];
imax = length(chi);
kmax = length(zeta);
x = [chi(1)+dx/2:dx:chi(imax)+dx/2];
z = [zeta(1)+dz/2:dz:zeta(kmax)+dz/2];


% Bathymetry
% The file H_of_x xontains the bathymetry at 1 km resolution. 
load H_of_x.mat;
dx0 = x0(2)-x0(1);
H0 = H0 - z0;
  
% Filter the topography if dx is larger than the 1 km data resolution.
if dx > dx0
    win = round(dx/dx0)
    if mod(win,2) == 0
        win = win + 1;
    end
    H0 = sgolayfilt(H0,0,win);  
end

% Interpolate the bathymetry at the grid size.
H0 = interp1(x0,H0,chi);


% Initialization
% flag_u is a matrix with 1 and 0 for wet and dry cells for u 
% flag_s is a matrix with 1 and 0 for wet and dry cells for scaler [O] 
flag_u(1:imax,1:kmax) = 0.0;
flag_s(1:imax,1:kmax) = 0.0;

u(1:imax,1:kmax) = 0.0;
w(1:imax,1:kmax) = 0.0;

% O     oxygen at timestep n
% O_1   oxygen at timestep n-1
% O_2   oxygen at timestep n-2
O(1:imax,1:kmax)   = O_ini;
O_1(1:imax,1:kmax) = O_ini;
O_2(1:imax,1:kmax) = O_ini;

% Boundary conditions as in Lefort et al.
O(imax,:)   = O_bnd;
O(:,1)      = O_bnd;
O_1(imax,:) = O_bnd;
O_1(:,1)    = O_bnd;
O_2(imax,:) = O_bnd;
O_2(:,1)    = O_bnd;


% Determine the depth at w grid. Variable Hw.
for i = 1:imax
    for k = 2:kmax
        if H0(i) < zeta(k)
            if abs(zeta(k)-H0(i)) <= abs(zeta(k-1)-H0(i))
                Hw(i) = zeta(k);
            else
                Hw(i) = zeta(k-1);
            end
            break
        else
            Hw(i) = H0(i);
        end
    end
end

% Determine the depth at u grid. Variable Hu. 
for i = 1:imax-1
    if Hw(i+1) == Hw(i)
        Hu(i) = Hw(i);
    elseif Hw(i) < Hw(i+1);
        Hu(i) = Hw(i);
    else
        Hu(i) = Hw(i+1);
    end
end
Hu(imax) = Hw(imax);

% Create the flag matrix to identify wet and dry cells given Hu and Hw.
for i = 1:imax
    for k = 1:kmax
        if z(k) < Hu(i)
            flag_u(i,k) = 1.0;
        end
        if z(k) < Hw(i)
            flag_s(i,k) = 1.0;
        end
    end
end

% Fin dht index of wet and dry cells.
index_wet_s = find(flag_s == 1.0);
index_wet_u = find(flag_u == 1.0);
index_dry_s = find(flag_s == 0.0);

% Put crazy values at dry cell to make sure that these cells are not
% manipulated otherwise it would cause obvious problems. 
O(index_dry_s)   = -999999;
O_1(index_dry_s) = -999999;
O_2(index_dry_s) = -999999;


% Determine the bottom cells of each water column. 
for i = 1:imax
    for k = 1:kmax-1
        if flag_u(i,k+1) == 0.0;
            kb_u(i) = k;
            break;
        end
    end
end
for i = 1:imax
    for k = 1:kmax-1
        if flag_s(i,k+1) == 0.0;
            kb_s(i) = k;
            break;
        end
    end
end

% Compute the horizontal velocity field u(x,z)
for i = 1:imax
    for k = 1:kb_u(i)
        u(i,k) = u0.*Hu(end)/Hu(i);
    end
end

% Compute the vertical velocity field w(x,z)
for i = 2:imax
    for k = kb_s(i):-1:2
        w(i,k) = dz*(u(i,k) - u(i-1,k))/dx + w(i,k+1);
    end
end

% Interpolation of u and w on the same grid. Only
% needed fot visualizarion of the results. Not used
% in the main algorithm.
[X,Z] = meshgrid(chi,z);
U = interp2(x,z,u',X,Z);
W = interp2(chi,zeta,w',X,Z,'linear');

% Main time loop
for n = 1:nmax
    for i = 2:imax-1
        for k = 2:kb_s(i)
            
            % terms for horizontal advection
            um = (u(i,k)+u(i-1,k))/2;
            dOdxb = flag_s(i-1,k)*(O_1(i,k)-O_1(i-1,k))/dx;
            dOdxf = flag_s(i+1,k)*(O_1(i+1,k)-O_1(i,k))/dx;
            
            % Terms for vertical advection
            wm = (w(i,k+1)+w(i,k))/2;
            dOdzb = flag_s(i,k-1)*(O_1(i,k)-O_1(i,k-1))/dz;
            dOdzf = flag_s(i,k+1)*(O_1(i,k+1)-O_1(i,k))/dz;

            % Advection 
            if upstream
                advx = (um+abs(um)/2)*dOdxb + (um-abs(um)/2)*dOdxf;
                advz = (wm+abs(wm)/2)*dOdzb + (wm-abs(wm)/2)*dOdzf;
            else
                advx = um*(dOdxf+dOdxb)/2.0;
                advz = wm*(dOdzf+dOdzb)/2.0;
            end

            % Horizontal diffusion
            flux_x_b = flag_s(i-1,k)*Kx*(O_2(i,k)-O_2(i-1,k))/dx;
            flux_x_f = flag_s(i+1,k)*Kx*(O_2(i+1,k)-O_2(i,k))/dx;
            diffx = (flux_x_f - flux_x_b)/dx;
            
            % Vertical diffusion
            flux_z_b = flag_s(i,k-1)*Kz*(O_2(i,k)-O_2(i,k-1))/dz;
            if flag_s(i,k+1) == 0.0
                flux_z_f = -SOD;  % This is where SOD is imposed
            else
                flux_z_f = flag_s(i,k+1)*Kz*(O_2(i,k+1)-O_2(i,k))/dz;
            end
            diffz = (flux_z_f - flux_z_b)/dz;
            
            % Solve for oxygen
            O(i,k) = 2*dt*(-advx - advz + diffx + diffz - R) + O_2(i,k);
            
        end
    end
    
    % Upstream boundary condition dO/dx = 0
    O(1,:) = O(2,:);
    
    % Update
    O_2 = O;
    O_1 = O;
    
    % Fin the min and max
    minO(n) = min(O(index_wet_s));
    maxO(n) = max(O(index_wet_s));
    
    % graph the results if requested
    if graph
        n_out = 100;
        if mod(n,n_out) == 0
            figure(1)
            clf
            %set(gcf,'Position',[100 1000 1000 300]);
            Ograph = O;
            Ograph(index_dry_s)   = NaN;
            colormap(flipud(jet));
            pcolor(chi,z,Ograph');
            shading('interp');
            set(gca,'ydir','reverse');
            xlabel('Distance (km)')
            ylabel('Depth (m)')
            title(sprintf('timestep %d', n))
            hold on;
            %quiver(X,Z,U,W,0.5,'ko');
            viztopo(chi,zeta,Hw);
            colorbar;
            pause(pause_in_sec)
        end
    end
    
end
O(index_dry_s)   = NaN;

save(output_fname,'x','z','chi','zeta','dx','dz','dt','Kx','Kz','Hu','Hw','u','w','minO','maxO','O');


figure(1);
pcolor(x/1000,z,flag_u');
hold;
plot(x/1000,H0);
plot(x/1000,Hw,'k');
shading('flat');
set(gca,'ydir','reverse');
xlabel('x (km)')
ylabel('z (m)')

figure(2);
subplot(2,1,1)
pcolor(x/1000,zeta,u');
shading('interp');
colorbar;
set(gca,'ydir','reverse','tickdir','out','XMinorTick','on','YMinorTick','on');
title('u (m/s)');
xlabel('Distance (km)');
ylabel('Depth (m)');

subplot(2,1,2)
pcolor(chi/1000,z,w');
shading('interp');
colorbar;
set(gca,'ydir','reverse','tickdir','out','XMinorTick','on','YMinorTick','on');
title('w (m/s)');
xlabel('x (km)');
ylabel('z (m)');


figure(3)
title ('Minimum and maximu O value');
plot(minO);
hold;
plot(maxO);


figure(4)
colormap(flipud(jet));
pcolor(chi./1000,z,O');
shading('interp');
set(gca,'ydir','reverse');
xlabel('x (km)')
ylabel('z (m)')
hold;
viztopo(chi/1000,zeta,Hw);
colorbar;


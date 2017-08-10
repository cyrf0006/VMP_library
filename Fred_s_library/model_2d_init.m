% model namelist
clear

dt = 600; %sec
dz = 1; %m
dx = 100; %m

load('H_of_x.mat');

x = x0(1):dx:x0(end);
H = round(interp1(x0, H0, x));
z0 = -150; % top of domain
zmax = -525;
z = z0:-dz:zmax;

Kx = 1.6e3; %m/s2
Kz = 1e-4; %m/s2
Rwc = -19.6; %mmol/m3/yr
Rsed = -3540; %mmol/m2/yr
u0 = -0.53; %m/s at x=L
init_cond = 170; %mmol/m3;
c_xL = 170; %mmol/m3 (b.c. at x = L)
c_z0 = 170; %mmol/m3 (b.c. at z = z0);


% Compute U W on the grid
% d(UH)/dx = 0 -> UH = cst
UH = u0*H(end);
U = UH./H; %U(x) = cst/H(x)
W = -U.*gradient(H)./dx;


% initializing the grid (c is the tracer)
c_grid = zeros(length(z), length(x));
for i = 1:length(x)
    c_grid(1:H(i), i) = 170;
    u_grid(1:H(i), i) = U(i);
    w_grid(1:H(i), i) = W(i);
end



figure(1) % initial condition
imagesc(x, z, c_grid)
set(gca, 'ydir', 'normal')



% Compute U W on the grid
% d(UH)/dx = 0 -> UH = cst
UH = u0*H(end);
U = UH./H; %U(x) = cst/H(x)

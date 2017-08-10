function REG = K_gotm(tprofile, p_t, gotmk, dt, DT, L)

% tprofile = initial Tprofile
% p_t = pressure profile, must corespond to gotmk
% gotmK = time-depth dependent K , ex: 'nuh_clim_highrelax2.dat'
% dt = timestep of the computation
% DT = timestep for the output (correspond to gotmk time dim)
% L = duration of the calculation
%
% ex: REG = K_gotm(T_climato_04(:,2), T_climato_04(:,1), 'nuh_clim_highrelax2.dat', 15, 3600*24*30, 3600*24*30*8);


% Time/depth-dependent K (from GOTM)
A = load(gotmk);
K_GOTM = flipdim(A(:,2:301)',1);
K_month = A(:,1);

figure(1) %idea of K
contourf(K_month, p_t, K_GOTM)
set(gca, 'ydir', 'reverse')
title('K as modeled by GOTM')
ylabel('depth')
xlabel('month')

% binned variables (Pb, Tb, Sb, DENSb)
Zmin = min(p_t);
Zmax = max(p_t);
dz = Zmax/length(p_t);

C = dt/dz^2; % coefficient

Pb = Zmin:dz:Zmax;
for i = 1:length(Pb)
    I = find(p_t >= Pb(i)-dz/2 &  p_t <= Pb(i)+dz/2); %find data both side of the bin center
    Tb(i) = mean(tprofile(I));
end


% 1st K
count = 1;
K = K_GOTM(1:length(p_t)-1,count)';
K1 = [K 0];  %boundary conditions
K2 = [0 K]; 

T = Tb;
Temp(count, :) = T; %

figure(2)
plot(T, Pb);
axis([-1 6 Zmin Zmax])
set(gca, 'ydir', 'reverse')
title(sprintf('timestep %d', count))
pause(0.25)  

hold on

count = count +1;

for t = 1:L/dt %loop on time

%    keyboard
    T_mk = [0 T(1:length(T)-1)]; % T_k-1
    T_pk = [T(2:length(T)) 0]; % T_k+1

    T_av = T + C.*(K1.*(T_pk - T) - K2.*(T - T_mk));
    
    if count == t*dt/DT % plot every DT (s)

        figure(2)
        plot(T_av, Pb);
        axis([-1 6 Zmin Zmax])
        set(gca, 'ydir', 'reverse')
        title(sprintf('timestep %d', t))
        pause(0.25)  
        
        Temp(count, :) = T;
%        Sal(count, :) = S;
        count = count +1;
        
        %New K
        K = K_GOTM(1:length(p_t)-1,count)';
        K1 = [K 0];  %boundary conditions
        K2 = [0 K]; 
        
    end
    
     T = T_av;
     %S = S_av;
end

hold off
% Regression of the core min T

for i = 1:length(Temp(:,1))
    monthly_min(i) = min(Temp(i,:));
end


figure(3)
clf
plot(1:length(Temp(:,1)), monthly_min, '.');
hold on

X = 1:length(Temp(:,1));
REG = polyfit(X, monthly_min, 1);

Y = REG(1).*X+REG(2);
plot(X, Y);



%Real regression
 temp = load('T_climato_04.dat');
 TT(1) = min(temp(:,2));
 temp = load('T_climato_05.dat');
 TT(2) = min(temp(:,2));
 temp = load('T_climato_06.dat');
 TT(3) = min(temp(:,2));
 temp = load('T_climato_07.dat');
 TT(4) = min(temp(:,2));
 temp = load('T_climato_08.dat');
 TT(5) = min(temp(:,2));
 temp = load('T_climato_09.dat');
 TT(6) = min(temp(:,2));
 temp = load('T_climato_10.dat');
 TT(7) = min(temp(:,2));
 temp = load('T_climato_11.dat');
 TT(8) = min(temp(:,2));

plot(1:length(Temp(:,1)), TT, 'k.')

REG2 = polyfit(X, TT, 1);
Y = REG2(1).*X+REG2(2);
plot(X, Y, 'k');
title({sprintf('K = %d', k_cst)})

hold off
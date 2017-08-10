%function Cd_fit(infile)
clear

data = load('/home/cyrf0006/WINDEX/data_processing/BMix_study/meanU.mat');
kappa = 0.41; 
zmax = 5;
zmin = 2;
zUb = 5;

load resVel.mat

for i = 1:length(data.reg_tide)
    

    global Z U Ub
    
    Z = data.meanU(:,1);
    U = data.meanU(:,2);
    U = data.U_tide(:,i);
    %U = nanmean(abs(data.u_tide),2);
% $$$     U = abs(res_along);
% $$$     Z = zhab;
    
    I = find(Z<=zmax & Z>=zmin);
    Z = Z(I);
    U = U(I);
 
    [Y, I] = min(abs(Z-zUb));
    Ub = U(I);
    Cd_0 = 1e-3;
    z0_0 = 0.0001;    

    % 1st guess
    a0 = [Cd_0 z0_0];
    U0 = sqrt(a0(1))./kappa.*Ub.*log(Z./a0(2));

    
    figure(1)
    clf
    plot(U, Z)
    set(gca, 'ydir','reverse')
    hold on
    %plot(U0, Z, 'g--')
    
    
    % Minimization
    [asol vals] = fminsearch ( 'fitFun3', a0);
    
    % best fit
    U1 = sqrt(asol(1))./kappa.*Ub.*log(Z./asol(2));
    plot(U1, Z, 'r--')
    hold off
    title(sprintf('%d', i))
     
    %disp(['Cd = ' num2str(asol(1)) ' | Z_0 = ' num2str(asol(2)) ' | Ub = ' num2str(asol(3))])
    %disp(['Cd = ' num2str(asol(1)) ' | Z_0 = ' num2str(asol(2))])
    
    Cd(i) = asol(1);
    Z0(i) = asol(2);
    errors(i) = vals;
    
end


I = find(errors>0);       
B = errors(I);
A = data.reg_tide(I);
clf
plot(A, B);

[Y, J] = min(B);

J = I(J);
U = data.U_tide(:,J);
Z = data.meanU(:,1);

[Y, I] = min(abs(Z-zUb));
Ub = U(I);

I = find(Z<=zmax & Z>=zmin);
Z = Z(I);
U = U(I);  

bestCd = Cd(J);
bestZ0 = Z0(J);


figure(1)
clf
plot(U, Z)
set(gca, 'ydir','reverse')
hold on
U1 = sqrt(Cd(J))./kappa.*Ub.*log(Z./Z0(J));
plot(U1, Z, '--r')
hold off
disp(['Best combination at time2 = ' num2str(data.reg_tide(J)) ' for profile 0-' num2str(zmax) 'm'])
title(['Cd = ' num2str(Cd(J)) ' | Z_0 = ' num2str(Z0(J)) ' at time2 = ' num2str(data.reg_tide(J)) ])
set(gca, 'ydir','normal')
xlabel('U (m/s)')
ylabel('hab (m)')
legend('observed vel. profile', 'best fit')
errors(J)


myZ = data.meanU(:,1);
myU = data.U_tide;
[Y, I] = min(abs(myZ-zUb));
myUb = myU(I,:);

eps = nan(size(myU));
for i = 1:length(myUb)
    eps(:,i) = bestCd.^(3/2)*myUb(i).^(3)./(kappa*myZ);
end

    
figure(2)
clf
imagesc(data.reg_tide, myZ, myU)
set(gca, 'ydir','normal')
title('Total current velocities')
colorbar
xlabel('Time to hightide (h)')
ylabel('hab (m)')


figure(3)
clf
imagesc(data.reg_tide, myZ, log10(eps))
set(gca, 'ydir','normal')
title('Predicted \epsilon (W/kg)')
colorbar
xlabel('Time to hightide (h)')
ylabel('hab (m)')



figure(4)
clf
plot(data.reg_tide, nanmean(eps,1))
ylabel('Predicted \epsilon (W/kg)')
xlabel('Time to hightide (h)')
set(gca, 'yscale', 'log')

reg_tide = data.reg_tide;
eps3 = nanmean(eps,1);
save predicted_eps.mat reg_tide eps3

save predicted_eps_contours2.mat reg_tide myZ eps
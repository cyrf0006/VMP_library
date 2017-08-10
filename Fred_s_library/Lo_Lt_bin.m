function Lo_Lt_bin(prof_files, eps_files)
    
% Lo_Lt_bin('Tprofiles_riki', 'Eprofiles_riki')

% few params:
g = 9.81; %m/s2
zmin = 0; % top of the fisrt bin
zmax = 50;
zbin = 5;
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';

% profiles*.mat
fid = fopen(prof_files);
C = textscan(fid, '%s', 'delimiter', '\n');
list = char(C{1});
no_profiles = size(list, 1);

% overturns
fid = fopen(eps_files);
C = textscan(fid, '%s', 'delimiter', '\n');
list_eps = char(C{1});
if no_profiles ~= size(list_eps, 1)
    disp('Error, no. profiles is not equal to no. eps_profiles')
    exit
end

count=1;
for i = 1:no_profiles
    
    disp(sprintf('profile %d', i));
    
    % ---- Work on eps profiles ---- %
    fname_in = list_eps(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    load(fname_in)
     
    
     % average of 2 epsilon profiles 
    if length(p_eps1) == length(p_eps2) % no problem, easy
        p_eps = p_eps1;
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        EPS2(:,2) = eps2;
    elseif min(p_eps1)<min(p_eps2) % p1 drive, nomatter size of p2!
        p_eps = p_eps1;        
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        ind = find(p_eps2 == p_eps1(1));
        EPS2(ind:length(p_eps2),2) = eps2;
    else % p2 drive, nomatter size of p1!
        p_eps = p_eps2;
        p_N = p_eps1;
        EPS2=nan(length(p_eps2),2);
        EPS2(:,2) = eps2;
        ind = find(p_eps1 == p_eps2(1));
        EPS2(ind:length(p_eps1),1) = eps1;
    end
    
    % "selection" average
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
   
    EPS = nanmean(EPS2,2); %MEAN EPSILON
                        
    % if only nans
    if sum(~isnan(EPS))==0;
        continue
    end
    
    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8);

    % uses home despike
    EPS = Y;
    
    
    % mean N2 
    N2=nan(length(p_eps),1);
    if length(p_N)==length(p_eps);
        N2(:) = N.^2;
    else
        N2(ind:length(p_N))=N.^2;
    end
    
    % bin profiles 1-m
    for k = 1:length(P_bin)
        I = find(p_eps >= P_bin(k)-zbin/2 & p_eps <= P_bin(k)+zbin/2);

        if ~isempty(I)==1
            eps_bin(k) = nanmean(EPS(I));
            N2(k) = nanmean(N2(I));
        else
            eps_bin(k) = NaN;
            N2(k)=NaN;
        end
    end
    % --------------------------------- %
    
    
    
    
    % ---- Work on T profiles ---- %
    fname_in = list(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    load(fname_in)
    DENS = sw_dens(SBS, SBT, P);

    I = find(isnan(eps_bin)==1);
    eps_bin(I)=[];
    P_bin(I)=[];
    
    for j = 1:length(P_bin)
        I = find(P >= P_bin(j)-zbin & P <= P_bin(j)+zbin);

        rho_raw = DENS(I); 
        rho_sort = sort(rho_raw);
        rho0 = nanmean(rho_raw);
        
        rms = abs(rho_raw - rho_sort);
        L_t = nanmean(rms);
        
        % Ozmidov scale
% $$$         % N = sqrt(N2(j)); %original
% $$$         p = polyfit(P(I),rho_raw,1);
        rho0 = mean(rho_raw);
        p = polyfit(P(I),rho_sort,1);
        rho0 = mean(rho_sort);
        
        rho0 = mean(rho_sort);
        N = sqrt((g/rho0)*p(1));
        
        
        L_o = sqrt(eps_bin(j)./N.^3);
        
        % save thorpe scale 
        thorpescale(count) = L_t; 
        % save Ozmidov scale
        ozmidovscale(count) = L_o;
        buo(count) = N;
        eps(count) = eps_bin(j);
        count = count + 1;

    end                
    % -------------------------- % 
    
    
end


% N vs epsilon
figure(1)
clf
I = find(~isnan(buo)==1);
A = buo(I);
B = eps(I);
I = find(~isnan(B)==1);  
B = B(I);
A = A(I);
[A, I] = sort(A);
B = B(I);
P=polyfit(log10(A), log10(B), 1);
C = log10(A).*P(1)+P(2);

plot(log10(A), log10(B), '.k')
hold on
plot(log10(A), C, 'r')

xlabel('N (s^{-1})')
ylabel('\epsilon (W kg^{-1})')



% Lo vs Lt
figure(2)
clf
I = find(~isnan(ozmidovscale)==1);
A = ozmidovscale(I);
B = thorpescale(I);
I = find(~isnan(B)==1);  
B = B(I);
A = A(I);
[A, I] = sort(A);
B = B(I);
P=polyfit(log10(A), log10(B), 1);
C = log10(A).*P(1)+P(2);

plot(log10(A), log10(B), '.k')
hold on
plot(log10(A), C, 'r')

xlabel('Lo')
ylabel('Lt')

% ---- test with bootstrap ----- %
dN = 0.1;
reg_N = -4.8:dN:-2;
N2_vector = buo.^2;
eps_vector = eps;

reg_eps = nan(round(length(eps_vector)/4), length(reg_N)); % /2 its just approx
for i = 1: length(reg_N)
    I = find(log10(N2_vector)>=reg_N(i)-dN/2 & log10(N2_vector)<=reg_N(i)+dN/2);
    reg_eps(1:length(I), i) = log10(eps_vector(I));
end

figure(4)
clf
imagesc(reg_N, 1:length(eps_vector)/4, sort(reg_eps, 1))
title('\epsilon distribution')
%ylim([0 120])
ylabel('\epsilon (W kg^{-1})')
colorbar

nboot = 500;
for i = 1:length(reg_N)
    I = find(~isnan(reg_eps(:,i)));
    samp = reg_eps(I,i);
    NN = length(samp);
    clear eps_boot_b
    % create random sampling
    for b = 1:nboot
        r = rand(NN,1);
        r = ceil(r*NN/1);
        eps_boot_b(b) = nanmean(samp(r));
    end

    % mean of random sampling
    eps_boot_dot(i) = nanmean(eps_boot_b);

    % compute error (EFRON & GONG 83)
    eps_error(i) = sqrt(sum((diff([eps_boot_b eps_boot_dot(i)],1, 2)).^2, 2)./(nboot-1));
end

figure(3)
clf
errorbar(reg_N, eps_boot_dot, eps_error, '.k')
xlabel('log10(N^2) (s^{-2})')
ylabel('log10(\epsilon) (W kg^{-1})')
P=polyfit(reg_N, eps_boot_dot, 1);
C = reg_N.*P(1)+P(2);
hold on
plot(reg_N, C, 'k')



keyboard
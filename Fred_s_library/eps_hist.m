function eps_hist(epsfile, prof, nboot)
    
    

mat_eps  = load(epsfile);
zvec = 1:4:301; 

[Y,zbin] = min(abs(prof-zvec));

% -- bootstrap -- %
N = size(mat_eps,2);

epsvec = mat_eps(zbin,:);

% create random sampling
for b = 1:nboot
    %r = rand(N,1);
    r = 0.5 + (26.4999-0.5).*rand(N,1);
    r = round(r);
    %r = ceil(r*N/1);
    %keyboard
    nu =  mean(log(epsvec(r)), 2);
    sig = std(log(epsvec(r)),0,2);
    eps_boot_b2(b) = exp(nu+sig.^2/2);
    eps_boot_b(b) =   mean(epsvec(r), 2);

end

% mean of random sampling
eps_boot_dot = nanmean(eps_boot_b,2);
figure(1)
hist(eps_boot_b)
figure(2)
hist(eps_boot_b2)

% compute error (EFRON & GONG 83)
eps_error = sqrt(sum((diff([eps_boot_b eps_boot_dot],1, 2)).^2, 2)./(nboot-1));

% compute 95% confidence interval
%eps_2p5 = 

keyboard
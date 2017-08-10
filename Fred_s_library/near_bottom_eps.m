function near_bottom_eps(eps_files, Pavg)

% function timedepth_eps(no_profile, epsMin, epsMax, varargin)
% was previously isopyc_tidecycle.m
% 
% ex: near_bottom_eps('eps_files', 10)
%
% "ls -1 eps_profile*.mat | sed 's/\.mat//' > eps_files"
%
% 
%


zmin = 0;
zmax = 350;
zbin = 1;
Pbin = [zmin+zbin/2:zbin:zmax-zbin/2]';
hab = Pbin;
%Pavg = 20;  % last Pavg meters of the water column average

g = 9.81;

% load eps_files file
fid = fopen(eps_files);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

no_profile = size(epsfiles, 1); %number of *.P files 

% find time min
fname = epsfiles(1, :);
I = find(fname==' '); %remove white space  
fname(I) = [];
load(fname)
T1 = mtime_eps(1);

% find time max
fname = epsfiles(end, :); 
I = find(fname==' '); %remove white space  
fname(I) = [];
load(fname)
T2 = mtime_eps(1);

% raw matrix to fill
epsMat = nan(length(Pbin), no_profile);
NMat = nan(length(Pbin), no_profile);

fid = fopen('./Matrices.mat');
if fid == -1 % doesnt exist

    for i = 1:no_profile
        
        disp(sprintf('profile %d', i))
        fname = epsfiles(i, :);
        I = find(fname==' ');   
        fname(I) = [];    
        load(fname)
        
        epsilon = nanmean([eps1', eps2'], 2);
        
        % Bin profiles
        epsbin = nan(length(Pbin),1);
        Nbin = nan(length(Pbin),1);
        
        for j = 1:length(Pbin)
            I = find(p_eps1 >= (Pbin(j) - zbin/2) & p_eps1 <= (Pbin(j) + zbin/2));
            if ~isempty(I) == 1
                epsbin(j) = nanmean(epsilon(I));
                Nbin(j) = nanmean(N(I));
            end
        end
        
        % rename variables and reference relative to hab
        epsilon = epsbin;
        N = Nbin;
        I = find(~isnan(N)==1); % we take N on purpose
        epsilon(1:I(end)) = flipud(epsilon(1:I(end))); % now hab
        N(1:I(end)) = flipud(N(1:I(end))); % now hab
        NMat(:, i) = N;
        epsMat(:, i) = epsilon;
        proftime(i) = mtime_eps(1);    
    end 
    save Matrices.mat epsMat NMat proftime

else
    load Matrices
end



timevec = proftime;
% clean N2 field 
% (linear itp in time to remove nans)
for i = 1:size(NMat,1)
    raw_vec = NMat(i,:);
    I = find(~isnan(raw_vec)==1);
    if ~isempty(I)==1 & length(I)>1
        raw_vec = raw_vec(I);
        vec_itp = interp1(timevec(I), raw_vec, timevec(I(1):I(end)));
        NMat(i,I(1):I(end)) = vec_itp;
    end
end


% Reduce vertical size and rename for plotting (copy-paste)
I = find(hab<=Pavg);
zvec = hab(I);
eps_mat = epsMat(I,:);
N_mat = NMat(I,:);




I = find(hab<=Pavg);

% mean profile
E = epsMat(I, :);
m = nanmean(log(E), 2);
s2 = nanvar(log(E), 2);
meanEps = exp(m+s2/2);

% mean value
E = E(:);
N = length(E);
nboot = 1000;
% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    eps_boot_b(b) = nanmean(E(r));
end

% mean of random sampling
eps_boot_dot = nanmean(eps_boot_b);
eps_sort = sort(eps_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);
eps_2p5 = eps_sort(:,CI_2p5);
eps_97p5 = eps_sort(:,CI_97p5);


m = nanmean(log(E));
s2 = nanvar(log(E));
meanEps = exp(m+s2/2);



keyboard


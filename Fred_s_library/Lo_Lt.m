function Lo_Lt(prof_files, ovt_files)
    
% Lo_Lt('Tprofiles_riki', 'fineovt_riki')

% few params:
g = 9.81; %m/s2
    
% profiles*.mat
fid = fopen(prof_files);
C = textscan(fid, '%s', 'delimiter', '\n');
list = char(C{1});
no_profiles = size(list, 1);

% overturns
fid = fopen(ovt_files);
C = textscan(fid, '%s', 'delimiter', '\n');
list_ovt = char(C{1});
if no_profiles ~= size(list_ovt, 1)
    disp('Error, no. profiles is not equal to no. overturns')
    exit
end

count=1;
for i = 1:no_profiles

    disp(sprintf('profile %d', i));
    % ---- Work on VMP profiles ---- %
    fname_in = list(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    load(fname_in)
     
    % --------------------------------- %
    
    % -------- Work on overturns ------ %
    fname_ovt = list_ovt(i, :); % remove white space
    I = find(fname_ovt==' ');
    fname_ovt(I)=[];
    
    % get overturns fotnhis profile
    command = sprintf('cut -d" " -f2 %s > /tmp/tmp', fname_ovt);
    system(command);
    ind1 = load('/tmp/tmp');

    command = sprintf('cut -d" " -f3 %s > /tmp/tmp', fname_ovt);
    system(command);
    ind2 = load('/tmp/tmp');

% $$$     command = sprintf('cut -d" " -f4 %s > /tmp/tmp', fname_ovt);
% $$$     system(command);
% $$$     z1 = load('/tmp/tmp');
% $$$     
% $$$     command = sprintf('cut -d" " -f5 %s > /tmp/tmp', fname_ovt);
% $$$     system(command);
% $$$     z2 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f6 %s > /tmp/tmp', fname_ovt);
    system(command);
    L_t = load('/tmp/tmp');

    
    if isempty(L_t)==1
        continue
    end
    
  
    % else, loop on overturns

    % few param.
    ax = ang2acc(pitch); ay = ang2acc(roll);
    DENS = sw_dens(SBS, SBT, P);
    [nu, mu] = viscosity(SBS, SBT, DENS);    
    fs = round(fs);
    FS = round(FS);

    no_ovt = length(L_t);
    for j = 1:no_ovt
        
        if ind2(j)-ind1(j) < 75 % overturn at least 10cm
            continue
        end
        
        II = [ind1(j):ind2(j)]'; %for P
        I = find(p > P(ind1(j)) & p < P(ind2(j))); % P -> p
        Wm = nanmean(W(II)); % mean falling speed for this ovt 
        visco = nanmean(nu(II)); % mean viscosity for this ovt   
        
        % preanalysis (cleaning, flagging, etc.)
        sh1 = shear_preanalysis(shear1(I), w(I), pitch(I), roll(I), p(I), 5);
        sh2 = shear_preanalysis(shear2(I), w(I), pitch(I), roll(I), p(I), 5);
        
        % spectral integration
        eps1 = spec_int(sh1, [ax(I) ay(I) az(I)], Wm, visco, fs);
        eps2 = spec_int(sh2, [ax(I) ay(I) az(I)], Wm, visco, fs);
        
        % mean epsilon
        if eps1 > 10*eps2
            eps2 = nan;
        elseif eps2 > 10*eps1
            eps1 = nan;
        end
        
        epsilon = nanmean([eps1 eps2]);
        
        % Ozmidov scale
        rho = sort(DENS(II));
        rho0 = nanmean(rho);
        pp = polyfit(P(II), rho, 1);
        N = sqrt((g/rho0)*pp(1));
        L_o = sqrt(epsilon./N.^3);
        
        % save thorpe scale 
        thorpescale(count) = L_t(j); 
        % save Ozmidov scale
        ozmidovscale(count) = L_o;
        buo(count) = N;
        eps(count) = epsilon;
        count = count + 1;

    end                
    % -------------------------- % 
    
    
end


keyboard
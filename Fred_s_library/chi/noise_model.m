function noise_model(profiles, zbin, probe)

% function noise_model(chifiles)
%
% 
% usage ex: noise_model('Tprofiles_all', 4, 1)
%   to be run in ~/WINDEX/data_processing/batchelor/chi
% or
% ex: noise_model('hit_bottom_2011_2', 4, 1)
%   to be run in ~/WINDEX/data_processing/batchelor/chi_hit_bottom
% F. Cyr - april 2012

       
fid = fopen(profiles);
C = textscan(fid, '%s', 'delimiter', '\n');
pro_files = char(C{1});

no_profiles = size(pro_files, 1);

clear Wm
for count = 1:no_profiles

    fname_in = pro_files(count, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    
    load(fname_in);

    
    d=sprintf('TREATMENT OF PROFILE %d ...', count);
    disp(d)
    
    % Angle and viscosity
    DENS = sw_dens(SBS, SBT, P);
    [nu, mu] = viscosity(SBS, SBT, DENS);

    % round frequencies
    fs = round(fs);
    FS = round(FS);

    dp = p(2)-p(1);
    
    if probe == 1
        dtdz = gradient(t1, dp);
    elseif probe == 2
        dtdz = gradient(t2, dp);
    else
        disp('Wrong probe number!')
        break
    end
    
    
    % ---------------------------- FFT calculation --------------------- %
    
    % Compute good indexes
    % interp from slow freq. (64Hz) to fast freq. (512Hz) (zs = zshear)
    z = p;
    %zs=interp1([1:length(z)],z,[1:length(dtdz)]'*(FS/fs)); % <------------------ERROR!!

% $$$     if isempty(zs)==1 % no shear data to process
% $$$         disp('skip bin')
% $$$         continue;
% $$$     end

    z1 = ceil(min(abs(z)));  % top of profile
    z2 = floor(max(abs(z)));  % bottom of profile
    zfft = z1+zbin/2:zbin:z2-zbin/2; %reg. z prof. for FFT

    if isempty(zfft) % no shear data to process
        disp('profile empty... skip!')
        continue;
    end % if

    I=find(z>=z1 & z<=z2); % we keep only full bins
    dtdz=dtdz(I);
    z = z(I); 
 
    for i = 1:length(zfft)

        % compute indexes on which FFT will be applied 
        I = find(z>zfft(i)-zbin/2 & z<zfft(i)+zbin/2); % I = Imicrosctruc
        Ifine = find(P>zfft(i)-zbin/2 & P<zfft(i)+zbin/2); % Ifinesctruc

        if (isempty(I)==1)
            continue
        end
        
        II = find(isnan(dtdz(I))==1);        
        if ~isempty(II) == 1
            disp('skip bin')
            continue % skip if at least one NaN in bin (full bin only)
        end

        
        % Compute FFT 
        % make sure that freq. vector is always the same
        if count == 1 & i == 1
            [ps0, f0] = pwelch(dtdz(I), [], [], [], 512); % Raw Power Spectrum
            Freq = f0;
            Nfft = length(f0)*2-1;
            noise = ps0;
        end

        [ps0, f0] = pwelch(dtdz(I), [], [], Nfft, 512); % Raw Power Spectrum        
                                                      
        %[size(f0), f0(1), f0(end)]
        
        % Conversion Freq -> Wavenumber
        Wm(count, i) = mean(w(Ifine)); % mean falling speed for this bin         
                                       %        k0 = f0/Wm;
        
        I = find(ps0 < noise);
        if ~isempty(I) == 1
            noise(I) = ps0(I);
        end
        
    end

 
end


keyboard


% ---- Apply a smoothing filter ---- %
a = 1;
window = 10; % no. of pts
b = ones(1,window);
b = b/length(b);
noisefilt = filtfiltm(b,a,noise', 2);
% ------------------------------------ %

Wm = Wm(:);
I = find(Wm<0.5);
Wm(I) = [];
fnoise = Freq;
knoise = fnoise/mean(Wm);



outfile = sprintf('noise_%dm_veryall.mat', zbin)
save(outfile, 'noise', 'noisefilt', 'fnoise', 'knoise');

